import { create } from 'zustand';
import * as THREE from 'three';
import { v4 as uuidv4 } from 'uuid';
import { subdivideSelectedFaces } from '../utils/meshUtils';
import { applyKnurling, applyHoneycomb, applyDecimate } from '../utils/texturizeUtils';
import type { KnurlPattern } from '../utils/texturizeUtils';

export interface ModelData {
    id: string;
    name: string;
    bufferGeometry: THREE.BufferGeometry;
    color: string;
    visible: boolean;
    position: [number, number, number];
    rotation: [number, number, number];
    scale: [number, number, number];
    meshVersion: number;
}

export type TransformMode = 'smart' | 'subdivide' | 'texturize' | null;

export type HistoryAction =
    | { type: 'IMPORT', model: ModelData }
    | { type: 'REMOVE', modelId: string }
    | { type: 'UPDATE', modelId: string, data: Partial<ModelData> }
    | { type: 'SUBDIVIDE', modelId: string, steps: number, selection: number[] }
    | { type: 'TEXTURIZE_KNURLING', modelId: string, params: { pitch: number, depth: number, angle: number, pattern: KnurlPattern }, selection: number[] }
    | { type: 'TEXTURIZE_HONEYCOMB', modelId: string, params: { cellSize: number, wallThickness: number, depth: number }, selection: number[] }
    | { type: 'TEXTURIZE_DECIMATE', modelId: string, params: { reduction: number }, selection: number[] }
    | { type: 'INITIAL' };

interface HistoryEntry {
    id: string;
    label: string;
    action: HistoryAction;
    // We still keep snapshots for performance, but they will be rebuilt on deletion
    models: ModelData[];
    smartSelection: Record<string, number[]>;
    timestamp: number;
}

interface AppState {
    models: ModelData[];
    selectedId: string | null;
    transformMode: TransformMode;
    showMesh: boolean;
    showEdges: boolean;
    showClassification: boolean;
    classificationAngle: number;
    smartSelection: Record<string, number[]>;

    // History State
    history: HistoryEntry[];
    historyIndex: number;
    selectedHistoryId: string | null;

    // History Actions
    undo: () => void;
    redo: () => void;
    jumpToHistory: (index: number) => void;
    recordHistory: (label: string, action: HistoryAction) => void;
    deleteHistoryEntry: (id: string) => void;
    setSelectedHistoryId: (id: string | null) => void;

    // Standard Actions
    addModel: (model: ModelData) => void;
    removeModel: (id: string) => void;
    selectModel: (id: string | null) => void;
    updateModel: (id: string, data: Partial<ModelData>) => void;
    setTransformMode: (mode: TransformMode) => void;
    toggleShowMesh: () => void;
    toggleShowEdges: () => void;
    toggleClassification: () => void;
    setClassificationAngle: (angle: number) => void;
    setSmartSelection: (modelId: string, indices: number[]) => void;
    addToSmartSelection: (modelId: string, indices: number[]) => void;
    removeFromSmartSelection: (modelId: string, indices: number[]) => void;
    clearSmartSelection: () => void;
    subdivideSelection: (modelId: string, steps: number) => void;
    applyTexturize: (modelId: string, params:
        | { type: 'knurling', pitch: number, depth: number, angle: number, pattern: KnurlPattern }
        | { type: 'honeycomb', cellSize: number, wallThickness: number, depth: number }
        | { type: 'decimate', reduction: number }
    ) => void;
    selectAllFaces: (modelId: string) => void;
}

const cloneModelsSnapshot = (models: ModelData[]): ModelData[] => {
    return models.map(m => ({
        ...m,
        bufferGeometry: m.bufferGeometry.clone()
    }));
};

export const useStore = create<AppState>((set, get) => ({
    models: [],
    selectedId: null,
    transformMode: null,
    showMesh: true,
    showEdges: true,
    showClassification: false,
    classificationAngle: 20,
    smartSelection: {},

    history: [{
        id: 'initial',
        label: 'Start Project',
        action: { type: 'INITIAL' },
        models: [],
        smartSelection: {},
        timestamp: Date.now()
    }],
    historyIndex: 0,
    selectedHistoryId: 'initial',

    setSelectedHistoryId: (id) => set({ selectedHistoryId: id }),

    recordHistory: (label, action) => {
        const { models, smartSelection, history, historyIndex } = get();
        const newHistory = history.slice(0, historyIndex + 1);

        const entry: HistoryEntry = {
            id: uuidv4(),
            label,
            action,
            models: cloneModelsSnapshot(models),
            smartSelection: JSON.parse(JSON.stringify(smartSelection)),
            timestamp: Date.now()
        };

        set({
            history: [...newHistory, entry],
            historyIndex: newHistory.length,
            selectedHistoryId: entry.id
        });
    },

    deleteHistoryEntry: (id) => {
        const { history } = get();
        const entryToDelete = history.find(h => h.id === id);
        if (!entryToDelete || entryToDelete.action.type === 'INITIAL') return;

        // 1. Remove from sequence
        const newHistoryActions = history.filter(h => h.id !== id).map(h => ({ label: h.label, action: h.action }));

        // 2. Re-calculate everything from scratch
        let currentModels: ModelData[] = [];
        let currentSelection: Record<string, number[]> = {};
        const rebuiltHistory: HistoryEntry[] = [];

        newHistoryActions.forEach((item) => {
            const action = item.action;

            // Apply action to current state
            if (action.type === 'INITIAL') {
                currentModels = [];
                currentSelection = {};
            } else if (action.type === 'IMPORT') {
                currentModels = [...currentModels, { ...action.model, bufferGeometry: action.model.bufferGeometry.clone() }];
            } else if (action.type === 'REMOVE') {
                currentModels = currentModels.filter(m => m.id !== action.modelId);
            } else if (action.type === 'UPDATE') {
                currentModels = currentModels.map(m => m.id === action.modelId ? { ...m, ...action.data } : m);
            } else if (action.type === 'SUBDIVIDE') {
                const model = currentModels.find(m => m.id === action.modelId);
                if (model) {
                    model.bufferGeometry = subdivideSelectedFaces(model.bufferGeometry, action.selection, action.steps);
                    model.meshVersion++;
                }
            } else if (action.type === 'TEXTURIZE_KNURLING') {
                const model = currentModels.find(m => m.id === action.modelId);
                if (model) {
                    model.bufferGeometry = applyKnurling(model.bufferGeometry, action.selection, action.params);
                    model.meshVersion++;
                }
            } else if (action.type === 'TEXTURIZE_HONEYCOMB') {
                const model = currentModels.find(m => m.id === action.modelId);
                if (model) {
                    model.bufferGeometry = applyHoneycomb(model.bufferGeometry, action.selection, action.params);
                    model.meshVersion++;
                }
            } else if (action.type === 'TEXTURIZE_DECIMATE') {
                const model = currentModels.find(m => m.id === action.modelId);
                if (model) {
                    model.bufferGeometry = applyDecimate(model.bufferGeometry, action.selection, action.params);
                    model.meshVersion++;
                }
            }

            rebuiltHistory.push({
                id: uuidv4(), // Generate new IDs for the rebuilt track to avoid conflicts
                label: item.label,
                action: action,
                models: cloneModelsSnapshot(currentModels),
                smartSelection: JSON.parse(JSON.stringify(currentSelection)),
                timestamp: Date.now()
            });
        });

        // 3. Update state to the last successful state in the new history
        const lastEntry = rebuiltHistory[rebuiltHistory.length - 1];
        set({
            history: rebuiltHistory,
            historyIndex: rebuiltHistory.length - 1,
            models: cloneModelsSnapshot(lastEntry.models),
            smartSelection: JSON.parse(JSON.stringify(lastEntry.smartSelection)),
            selectedHistoryId: lastEntry.id
        });
    },

    undo: () => {
        const { historyIndex, history, deleteHistoryEntry } = get();
        // Don't delete the initial state
        if (historyIndex > 0) {
            const entryToDelete = history[historyIndex];
            deleteHistoryEntry(entryToDelete.id);
        }
    },

    redo: () => {
        const { historyIndex, history } = get();
        if (historyIndex < history.length - 1) {
            const nextIndex = historyIndex + 1;
            const state = history[nextIndex];
            set({
                models: cloneModelsSnapshot(state.models),
                smartSelection: JSON.parse(JSON.stringify(state.smartSelection)),
                historyIndex: nextIndex,
                selectedHistoryId: state.id
            });
        }
    },

    jumpToHistory: (index) => {
        const { history } = get();
        if (index >= 0 && index < history.length) {
            const state = history[index];
            set({
                models: cloneModelsSnapshot(state.models),
                smartSelection: JSON.parse(JSON.stringify(state.smartSelection)),
                historyIndex: index,
                selectedHistoryId: state.id
            });
        }
    },

    addModel: (model) => {
        set((state) => ({ models: [...state.models, model] }));
        get().recordHistory(`Import ${model.name}`, { type: 'IMPORT', model });
    },

    removeModel: (id) => {
        const model = get().models.find(m => m.id === id);
        set((state) => ({
            models: state.models.filter((m) => m.id !== id),
            selectedId: state.selectedId === id ? null : state.selectedId
        }));
        get().recordHistory(`Remove ${model?.name || 'Object'}`, { type: 'REMOVE', modelId: id });
    },

    selectModel: (id) => set({ selectedId: id }),

    updateModel: (id, data) => {
        set((state) => ({
            models: state.models.map((m) => (m.id === id ? { ...m, ...data } : m))
        }));
        if (data.color || data.name) {
            get().recordHistory(`Edit ${data.name ? 'Name' : 'Color'}`, { type: 'UPDATE', modelId: id, data });
        }
    },

    setTransformMode: (mode) => set({ transformMode: mode }),
    toggleShowMesh: () => set((state) => ({ showMesh: !state.showMesh })),
    toggleShowEdges: () => set((state) => ({ showEdges: !state.showEdges })),
    toggleClassification: () => set((state) => ({ showClassification: !state.showClassification })),
    setClassificationAngle: (angle) => set({ classificationAngle: angle }),
    setSmartSelection: (modelId: string, indices: number[]) => set((state) => ({
        smartSelection: { ...state.smartSelection, [modelId]: indices }
    })),
    addToSmartSelection: (modelId: string, indices: number[]) => set((state) => {
        const current = state.smartSelection[modelId] || [];
        const merged = Array.from(new Set([...current, ...indices]));
        return { smartSelection: { ...state.smartSelection, [modelId]: merged } };
    }),
    removeFromSmartSelection: (modelId: string, indices: number[]) => set((state) => {
        const current = state.smartSelection[modelId] || [];
        const toRemove = new Set(indices);
        const filtered = current.filter(idx => !toRemove.has(idx));
        return { smartSelection: { ...state.smartSelection, [modelId]: filtered } };
    }),
    clearSmartSelection: () => set({ smartSelection: {} }),

    subdivideSelection: (modelId, steps) => {
        const model = get().models.find(m => m.id === modelId);
        if (!model) return;
        const selection = get().smartSelection[modelId] || [];
        if (selection.length === 0) return;

        const newGeometry = subdivideSelectedFaces(model.bufferGeometry, selection, steps);

        set((state) => ({
            models: state.models.map(m => m.id === modelId ? { ...m, bufferGeometry: newGeometry, meshVersion: m.meshVersion + 1 } : m),
            smartSelection: { ...state.smartSelection, [modelId]: [] },
            showMesh: true
        }));

        get().recordHistory(`Subdivide (${steps} steps)`, {
            type: 'SUBDIVIDE',
            modelId,
            steps,
            selection: [...selection]
        });
    },

    applyTexturize: (modelId, params) => {
        const model = get().models.find(m => m.id === modelId);
        if (!model) return;
        const selection = get().smartSelection[modelId] || [];
        if (selection.length === 0) return;

        let newGeometry = model.bufferGeometry;
        if (params.type === 'knurling') {
            newGeometry = applyKnurling(model.bufferGeometry, selection, {
                pitch: params.pitch,
                depth: params.depth,
                angle: params.angle,
                pattern: params.pattern
            });
        } else if (params.type === 'honeycomb') {
            newGeometry = applyHoneycomb(model.bufferGeometry, selection, {
                cellSize: params.cellSize,
                wallThickness: params.wallThickness,
                depth: params.depth
            });
        } else if (params.type === 'decimate') {
            newGeometry = applyDecimate(model.bufferGeometry, selection, {
                reduction: params.reduction
            });
        }

        set((state) => ({
            models: state.models.map(m => m.id === modelId ? { ...m, bufferGeometry: newGeometry, meshVersion: m.meshVersion + 1 } : m),
            smartSelection: { ...state.smartSelection, [modelId]: [] },
            showMesh: true
        }));

        get().recordHistory(`Texturize: ${params.type}`,
            params.type === 'knurling'
                ? { type: 'TEXTURIZE_KNURLING', modelId, params: { pitch: params.pitch, depth: params.depth, angle: params.angle, pattern: params.pattern }, selection: [...selection] }
                : params.type === 'honeycomb'
                    ? { type: 'TEXTURIZE_HONEYCOMB', modelId, params: { cellSize: params.cellSize, wallThickness: params.wallThickness, depth: params.depth }, selection: [...selection] }
                    : { type: 'TEXTURIZE_DECIMATE', modelId, params: { reduction: params.reduction }, selection: [...selection] }
        );
    },

    selectAllFaces: (modelId) => {
        const model = get().models.find(m => m.id === modelId);
        if (!model) return;
        const faceCount = model.bufferGeometry.attributes.position.count / 3;
        const allIndices = Array.from({ length: faceCount }, (_, i) => i);
        set((state) => ({
            smartSelection: { ...state.smartSelection, [modelId]: allIndices }
        }));
    },
}));
