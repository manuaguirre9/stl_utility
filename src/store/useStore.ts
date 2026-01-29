import { create } from 'zustand';
import * as THREE from 'three';
import { v4 as uuidv4 } from 'uuid';
import { subdivideSelectedFaces, getContiguousIslands } from '../utils/meshUtils';
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
    historyPreviewId: string | null;

    // Re-edit State
    reEditParams: any | null;
    preReEditIndex: number | null;

    // History Actions
    undo: () => void;
    redo: () => void;
    jumpToHistory: (index: number) => void;
    recordHistory: (label: string, action: HistoryAction) => void;
    deleteHistoryEntry: (id: string | string[]) => void;
    selectedHistoryIds: string[];
    setSelectedHistoryIds: (ids: string[]) => void;
    setHistoryPreviewId: (id: string | null) => void;
    reEditHistoryItem: (index: number) => void;
    recalculateHistoryItem: (id: string, newParams: any) => void;

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
    historyPreviewId: null,
    reEditParams: null,
    preReEditIndex: null,
    selectedHistoryIds: ['initial'],

    setSelectedHistoryIds: (ids) => {
        set({ selectedHistoryIds: ids });

        // Auto-populate re-edit params if a tool-related action is selected
        if (ids.length === 1) {
            const entry = get().history.find(h => h.id === ids[0]);
            if (entry && entry.id !== 'initial') {
                const action = entry.action;
                let toolMode: TransformMode = null;
                let params: any = null;

                if (action.type === 'SUBDIVIDE') {
                    toolMode = 'subdivide'; params = { steps: action.steps };
                } else if (action.type === 'TEXTURIZE_KNURLING') {
                    toolMode = 'texturize'; params = { type: 'knurling', ...action.params };
                } else if (action.type === 'TEXTURIZE_HONEYCOMB') {
                    toolMode = 'texturize'; params = { type: 'honeycomb', ...action.params };
                } else if (action.type === 'TEXTURIZE_DECIMATE') {
                    toolMode = 'texturize'; params = { type: 'decimate', ...action.params };
                }

                if (toolMode) {
                    set({ transformMode: toolMode, reEditParams: params });
                }
            }
        }
    },

    setHistoryPreviewId: (id) => set({ historyPreviewId: id }),

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
            selectedHistoryIds: [entry.id],
            reEditParams: null
        });
    },

    deleteHistoryEntry: (ids) => {
        const { history } = get();
        const idList = Array.isArray(ids) ? ids : [ids];
        const idsToRemove = new Set(idList);

        if (idsToRemove.has('initial')) idsToRemove.delete('initial');
        if (idsToRemove.size === 0) return;

        const newHistoryActions = history.filter(h => !idsToRemove.has(h.id)).map(h => ({ label: h.label, action: h.action }));

        let currentModels: ModelData[] = [];
        let currentSelection: Record<string, number[]> = {};
        const rebuiltHistory: HistoryEntry[] = [];

        newHistoryActions.forEach((item) => {
            const action = item.action;
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
                id: uuidv4(),
                label: item.label,
                action: action,
                models: cloneModelsSnapshot(currentModels),
                smartSelection: JSON.parse(JSON.stringify(currentSelection)),
                timestamp: Date.now()
            });
        });

        const lastEntry = rebuiltHistory[rebuiltHistory.length - 1];
        set({
            history: rebuiltHistory,
            historyIndex: rebuiltHistory.length - 1,
            models: cloneModelsSnapshot(lastEntry.models),
            smartSelection: JSON.parse(JSON.stringify(lastEntry.smartSelection)),
            selectedHistoryIds: [lastEntry.id]
        });
    },

    undo: () => {
        const { historyIndex, history, deleteHistoryEntry } = get();
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
                selectedHistoryIds: [state.id]
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
                selectedHistoryIds: [state.id],
                preReEditIndex: null // Clear jump-back memory on manual navigation
            });
        }
    },

    reEditHistoryItem: (index) => {
        const { history, historyIndex } = get();
        const entry = history[index];
        if (!entry || entry.id === 'initial') return;

        const prevState = history[index - 1];
        const action = entry.action;
        let toolMode: TransformMode = null;
        let params: any = null;
        let selection: number[] = [];
        let modelId: string | null = null;

        if (action.type === 'SUBDIVIDE') {
            toolMode = 'subdivide'; params = { steps: action.steps }; selection = action.selection; modelId = action.modelId;
        } else if (action.type === 'TEXTURIZE_KNURLING') {
            toolMode = 'texturize'; params = { type: 'knurling', ...action.params }; selection = action.selection; modelId = action.modelId;
        } else if (action.type === 'TEXTURIZE_HONEYCOMB') {
            toolMode = 'texturize'; params = { type: 'honeycomb', ...action.params }; selection = action.selection; modelId = action.modelId;
        } else if (action.type === 'TEXTURIZE_DECIMATE') {
            toolMode = 'texturize'; params = { type: 'decimate', ...action.params }; selection = action.selection; modelId = action.modelId;
        }

        if (toolMode && modelId) {
            set({
                models: cloneModelsSnapshot(prevState.models),
                smartSelection: { [modelId]: selection },
                historyIndex: index, // Scrubber to the right of the item
                selectedId: modelId,
                transformMode: toolMode,
                reEditParams: params,
                selectedHistoryIds: [entry.id],
                preReEditIndex: historyIndex // Save where we were
            });
        }
    },

    recalculateHistoryItem: (id, newParams) => {
        const { history, historyIndex, preReEditIndex, selectedId } = get();
        const editIndex = history.findIndex(h => h.id === id);
        if (editIndex <= 0) return;

        // Determine which index to restore to. 
        // If we have a saved preReEditIndex, we go back there. 
        // Otherwise we stay at the current historyIndex.
        const originalScrubberIndex = preReEditIndex !== null ? preReEditIndex : historyIndex;

        // 1. Prepare ALL actions for replay
        const actionsToReplay = history.map((h) => {
            if (h.id === id) {
                const oldAction = h.action;
                let updatedAction = { ...oldAction };
                if (oldAction.type === 'SUBDIVIDE') {
                    (updatedAction as any).steps = newParams.steps;
                } else if (oldAction.type.startsWith('TEXTURIZE')) {
                    (updatedAction as any).params = { ...newParams };
                }
                return { label: h.label, action: updatedAction };
            }
            return { label: h.label, action: h.action };
        });

        // 2. Re-apply everything from scratch
        let currentModels: ModelData[] = [];
        let currentSelection: Record<string, number[]> = {};
        const rebuiltHistory: HistoryEntry[] = [];

        actionsToReplay.forEach((item) => {
            const action = item.action as any;
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
                id: uuidv4(),
                label: item.label,
                action: action,
                models: cloneModelsSnapshot(currentModels),
                smartSelection: JSON.parse(JSON.stringify(currentSelection)),
                timestamp: Date.now()
            });
        });

        // 3. Restore User Context
        const finalRestoreIndex = Math.min(originalScrubberIndex, rebuiltHistory.length - 1);
        const restoredState = rebuiltHistory[finalRestoreIndex];

        set({
            history: rebuiltHistory,
            historyIndex: finalRestoreIndex,
            models: cloneModelsSnapshot(restoredState.models),
            smartSelection: JSON.parse(JSON.stringify(restoredState.smartSelection)),
            selectedHistoryIds: [rebuiltHistory[editIndex].id],
            selectedId: selectedId,
            preReEditIndex: null // Clear the memory
        });
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

    setTransformMode: (mode) => set({ transformMode: mode, reEditParams: null }),
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

        const islands = getContiguousIslands(model.bufferGeometry, selection)
            .sort((a, b) => Math.min(...a) - Math.min(...b));

        let currentGeometry = model.bufferGeometry;
        const processedOriginals = new Set<number>();

        set({ smartSelection: { ...get().smartSelection, [modelId]: [] }, showMesh: true });

        islands.forEach((island, idx) => {
            const translatedIsland = island.map(origIdx => {
                let shift = 0;
                processedOriginals.forEach(o => { if (o < origIdx) shift++; });
                return origIdx - shift;
            });

            const nextGeometry = subdivideSelectedFaces(currentGeometry, translatedIsland, steps);
            island.forEach(i => processedOriginals.add(i));
            currentGeometry = nextGeometry;

            set((state) => ({
                models: state.models.map(m => m.id === modelId ? { ...m, bufferGeometry: currentGeometry, meshVersion: m.meshVersion + 1 } : m)
            }));

            get().recordHistory(`Subdivide Surface ${idx + 1}`, {
                type: 'SUBDIVIDE', modelId, steps, selection: translatedIsland
            });
        });
    },

    applyTexturize: (modelId, params) => {
        const model = get().models.find(m => m.id === modelId);
        if (!model) return;
        const selection = get().smartSelection[modelId] || [];
        if (selection.length === 0) return;

        const islands = getContiguousIslands(model.bufferGeometry, selection)
            .sort((a, b) => Math.min(...a) - Math.min(...b));

        let currentGeometry = model.bufferGeometry;
        const processedOriginals = new Set<number>();

        set({ smartSelection: { ...get().smartSelection, [modelId]: [] }, showMesh: true });

        islands.forEach((island, idx) => {
            const translatedIsland = island.map(origIdx => {
                let shift = 0;
                processedOriginals.forEach(o => { if (o < origIdx) shift++; });
                return origIdx - shift;
            });

            let nextGeometry = currentGeometry;
            if (params.type === 'knurling') {
                nextGeometry = applyKnurling(currentGeometry, translatedIsland, params);
            } else if (params.type === 'honeycomb') {
                nextGeometry = applyHoneycomb(currentGeometry, translatedIsland, params);
            } else if (params.type === 'decimate') {
                nextGeometry = applyDecimate(currentGeometry, translatedIsland, params);
            }

            island.forEach(i => processedOriginals.add(i));
            currentGeometry = nextGeometry;

            set((state) => ({
                models: state.models.map(m => m.id === modelId ? { ...m, bufferGeometry: currentGeometry, meshVersion: m.meshVersion + 1 } : m)
            }));

            get().recordHistory(`Texturize ${params.type} - Surface ${idx + 1}`,
                params.type === 'knurling'
                    ? { type: 'TEXTURIZE_KNURLING', modelId, params: { ...params }, selection: translatedIsland }
                    : params.type === 'honeycomb'
                        ? { type: 'TEXTURIZE_HONEYCOMB', modelId, params: { ...params }, selection: translatedIsland }
                        : { type: 'TEXTURIZE_DECIMATE', modelId, params: { ...params }, selection: translatedIsland }
            );
        });
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
