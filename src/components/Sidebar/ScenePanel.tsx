import React, { useState, useRef } from 'react';
import { useStore } from '../../store/useStore';
import type { KnurlPattern } from '../../utils/texturizeUtils';
import { Eye, EyeOff, Box, Trash2, Scissors, Grid3X3, X, RotateCcw, Zap, Plus } from 'lucide-react';
import { Tooltip } from '../UI/Tooltip';
import { v4 as uuidv4 } from 'uuid';
import { STLLoader } from 'three/examples/jsm/loaders/STLLoader.js';
import { fitCylinderToSelection } from '../../utils/meshUtils';

interface ScenePanelProps {
    onClose?: () => void;
}

export const ScenePanel: React.FC<ScenePanelProps> = ({ onClose }) => {
    const models = useStore((state) => state.models);
    const selectedId = useStore((state) => state.selectedId);
    const selectModel = useStore((state) => state.selectModel);
    const updateModel = useStore((state) => state.updateModel);
    const removeModel = useStore((state) => state.removeModel);
    const addModel = useStore((state) => state.addModel);
    const subdivideSelection = useStore((state) => state.subdivideSelection);
    const applyTexturize = useStore((state) => state.applyTexturize);
    const selectAllFaces = useStore((state) => state.selectAllFaces);
    const smartSelection = useStore((state) => state.smartSelection);
    const transformMode = useStore((state) => state.transformMode);
    const reEditParams = useStore((state) => state.reEditParams);
    const selectedHistoryIds = useStore((state) => state.selectedHistoryIds);
    const recalculateHistoryItem = useStore((state) => state.recalculateHistoryItem);
    const textureType = useStore((state) => state.textureType);
    const setTextureType = useStore((state) => state.setTextureType);
    const isProcessing = useStore((state) => state.isProcessing);
    const applyManifoldRepair = useStore((state) => state.applyManifoldRepair);

    const fileInputRef = useRef<HTMLInputElement>(null);
    const [subdivideSteps, setSubdivideSteps] = useState(1);

    // Texturize states
    const [pitch, setPitch] = useState(2.0);
    const [referencePitch, setReferencePitch] = useState(2.0);
    const [depth, setDepth] = useState(0.5);
    const [angle, setAngle] = useState(0);
    const [cellSize, setCellSize] = useState(5.0);
    const [referenceCellSize, setReferenceCellSize] = useState(5.0);
    const [wallThickness, setWallThickness] = useState(0.5);
    const [knurlPattern, setKnurlPattern] = useState<KnurlPattern>('diamond');
    const [reductionRatio, setReductionRatio] = useState(0.5);
    const [direction, setDirection] = useState<'inward' | 'outward'>('inward');
    const [holeFillEnabled, setHoleFillEnabled] = useState(true);
    const [fuzzyThickness, setFuzzyThickness] = useState(0.2);
    const [pointDistance, setPointDistance] = useState(0.2);
    const [repairEdgeMult] = useState(50);
    const [repairHoleMult] = useState(200);

    // Sync with re-edit params
    React.useEffect(() => {
        if (reEditParams) {
            if (transformMode === 'subdivide') {
                if (reEditParams.steps !== undefined) setSubdivideSteps(reEditParams.steps);
            } else if (transformMode === 'texturize') {
                if (reEditParams.type) setTextureType(reEditParams.type);
                if (reEditParams.pitch !== undefined) { setPitch(reEditParams.pitch); setReferencePitch(reEditParams.pitch); }
                if (reEditParams.depth !== undefined) setDepth(reEditParams.depth);
                if (reEditParams.angle !== undefined) setAngle(reEditParams.angle);
                if (reEditParams.pattern) setKnurlPattern(reEditParams.pattern);
                if (reEditParams.cellSize !== undefined) { setCellSize(reEditParams.cellSize); setReferenceCellSize(reEditParams.cellSize); }
                if (reEditParams.wallThickness !== undefined) setWallThickness(reEditParams.wallThickness);
                if (reEditParams.reduction !== undefined) setReductionRatio(reEditParams.reduction);
                if (reEditParams.direction) setDirection(reEditParams.direction);
                if (reEditParams.holeFillEnabled !== undefined) setHoleFillEnabled(reEditParams.holeFillEnabled);
                if (reEditParams.thickness !== undefined) setFuzzyThickness(reEditParams.thickness);
                if (reEditParams.pointDistance !== undefined) setPointDistance(reEditParams.pointDistance);
            }
        }
    }, [reEditParams, transformMode]);

    const isHistoryReEdit = selectedHistoryIds.length === 1 && selectedHistoryIds[0] !== 'initial' && reEditParams !== null;
    const selectedModel = models.find((m) => m.id === selectedId);

    // Stitching State
    const [stitchStats, setStitchStats] = React.useState<{ openEdges: number, nonManifoldEdges: number, avgEdgeLength: number } | null>(null);
    const [analyzedHoles, setAnalyzedHoles] = React.useState<any[] | null>(null);
    const [rawBoundaries, setRawBoundaries] = React.useState<Float32Array | null>(null);
    const [isRepairing, setIsRepairing] = React.useState(false);
    const setBoundaryEdges = useStore((state) => state.setBoundaryEdges);

    // Auto-analyze when entering stitching mode
    React.useEffect(() => {
        if (transformMode === 'stitching' && selectedModel) {
            import('../../utils/meshUtils').then(({ analyzeMesh }) => {
                const stats = analyzeMesh(selectedModel.bufferGeometry);
                setStitchStats(stats);
                setBoundaryEdges(selectedModel.id, stats.openEdgePositions, stats.nonManifoldEdgePositions);
            });
            import('../../utils/manifoldRepairService').then(({ analyzeRepairableHoles }) => {
                const results = analyzeRepairableHoles(selectedModel.bufferGeometry);
                setAnalyzedHoles(results.holes);
                setRawBoundaries(results.rawBoundaries);
            });
        } else {
            setStitchStats(null);
            setAnalyzedHoles(null);
            setRawBoundaries(null);
            if (selectedId) {
                setBoundaryEdges(selectedId, null, null);
                useStore.getState().setFillPreviewEdges(selectedId, null, null);
            }
        }
    }, [transformMode, selectedModel?.meshVersion, selectedId]);

    // Live update hole preview based on sliders
    React.useEffect(() => {
        if (!selectedId || !stitchStats || !analyzedHoles) return;
        const avgEdge = stitchStats.avgEdgeLength;
        const maxTriSize = avgEdge * repairEdgeMult;
        const maxHoleDiameter = avgEdge * repairHoleMult;

        const fillable: number[] = [];
        const unfillable: number[] = [];

        analyzedHoles.forEach(hole => {
            if (hole.diameter <= maxHoleDiameter && hole.maxEdge <= maxTriSize) {
                for (let i = 0; i < hole.positions.length; i++) fillable.push(hole.positions[i]);
            } else {
                for (let i = 0; i < hole.positions.length; i++) unfillable.push(hole.positions[i]);
            }
        });

        if (rawBoundaries) {
            for (let i = 0; i < rawBoundaries.length; i++) {
                unfillable.push(rawBoundaries[i]);
            }
        }

        useStore.getState().setFillPreviewEdges(
            selectedId,
            fillable.length > 0 ? new Float32Array(fillable) : null,
            unfillable.length > 0 ? new Float32Array(unfillable) : null
        );
    }, [analyzedHoles, rawBoundaries, repairEdgeMult, repairHoleMult, stitchStats, selectedId]);

    const handleRepair = async () => {
        if (!selectedId) return;
        setIsRepairing(true);
        try {
            await applyManifoldRepair(selectedId, {
                edgeMultiplier: 50,
                holeDiameterMultiplier: 200
            });
        } finally {
            setIsRepairing(false);
        }
    };

    const selection = selectedId ? smartSelection[selectedId] || [] : [];

    const cylinderFit = React.useMemo(() => {
        if (!selectedModel || selection.length === 0 || transformMode !== 'texturize') return null;
        const fit = fitCylinderToSelection(selectedModel.bufferGeometry, selection);
        return (fit && !fit.isPlanar && fit.avgR) ? fit : null;
    }, [selectedId, selection.length, transformMode]);

    const circumference = cylinderFit ? 2 * Math.PI * cylinderFit.avgR! : 0;

    // In history re-edit mode, we might not have a face selection in the current live mesh, 
    // but the button should still show if we are browsing history.
    const selectionCount = selection.length;
    const hasSelection = selectionCount > 0;

    const showSubdividePanel = transformMode === 'subdivide' && (hasSelection || isHistoryReEdit);
    const showTexturizePanel = transformMode === 'texturize' && (hasSelection || isHistoryReEdit);


    const handleApply = () => {
        if (isHistoryReEdit) {
            // Recalculate existing operation
            const params = transformMode === 'subdivide'
                ? { steps: subdivideSteps }
                : { type: textureType, pitch, depth, angle, pattern: knurlPattern, cellSize, wallThickness, reduction: reductionRatio, direction, holeFillEnabled, thickness: fuzzyThickness, pointDistance };
            recalculateHistoryItem(selectedHistoryIds[0], params);
        } else if (selectedModel) {
            // Apply new operation
            if (transformMode === 'subdivide') {
                subdivideSelection(selectedModel.id, subdivideSteps);
            } else if (transformMode === 'texturize') {
                if (textureType === 'knurling') {
                    applyTexturize(selectedModel.id, { type: 'knurling', pitch, depth, angle, pattern: knurlPattern, holeFillEnabled });
                } else if (textureType === 'honeycomb') {
                    applyTexturize(selectedModel.id, { type: 'honeycomb', cellSize, wallThickness, depth, angle, direction, holeFillEnabled });
                } else if (textureType === 'fuzzy') {
                    applyTexturize(selectedModel.id, { type: 'fuzzy', thickness: fuzzyThickness, pointDistance, holeFillEnabled });
                } else if (textureType === 'decimate') {
                    applyTexturize(selectedModel.id, { type: 'decimate', reduction: reductionRatio });
                }
            }
        }
    };

    const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
        const file = event.target.files?.[0];
        if (!file) return;

        const reader = new FileReader();
        reader.onload = (e) => {
            const contents = e.target?.result as ArrayBuffer;
            const loader = new STLLoader();
            const geometry = loader.parse(contents);

            geometry.computeBoundingSphere();
            geometry.computeBoundingBox();

            const sphere = geometry.boundingSphere;
            const box = geometry.boundingBox;

            if (!sphere || !box || isNaN(sphere.radius) || !isFinite(sphere.radius) || sphere.radius <= 0) {
                alert("Invalid or corrupted STL file. Geometry bounds are infinite or zero.");
                return;
            }

            geometry.center();

            addModel({
                id: uuidv4(),
                name: file.name,
                bufferGeometry: geometry,
                color: '#ff6b00',
                visible: true,
                position: [0, 0, 0],
                rotation: [0, 0, 0],
                scale: [1, 1, 1],
                meshVersion: 0,
            });
        };
        reader.readAsArrayBuffer(file);
        event.target.value = '';
    };

    const showStitchingPanel = transformMode === 'stitching' && selectedId;

    return (
        <div style={{
            width: 'var(--sidebar-width)',
            height: '100%',
            backgroundColor: 'var(--bg-panel)',
            borderLeft: '1px solid var(--border-color)',
            display: 'flex',
            flexDirection: 'column',
            overflowY: 'auto'
        }}>
            {/* Scene Hierarchy */}
            <div style={{ padding: 'var(--spacing-md)', borderBottom: '1px solid var(--border-color)' }}>
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '12px' }}>
                    <h3 style={{ fontSize: '12px', color: 'var(--text-secondary)', textTransform: 'uppercase' }}>Scene Hierarchy</h3>
                    <div style={{ display: 'flex', gap: '8px', alignItems: 'center' }}>
                        <button
                            className="mobile-only"
                            onClick={onClose}
                            style={{ color: 'var(--text-muted)' }}
                        >
                            <X size={16} />
                        </button>
                        <input
                            type="file"
                            ref={fileInputRef}
                            onChange={handleFileUpload}
                            accept=".stl"
                            style={{ display: 'none' }}
                        />
                        <Tooltip content="Import STL File" icon={<Plus size={12} />} align="left">
                            <button
                                onClick={() => fileInputRef.current?.click()}
                                style={{
                                    color: 'var(--accent-primary)',
                                    fontSize: '20px',
                                    fontWeight: 'bold',
                                    cursor: 'pointer',
                                    padding: '4px',
                                    display: 'flex',
                                    alignItems: 'center',
                                    justifyContent: 'center'
                                }}
                            >
                                <Plus size={18} />
                            </button>
                        </Tooltip>
                    </div>
                </div>

                <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
                    {models.map((model) => (
                        <div
                            key={model.id}
                            onClick={() => {
                                selectModel(model.id);
                                selectAllFaces(model.id);
                            }}
                            style={{
                                display: 'flex',
                                alignItems: 'center',
                                padding: '8px',
                                borderRadius: 'var(--radius-sm)',
                                backgroundColor: selectedId === model.id ? 'var(--accent-primary)' : 'transparent',
                                cursor: 'pointer',
                                color: selectedId === model.id ? '#fff' : 'var(--text-primary)'
                            }}
                        >
                            <Box size={14} style={{ marginRight: '8px', opacity: 0.7 }} />
                            <span style={{ flex: 1, fontSize: '13px', whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis' }}>
                                {model.name}
                            </span>
                            <Tooltip content={model.visible ? "Hide Model" : "Show Model"} align="left">
                                <button
                                    onClick={(e) => {
                                        e.stopPropagation();
                                        updateModel(model.id, { visible: !model.visible });
                                    }}
                                    style={{ padding: '2px', marginRight: '4px' }}
                                >
                                    {model.visible ? <Eye size={14} /> : <EyeOff size={14} />}
                                </button>
                            </Tooltip>
                            <Tooltip content="Delete Model" align="left">
                                <button
                                    onClick={(e) => {
                                        e.stopPropagation();
                                        removeModel(model.id);
                                    }}
                                    style={{ padding: '2px', color: 'var(--text-muted)' }}
                                    onMouseEnter={(e) => e.currentTarget.style.color = '#ff4444'}
                                    onMouseLeave={(e) => e.currentTarget.style.color = 'var(--text-muted)'}
                                >
                                    <Trash2 size={14} />
                                </button>
                            </Tooltip>
                        </div>
                    ))}
                    {models.length === 0 && (
                        <div style={{ padding: '20px', textAlign: 'center', color: 'var(--text-muted)', fontSize: '12px' }}>
                            No models imported
                        </div>
                    )}
                </div>
            </div>

            {/* Properties Panel (SHARED) */}
            {selectedModel && (
                <div style={{ padding: 'var(--spacing-md)' }}>
                    {/* Stitching Panel Header & Content */}
                    {showStitchingPanel ? (
                        <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
                            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '8px' }}>
                                <h3 style={{ fontSize: '12px', color: 'var(--text-secondary)', textTransform: 'uppercase' }}>Stitching Repair</h3>
                            </div>

                            <div style={{ padding: '12px', backgroundColor: 'var(--bg-panel-light)', borderRadius: '8px', border: '1px solid var(--border-color)' }}>
                                <h3 style={{ marginTop: 0, fontSize: '14px', color: 'var(--text-secondary)', marginBottom: '12px' }}>Mesh Analysis</h3>
                                {stitchStats ? (
                                    <div style={{ display: 'flex', flexDirection: 'column', gap: '8px' }}>
                                        <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                                            <span style={{ color: 'var(--text-secondary)' }}>Open Edges (Holes):</span>
                                            <span style={{ fontWeight: 'bold', color: stitchStats.openEdges > 0 ? 'var(--accent-error)' : 'var(--accent-success)' }}>
                                                {stitchStats.openEdges}
                                            </span>
                                        </div>
                                        <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                                            <span style={{ color: 'var(--text-secondary)' }}>Non-Manifold Edges:</span>
                                            <span style={{ fontWeight: 'bold', color: stitchStats.nonManifoldEdges > 0 ? 'var(--accent-warning)' : 'var(--accent-success)' }}>
                                                {stitchStats.nonManifoldEdges}
                                            </span>
                                        </div>
                                    </div>
                                ) : (
                                    <span style={{ color: 'var(--text-secondary)', fontStyle: 'italic' }}>Analyzing...</span>
                                )}
                            </div>

                            <Tooltip content="Launch Advanced Mesh Repair (Watertight)" icon={<Zap size={12} />} align="left">
                                <button
                                    onClick={handleRepair}
                                    disabled={!stitchStats || isRepairing || isProcessing}
                                    style={{
                                        padding: '12px',
                                        backgroundColor: 'var(--accent-primary)',
                                        color: 'white',
                                        border: 'none',
                                        borderRadius: '8px',
                                        fontWeight: 600,
                                        cursor: 'pointer',
                                        opacity: (!stitchStats || isRepairing) ? 0.5 : 1,
                                        display: 'flex',
                                        alignItems: 'center',
                                        justifyContent: 'center',
                                        gap: '8px',
                                        width: '100%'
                                    }}
                                >
                                    {isRepairing ? (
                                        <>
                                            <RotateCcw className="spin" size={16} /> Repairing...
                                        </>
                                    ) : (
                                        <>
                                            <Zap size={16} /> Repair & Stitch (Watertight)
                                        </>
                                    )}
                                </button>
                            </Tooltip>

                            {stitchStats && (
                                <div style={{
                                    padding: '12px',
                                    backgroundColor: 'rgba(0, 210, 255, 0.05)',
                                    borderRadius: '8px',
                                    border: '1px solid rgba(0, 210, 255, 0.2)',
                                    marginBottom: '16px'
                                }}>
                                    <div style={{ fontSize: '11px', color: 'var(--text-muted)', marginBottom: '8px', display: 'flex', justifyContent: 'space-between' }}>
                                        <span>REPAIR STATUS</span>
                                        <span style={{ color: '#00d2ff', fontWeight: 'bold' }}>
                                            {analyzedHoles ? `${analyzedHoles.length} Holes detected` : 'Analyzing...'}
                                        </span>
                                    </div>
                                    <p style={{ fontSize: '10px', color: 'var(--text-muted)', fontStyle: 'italic', margin: 0 }}>
                                        The automatic engine will analyze and close all geometric gaps to ensure the model is watertight for 3D printing.
                                    </p>
                                </div>
                            )}

                            <div style={{ display: 'flex', flexDirection: 'column', gap: '8px' }}>
                                <p style={{ fontSize: '11px', color: 'var(--text-secondary)', lineHeight: '1.4', margin: 0 }}>
                                    • <strong>Stitch & Repair:</strong> Legacy JS-based method for small cracks.
                                </p>
                                <p style={{ fontSize: '11px', color: 'var(--text-secondary)', lineHeight: '1.4', margin: 0 }}>
                                    • <strong>Fill Gaps:</strong> Conservative repair — fills small gaps with tiny triangles. Doesn't remove or modify existing geometry.
                                </p>
                                <p style={{ fontSize: '10px', color: 'var(--text-muted)', lineHeight: '1.4', marginTop: '4px', fontStyle: 'italic' }}>
                                    Increase the sliders if holes aren't being filled. The repair will never modify existing faces.
                                </p>
                            </div>
                        </div>
                    ) : (
                        /* Standard Properties Panel Content */
                        <>
                            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '16px' }}>
                                <h3 style={{ fontSize: '12px', color: 'var(--text-secondary)', textTransform: 'uppercase' }}>
                                    {isHistoryReEdit ? 'Operation Re-edit' : 'Properties'}
                                </h3>
                                {isHistoryReEdit && (
                                    <span style={{ fontSize: '10px', color: '#ff6b00', fontWeight: 'bold', textTransform: 'uppercase' }}>Past Action</span>
                                )}
                            </div>

                            {!isHistoryReEdit && (
                                <div style={{ display: 'flex', gap: '8px', alignItems: 'flex-end', marginBottom: '16px' }}>
                                    <div style={{ flex: 1 }}>
                                        <label style={{ display: 'block', fontSize: '11px', color: 'var(--text-muted)', marginBottom: '4px' }}>Name</label>
                                        <input
                                            type="text"
                                            value={selectedModel.name}
                                            onChange={(e) => updateModel(selectedModel.id, { name: e.target.value })}
                                            style={{
                                                width: '100%',
                                                backgroundColor: 'var(--bg-input)',
                                                border: '1px solid var(--border-color)',
                                                color: 'var(--text-primary)',
                                                padding: '6px 8px',
                                                borderRadius: 'var(--radius-sm)',
                                                fontSize: '13px',
                                                height: '32px'
                                            }}
                                        />
                                    </div>
                                    <div style={{ width: '40px' }}>
                                        <label style={{ display: 'block', fontSize: '11px', color: 'var(--text-muted)', marginBottom: '4px' }}>Color</label>
                                        <div style={{ position: 'relative', width: '32px', height: '32px' }}>
                                            <input
                                                type="color"
                                                value={selectedModel.color}
                                                onChange={(e) => updateModel(selectedModel.id, { color: e.target.value })}
                                                style={{
                                                    position: 'absolute', top: 0, left: 0,
                                                    width: '100%', height: '100%', padding: '0',
                                                    border: '1px solid var(--border-color)',
                                                    borderRadius: 'var(--radius-sm)',
                                                    cursor: 'pointer', opacity: 0
                                                }}
                                            />
                                            <div style={{
                                                width: '100%', height: '100%',
                                                backgroundColor: selectedModel.color,
                                                border: '1px solid var(--border-color)',
                                                borderRadius: 'var(--radius-sm)',
                                                pointerEvents: 'none'
                                            }} />
                                        </div>
                                    </div>
                                </div>
                            )}

                            {/* Subdivide Tool Section */}
                            {showSubdividePanel && (
                                <div style={{
                                    marginTop: '24px', padding: '16px',
                                    backgroundColor: isHistoryReEdit ? 'rgba(255, 107, 0, 0.1)' : 'rgba(255, 107, 0, 0.05)',
                                    borderRadius: 'var(--radius-md)',
                                    border: isHistoryReEdit ? '1px solid #ff6b00' : '1px solid rgba(255, 107, 0, 0.2)'
                                }}>
                                    <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '16px' }}>
                                        <Scissors size={18} style={{ color: isHistoryReEdit ? '#ff6b00' : 'var(--accent-primary)' }} />
                                        <h3 style={{ fontSize: '14px', fontWeight: '600', color: 'var(--text-primary)' }}>
                                            {isHistoryReEdit ? 'Re-calculate Subdivision' : 'Subdivide Tool'}
                                        </h3>
                                    </div>
                                    <div style={{ marginBottom: '16px' }}>
                                        <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                            <label style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>Steps</label>
                                            <span style={{ fontSize: '12px', color: isHistoryReEdit ? '#ff6b00' : 'var(--accent-primary)', fontWeight: '500' }}>{subdivideSteps}</span>
                                        </div>
                                        <input type="range" min="1" max="3" step="1" value={subdivideSteps}
                                            onChange={(e) => setSubdivideSteps(parseInt(e.target.value))}
                                            style={{ width: '100%', accentColor: isHistoryReEdit ? '#ff6b00' : 'var(--accent-primary)' }}
                                        />
                                    </div>
                                    <Tooltip content="Apply Loop Subdivision to selected faces" icon={<Scissors size={12} />} align="left">
                                        <button onClick={handleApply} disabled={isProcessing} style={{
                                            width: '100%', padding: '10px',
                                            backgroundColor: isHistoryReEdit ? '#ff6b00' : 'var(--accent-primary)',
                                            color: 'white', borderRadius: 'var(--radius-sm)',
                                            display: 'flex', alignItems: 'center', justifyContent: 'center',
                                            gap: '8px', fontSize: '13px', border: 'none',
                                            cursor: isProcessing ? 'not-allowed' : 'pointer',
                                            fontWeight: 'bold',
                                            opacity: isProcessing ? 0.6 : 1
                                        }}>
                                            {isHistoryReEdit ? <RotateCcw size={14} /> : <Scissors size={14} />}
                                            {isHistoryReEdit ? 'Recalculate & Reprocess' : 'Apply Subdivision'}
                                        </button>
                                    </Tooltip>
                                </div>
                            )}

                            {/* Texturize Tool Section */}
                            {showTexturizePanel && (
                                <div style={{
                                    marginTop: '24px', padding: '16px',
                                    backgroundColor: isHistoryReEdit ? 'rgba(0, 210, 255, 0.1)' : 'rgba(0, 210, 255, 0.05)',
                                    borderRadius: 'var(--radius-md)',
                                    border: isHistoryReEdit ? '1px solid #00d2ff' : '1px solid rgba(0, 210, 255, 0.2)'
                                }}>
                                    <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '16px' }}>
                                        <Grid3X3 size={18} style={{ color: isHistoryReEdit ? '#00d2ff' : 'var(--accent-primary)' }} />
                                        <h3 style={{ fontSize: '14px', fontWeight: '600', color: 'var(--text-primary)' }}>
                                            {isHistoryReEdit ? 'Update Texture' : 'Texturize Tool'}
                                        </h3>
                                    </div>

                                    <div style={{ marginBottom: '16px' }}>
                                        <div style={{
                                            fontSize: '11px', color: 'white', textTransform: 'uppercase',
                                            letterSpacing: '0.05em', backgroundColor: 'var(--accent-primary)',
                                            padding: '4px 8px', borderRadius: '4px', display: 'inline-block', fontWeight: 'bold'
                                        }}>
                                            Active: {textureType === 'decimate' ? 'Simplify' : textureType}
                                        </div>
                                    </div>

                                    <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
                                        {textureType === 'knurling' && (
                                            <>
                                                <div>
                                                    <label style={{ display: 'block', fontSize: '11px', color: 'var(--text-muted)', marginBottom: '8px', textTransform: 'uppercase', letterSpacing: '0.05em' }}>1. Knurl Pattern</label>
                                                    <div style={{ display: 'flex', gap: '4px' }}>
                                                        {(['diamond', 'straight', 'diagonal', 'square'] as const).map((p) => (
                                                            <button key={p} onClick={() => {
                                                                setKnurlPattern(p);
                                                                if (p === 'straight' && angle !== 0 && angle !== 90) setAngle(0);
                                                                if (p !== 'straight' && angle === 90) setAngle(45);
                                                            }} style={{
                                                                flex: 1, padding: '8px 4px', fontSize: '10px', borderRadius: '4px',
                                                                backgroundColor: knurlPattern === p ? 'var(--accent-primary)' : 'var(--bg-input)',
                                                                color: 'white', border: '1px solid var(--border-color)',
                                                                cursor: 'pointer', textTransform: 'capitalize', transition: 'all 0.2s'
                                                            }}>
                                                                {p}
                                                            </button>
                                                        ))}
                                                    </div>
                                                </div>

                                                <div>
                                                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                        <label style={{ fontSize: '11px', color: 'var(--text-muted)', textTransform: 'uppercase', letterSpacing: '0.05em' }}>2. Knurling Angle</label>
                                                        <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '600' }}>{angle}°</span>
                                                    </div>
                                                    <div style={{ display: 'flex', gap: '4px' }}>
                                                        {(knurlPattern === 'straight' ? [0, 90] : [0, 30, 45]).map((a) => (
                                                            <button key={a} onClick={() => setAngle(a)} style={{
                                                                flex: 1, padding: '8px', fontSize: '11px', borderRadius: '4px',
                                                                backgroundColor: angle === a ? (isHistoryReEdit ? '#00d2ff' : 'var(--accent-primary)') : 'var(--bg-input)',
                                                                color: 'white', border: '1px solid var(--border-color)',
                                                                cursor: 'pointer', transition: 'all 0.2s'
                                                            }}>
                                                                {a === 0 && knurlPattern === 'straight' ? 'Vertical' : (a === 90 ? 'Horizontal' : `${a}°`)}
                                                            </button>
                                                        ))}
                                                    </div>
                                                </div>

                                                <div>
                                                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                        <label style={{ fontSize: '11px', color: 'var(--text-muted)', textTransform: 'uppercase', letterSpacing: '0.05em' }}>3. Pitch (mm)</label>
                                                        <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '600' }}>{pitch.toFixed(2)}</span>
                                                    </div>
                                                    <input type="range" min="0.1" max="10" step="0.0001" value={pitch}
                                                        onChange={(e) => { const val = parseFloat(e.target.value); setPitch(val); setReferencePitch(val); }}
                                                        style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                                    />
                                                    {cylinderFit && (
                                                        <div style={{ marginTop: '8px', padding: '10px', backgroundColor: 'rgba(0,0,0,0.2)', borderRadius: '6px', border: '1px solid rgba(0,210,255,0.3)' }}>
                                                            <div style={{ fontSize: '11px', color: '#00d2ff', fontWeight: 'bold', marginBottom: '6px', letterSpacing: '0.05em' }}>PERFECT PITCH</div>
                                                            <div style={{ fontSize: '10px', color: 'var(--text-muted)', marginBottom: '8px', opacity: 0.8 }}>CYLINDER DETECTED (C: {circumference.toFixed(1)}mm)</div>
                                                            <div style={{ display: 'flex', gap: '6px', flexWrap: 'wrap' }}>
                                                                {(() => {
                                                                    const isEvenPattern = knurlPattern === 'diamond' || knurlPattern === 'square';
                                                                    const n_nom = Math.round(circumference / referencePitch);
                                                                    let n_base = n_nom;
                                                                    if (isEvenPattern && n_base % 2 !== 0) n_base += 1;
                                                                    const candidates = new Set<number>();
                                                                    candidates.add(n_base);
                                                                    candidates.add(n_base + (isEvenPattern ? 2 : 1));
                                                                    candidates.add(n_base - (isEvenPattern ? 2 : 1));
                                                                    const sorted = Array.from(candidates).filter(n => n >= 2).sort((a, b) => a - b);
                                                                    return sorted.map(n => {
                                                                        const sugP = circumference / n;
                                                                        const isSelected = Math.abs(pitch - sugP) < 0.001;
                                                                        return (
                                                                            <button key={n} onClick={() => setPitch(sugP)} style={{
                                                                                padding: '4px 8px', fontSize: '10px',
                                                                                backgroundColor: isSelected ? 'var(--accent-primary)' : 'rgba(0, 210, 255, 0.2)',
                                                                                color: isSelected ? 'white' : '#00d2ff',
                                                                                border: isSelected ? '1px solid var(--accent-primary)' : '1px solid rgba(0, 210, 255, 0.4)',
                                                                                borderRadius: '4px', cursor: 'pointer', fontWeight: '500', transition: 'all 0.2s'
                                                                            }}>
                                                                                N={n} → {sugP.toFixed(2)}mm
                                                                            </button>
                                                                        );
                                                                    });
                                                                })()}
                                                            </div>
                                                        </div>
                                                    )}
                                                </div>

                                                <div>
                                                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                        <label style={{ fontSize: '11px', color: 'var(--text-muted)', textTransform: 'uppercase', letterSpacing: '0.05em' }}>4. Depth (mm)</label>
                                                        <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '600' }}>{depth.toFixed(2)}</span>
                                                    </div>
                                                    <input type="range" min="0.1" max="5" step="0.1" value={depth}
                                                        onChange={(e) => setDepth(parseFloat(e.target.value))}
                                                        style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                                    />
                                                </div>
                                            </>
                                        )}

                                        {textureType === 'honeycomb' && (
                                            <>
                                                <div>
                                                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                        <label style={{ fontSize: '11px', color: 'var(--text-muted)', textTransform: 'uppercase', letterSpacing: '0.05em' }}>1. Pattern Angle</label>
                                                        <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '600' }}>{angle}°</span>
                                                    </div>
                                                    <div style={{ display: 'flex', gap: '4px' }}>
                                                        {[0, 30, 45, 90].map((a) => (
                                                            <button key={a} onClick={() => setAngle(a)} style={{
                                                                flex: 1, padding: '8px', fontSize: '11px', borderRadius: '4px',
                                                                backgroundColor: angle === a ? 'var(--accent-primary)' : 'var(--bg-input)',
                                                                color: 'white', border: '1px solid var(--border-color)', cursor: 'pointer'
                                                            }}>
                                                                {a}°
                                                            </button>
                                                        ))}
                                                    </div>
                                                </div>

                                                <div>
                                                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                        <label style={{ fontSize: '11px', color: 'var(--text-muted)', textTransform: 'uppercase', letterSpacing: '0.05em' }}>2. Cell Size (mm)</label>
                                                        <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '600' }}>{cellSize.toFixed(2)}</span>
                                                    </div>
                                                    <input type="range" min="1" max="20" step="0.1" value={cellSize}
                                                        onChange={(e) => { const val = parseFloat(e.target.value); setCellSize(val); setReferenceCellSize(val); }}
                                                        style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                                    />
                                                    {cylinderFit && (
                                                        <div style={{ marginTop: '8px', padding: '10px', backgroundColor: 'rgba(0,0,0,0.2)', borderRadius: '6px', border: '1px solid rgba(0,210,255,0.3)' }}>
                                                            <div style={{ fontSize: '11px', color: '#00d2ff', fontWeight: 'bold', marginBottom: '6px', letterSpacing: '0.05em' }}>PERFECT CELL SIZE</div>
                                                            <div style={{ fontSize: '10px', color: 'var(--text-muted)', marginBottom: '8px', opacity: 0.8 }}>CYLINDER DETECTED (C: {circumference.toFixed(1)}mm)</div>
                                                            <div style={{ display: 'flex', gap: '6px', flexWrap: 'wrap' }}>
                                                                {[0, 1].map(isAlt => {
                                                                    const ang = ((angle % 180) + 180) % 180;
                                                                    const sym = [0, 30, 60, 90, 120, 150, 180].reduce((p, c) => Math.abs(c - ang) < Math.abs(p - ang) ? c : p);
                                                                    const factor = (sym % 60 === 0) ? 1.0 : Math.sqrt(3);
                                                                    const n = Math.round(circumference / (referenceCellSize * factor)) + (isAlt ? 1 : 0);
                                                                    const sugW = circumference / (Math.max(1, n) * factor);
                                                                    if (sugW < 1 || sugW > 25) return null;
                                                                    const isSelected = Math.abs(cellSize - sugW) < 0.001;
                                                                    return (
                                                                        <button key={isAlt} onClick={() => setCellSize(sugW)} style={{
                                                                            padding: '4px 8px', fontSize: '10px',
                                                                            backgroundColor: isSelected ? 'var(--accent-primary)' : 'rgba(0, 210, 255, 0.2)',
                                                                            color: isSelected ? 'white' : '#00d2ff',
                                                                            border: isSelected ? '1px solid var(--accent-primary)' : '1px solid rgba(0, 210, 255, 0.4)',
                                                                            borderRadius: '4px', cursor: 'pointer', fontWeight: '500', transition: 'all 0.2s'
                                                                        }}>
                                                                            N={n} → {sugW.toFixed(2)}mm
                                                                        </button>
                                                                    );
                                                                })}
                                                            </div>
                                                        </div>
                                                    )}
                                                </div>

                                                <div>
                                                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                        <label style={{ fontSize: '11px', color: 'var(--text-muted)', textTransform: 'uppercase', letterSpacing: '0.05em' }}>3. Wall Thickness (mm)</label>
                                                        <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '600' }}>{wallThickness.toFixed(2)}</span>
                                                    </div>
                                                    <input type="range" min="0.05" max="2.0" step="0.05" value={wallThickness}
                                                        onChange={(e) => setWallThickness(parseFloat(e.target.value))}
                                                        style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                                    />
                                                </div>

                                                <div>
                                                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                        <label style={{ fontSize: '11px', color: 'var(--text-muted)', textTransform: 'uppercase', letterSpacing: '0.05em' }}>4. Depth (mm)</label>
                                                        <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '600' }}>{depth.toFixed(2)}</span>
                                                    </div>
                                                    <input type="range" min="0.1" max="5" step="0.1" value={depth}
                                                        onChange={(e) => setDepth(parseFloat(e.target.value))}
                                                        style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                                    />
                                                </div>

                                                <div>
                                                    <label style={{ display: 'block', fontSize: '11px', color: 'var(--text-muted)', marginBottom: '8px', textTransform: 'uppercase', letterSpacing: '0.05em' }}>5. Direction</label>
                                                    <div style={{ display: 'flex', gap: '4px' }}>
                                                        {(['inward', 'outward'] as const).map((d) => (
                                                            <button key={d} onClick={() => setDirection(d)} style={{
                                                                flex: 1, padding: '8px', fontSize: '11px', borderRadius: '4px',
                                                                backgroundColor: direction === d ? 'var(--accent-primary)' : 'var(--bg-input)',
                                                                color: 'white', border: '1px solid var(--border-color)',
                                                                cursor: 'pointer', textTransform: 'capitalize'
                                                            }}>
                                                                {d}
                                                            </button>
                                                        ))}
                                                    </div>
                                                </div>
                                            </>
                                        )}

                                        {textureType === 'fuzzy' && (
                                            <>
                                                <div>
                                                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                        <label style={{ fontSize: '11px', color: 'var(--text-muted)', textTransform: 'uppercase', letterSpacing: '0.05em' }}>1. Noise Thickness (mm)</label>
                                                        <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '600' }}>{fuzzyThickness.toFixed(2)}</span>
                                                    </div>
                                                    <input type="range" min="0.2" max="5.0" step="0.05" value={fuzzyThickness}
                                                        onChange={(e) => setFuzzyThickness(parseFloat(e.target.value))}
                                                        style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                                    />
                                                </div>
                                                <div>
                                                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                        <label style={{ fontSize: '11px', color: 'var(--text-muted)', textTransform: 'uppercase', letterSpacing: '0.05em' }}>2. Noise Density (mm dist)</label>
                                                        <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '600' }}>{pointDistance.toFixed(2)}</span>
                                                    </div>
                                                    <input type="range" min="0.2" max="5.0" step="0.05" value={pointDistance}
                                                        onChange={(e) => setPointDistance(parseFloat(e.target.value))}
                                                        style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                                    />
                                                    <div style={{ fontSize: '10px', color: 'var(--text-muted)', marginTop: '4px', fontStyle: 'italic' }}>Lower values = Higher resolution noise</div>

                                                    <div style={{
                                                        marginTop: '12px',
                                                        padding: '8px',
                                                        backgroundColor: 'rgba(255, 165, 0, 0.1)',
                                                        borderRadius: '4px',
                                                        border: '1px solid rgba(255, 165, 0, 0.3)',
                                                        display: 'flex',
                                                        gap: '8px',
                                                        alignItems: 'flex-start'
                                                    }}>
                                                        <Zap size={12} style={{ color: '#ffa500', marginTop: '2px', flexShrink: 0 }} />
                                                        <p style={{ fontSize: '10px', color: '#ffa500', margin: 0, lineHeight: '1.3' }}>
                                                            <strong>Note:</strong> High resolution noise (low density) significantly increases triangle count, which will result in larger file sizes and may affect real-time performance.
                                                        </p>
                                                    </div>
                                                </div>
                                            </>
                                        )}

                                        {(textureType === 'knurling' || textureType === 'honeycomb' || textureType === 'fuzzy') && (
                                            <div style={{ padding: '12px', backgroundColor: 'rgba(255, 107, 0, 0.05)', borderRadius: '4px', border: '1px dashed rgba(255, 107, 0, 0.2)' }}>
                                                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                                                    <label style={{ fontSize: '11px', color: 'var(--text-muted)', textTransform: 'uppercase', letterSpacing: '0.05em' }}>Watertight Seal</label>
                                                    <button onClick={() => setHoleFillEnabled(!holeFillEnabled)} style={{
                                                        padding: '4px 8px', fontSize: '10px',
                                                        backgroundColor: holeFillEnabled ? 'var(--accent-primary)' : 'var(--bg-tertiary)',
                                                        color: 'white', border: 'none', borderRadius: '4px', cursor: 'pointer'
                                                    }}>
                                                        {holeFillEnabled ? 'ON' : 'OFF'}
                                                    </button>
                                                </div>
                                            </div>
                                        )}

                                        {textureType === 'decimate' && (
                                            <div>
                                                <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                    <label style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>Reduction Ratio</label>
                                                    <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '500' }}>{(reductionRatio * 100).toFixed(0)}%</span>
                                                </div>
                                                <input type="range" min="0.01" max="0.99" step="0.01" value={reductionRatio}
                                                    onChange={(e) => setReductionRatio(parseFloat(e.target.value))}
                                                    style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                                />
                                            </div>
                                        )}

                                        <Tooltip content="Bake texture into mesh geometry" icon={<Grid3X3 size={12} />} align="left">
                                            <button onClick={handleApply} disabled={isProcessing} style={{
                                                width: '100%', padding: '12px',
                                                backgroundColor: isHistoryReEdit ? '#00d2ff' : 'var(--accent-primary)',
                                                color: 'white', border: 'none', borderRadius: 'var(--radius-sm)',
                                                cursor: isProcessing ? 'not-allowed' : 'pointer',
                                                fontWeight: 'bold', fontSize: '13px', marginTop: '16px',
                                                opacity: isProcessing ? 0.6 : 1
                                            }}>
                                                {isHistoryReEdit ? <RotateCcw size={14} style={{ marginRight: '8px' }} /> : null}
                                                {isHistoryReEdit ? 'Recalculate & Redo' : 'Apply Texture'}
                                            </button>
                                        </Tooltip>
                                    </div>
                                </div>
                            )}
                        </>
                    )}
                </div>
            )}
        </div>
    );
};
