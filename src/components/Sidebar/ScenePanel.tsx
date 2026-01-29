import React, { useState, useRef } from 'react';
import { useStore } from '../../store/useStore';
import type { KnurlPattern } from '../../utils/texturizeUtils';
import { Eye, EyeOff, Box, Trash2, Scissors, Grid3X3, X, Info, RotateCcw } from 'lucide-react';
import { estimateSelectionCircumference } from '../../utils/meshUtils';
import { v4 as uuidv4 } from 'uuid';
import { STLLoader } from 'three/examples/jsm/loaders/STLLoader.js';

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

    const fileInputRef = useRef<HTMLInputElement>(null);
    const [subdivideSteps, setSubdivideSteps] = useState(1);

    // Texturize states
    const [textureType, setTextureType] = useState<'knurling' | 'honeycomb' | 'decimate'>('knurling');
    const [pitch, setPitch] = useState(2.0);
    const [depth, setDepth] = useState(0.5);
    const [angle, setAngle] = useState(45);
    const [cellSize, setCellSize] = useState(5.0);
    const [wallThickness, setWallThickness] = useState(0.5);
    const [knurlPattern, setKnurlPattern] = useState<KnurlPattern>('diamond');
    const [reductionRatio, setReductionRatio] = useState(0.5);

    // Sync with re-edit params
    React.useEffect(() => {
        if (reEditParams) {
            if (transformMode === 'subdivide') {
                if (reEditParams.steps !== undefined) setSubdivideSteps(reEditParams.steps);
            } else if (transformMode === 'texturize') {
                if (reEditParams.type) setTextureType(reEditParams.type);
                if (reEditParams.pitch !== undefined) setPitch(reEditParams.pitch);
                if (reEditParams.depth !== undefined) setDepth(reEditParams.depth);
                if (reEditParams.angle !== undefined) setAngle(reEditParams.angle);
                if (reEditParams.pattern) setKnurlPattern(reEditParams.pattern);
                if (reEditParams.cellSize !== undefined) setCellSize(reEditParams.cellSize);
                if (reEditParams.wallThickness !== undefined) setWallThickness(reEditParams.wallThickness);
                if (reEditParams.reduction !== undefined) setReductionRatio(reEditParams.reduction);
            }
        }
    }, [reEditParams, transformMode]);

    const isHistoryReEdit = selectedHistoryIds.length === 1 && selectedHistoryIds[0] !== 'initial' && reEditParams !== null;
    const selectedModel = models.find((m) => m.id === selectedId);

    // In history re-edit mode, we might not have a face selection in the current live mesh, 
    // but the button should still show if we are browsing history.
    const selection = selectedId ? smartSelection[selectedId] || [] : [];
    const selectionCount = selection.length;
    const hasSelection = selectionCount > 0;

    const showSubdividePanel = transformMode === 'subdivide' && (hasSelection || isHistoryReEdit);
    const showTexturizePanel = transformMode === 'texturize' && (hasSelection || isHistoryReEdit);

    const circumference = React.useMemo(() => {
        if (!selectedModel || selection.length === 0) return 0;
        return estimateSelectionCircumference(selectedModel.bufferGeometry, selection);
    }, [selectedModel, selection]);

    const suggestedPitches = React.useMemo(() => {
        if (circumference <= 0) return [];
        const suggestions = [];
        const angRad = (angle * Math.PI) / 180;
        const factor = Math.cos(angRad) + Math.sin(angRad);
        for (let n = 20; n <= 100; n += 10) {
            const hPitch = circumference / n;
            suggestions.push(hPitch / factor);
        }
        return suggestions;
    }, [circumference, angle]);

    const handleApply = () => {
        if (isHistoryReEdit) {
            // Recalculate existing operation
            const params = transformMode === 'subdivide'
                ? { steps: subdivideSteps }
                : { type: textureType, pitch, depth, angle, pattern: knurlPattern, cellSize, wallThickness, reduction: reductionRatio };
            recalculateHistoryItem(selectedHistoryIds[0], params);
        } else if (selectedModel) {
            // Apply new operation
            if (transformMode === 'subdivide') {
                subdivideSelection(selectedModel.id, subdivideSteps);
            } else if (transformMode === 'texturize') {
                if (textureType === 'knurling') {
                    applyTexturize(selectedModel.id, { type: 'knurling', pitch, depth, angle, pattern: knurlPattern });
                } else if (textureType === 'honeycomb') {
                    applyTexturize(selectedModel.id, { type: 'honeycomb', cellSize, wallThickness, depth });
                } else {
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

    return (
        <div style={{
            width: 'var(--sidebar-width)',
            backgroundColor: 'var(--bg-panel)',
            borderLeft: '1px solid var(--border-color)',
            display: 'flex',
            flexDirection: 'column',
            height: '100%',
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
                            +
                        </button>
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
                            <button
                                onClick={(e) => {
                                    e.stopPropagation();
                                    updateModel(model.id, { visible: !model.visible });
                                }}
                                style={{ padding: '2px', marginRight: '4px' }}
                            >
                                {model.visible ? <Eye size={14} /> : <EyeOff size={14} />}
                            </button>
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
                        </div>
                    ))}
                    {models.length === 0 && (
                        <div style={{ padding: '20px', textAlign: 'center', color: 'var(--text-muted)', fontSize: '12px' }}>
                            No models imported
                        </div>
                    )}
                </div>
            </div>

            {/* Properties Panel */}
            {selectedModel && (
                <div style={{ padding: 'var(--spacing-md)' }}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '16px' }}>
                        <h3 style={{ fontSize: '12px', color: 'var(--text-secondary)', textTransform: 'uppercase' }}>
                            {isHistoryReEdit ? 'Operation Re-edit' : 'Properties'}
                        </h3>
                        {isHistoryReEdit && (
                            <span style={{ fontSize: '10px', color: '#ff6b00', fontWeight: 'bold', textTransform: 'uppercase' }}>Past Action</span>
                        )}
                    </div>

                    {!isHistoryReEdit && (
                        <div style={{ display: 'flex', flexDirection: 'column', gap: '12px' }}>
                            <div>
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
                                        fontSize: '13px'
                                    }}
                                />
                            </div>
                            <div>
                                <label style={{ display: 'block', fontSize: '11px', color: 'var(--text-muted)', marginBottom: '4px' }}>Color</label>
                                <input
                                    type="color"
                                    value={selectedModel.color}
                                    onChange={(e) => updateModel(selectedModel.id, { color: e.target.value })}
                                    style={{
                                        width: '100%',
                                        height: '30px',
                                        padding: '0',
                                        border: '1px solid var(--border-color)',
                                        borderRadius: 'var(--radius-sm)',
                                        cursor: 'pointer'
                                    }}
                                />
                            </div>
                        </div>
                    )}

                    {/* Subdivide Tool Section */}
                    {showSubdividePanel && (
                        <div style={{
                            marginTop: '24px',
                            padding: '16px',
                            backgroundColor: isHistoryReEdit ? 'rgba(255, 107, 0, 0.1)' : 'rgba(255, 107, 0, 0.05)',
                            borderRadius: 'var(--radius-md)',
                            border: isHistoryReEdit ? '1px solid #ff6b00' : '1px solid rgba(255, 107, 0, 0.2)'
                        }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '16px' }}>
                                <div style={{ color: isHistoryReEdit ? '#ff6b00' : 'var(--accent-primary)' }}>
                                    <Scissors size={18} />
                                </div>
                                <h3 style={{ fontSize: '14px', fontWeight: '600', color: 'var(--text-primary)' }}>
                                    {isHistoryReEdit ? 'Re-calculate Subdivision' : 'Subdivide Tool'}
                                </h3>
                            </div>

                            <div style={{ marginBottom: '16px' }}>
                                <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                    <label style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>Steps</label>
                                    <span style={{ fontSize: '12px', color: isHistoryReEdit ? '#ff6b00' : 'var(--accent-primary)', fontWeight: '500' }}>{subdivideSteps}</span>
                                </div>
                                <input
                                    type="range"
                                    min="1"
                                    max="3"
                                    step="1"
                                    value={subdivideSteps}
                                    onChange={(e) => setSubdivideSteps(parseInt(e.target.value))}
                                    style={{ width: '100%', accentColor: isHistoryReEdit ? '#ff6b00' : 'var(--accent-primary)' }}
                                />
                            </div>

                            <button
                                onClick={handleApply}
                                style={{
                                    width: '100%',
                                    padding: '10px',
                                    backgroundColor: isHistoryReEdit ? '#ff6b00' : 'var(--accent-primary)',
                                    color: 'white',
                                    borderRadius: 'var(--radius-sm)',
                                    display: 'flex',
                                    alignItems: 'center',
                                    justifyContent: 'center',
                                    gap: '8px',
                                    fontSize: '13px',
                                    border: 'none',
                                    cursor: 'pointer',
                                    fontWeight: 'bold'
                                }}
                            >
                                {isHistoryReEdit ? <RotateCcw size={14} /> : <Scissors size={14} />}
                                {isHistoryReEdit ? 'Recalculate & Reprocess' : 'Apply Subdivision'}
                            </button>
                        </div>
                    )}

                    {/* Texturize Tool Section */}
                    {showTexturizePanel && (
                        <div style={{
                            marginTop: '24px',
                            padding: '16px',
                            backgroundColor: isHistoryReEdit ? 'rgba(0, 210, 255, 0.1)' : 'rgba(0, 210, 255, 0.05)',
                            borderRadius: 'var(--radius-md)',
                            border: isHistoryReEdit ? '1px solid #00d2ff' : '1px solid rgba(0, 210, 255, 0.2)'
                        }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '16px' }}>
                                <div style={{ color: isHistoryReEdit ? '#00d2ff' : 'var(--accent-primary)' }}>
                                    <Grid3X3 size={18} />
                                </div>
                                <h3 style={{ fontSize: '14px', fontWeight: '600', color: 'var(--text-primary)' }}>
                                    {isHistoryReEdit ? 'Update Texture' : 'Texturize Tool'}
                                </h3>
                            </div>

                            <div style={{ display: 'flex', gap: '4px', marginBottom: '16px', backgroundColor: 'var(--bg-input)', padding: '2px', borderRadius: '4px' }}>
                                {(['knurling', 'honeycomb', 'decimate'] as const).map((t) => (
                                    <button
                                        key={t}
                                        onClick={() => setTextureType(t)}
                                        style={{
                                            flex: 1,
                                            padding: '6px',
                                            fontSize: '11px',
                                            borderRadius: '3px',
                                            backgroundColor: textureType === t ? (isHistoryReEdit ? '#00d2ff' : 'var(--accent-primary)') : 'transparent',
                                            color: 'white',
                                            border: 'none',
                                            cursor: 'pointer',
                                            transition: 'all 0.2s',
                                            textTransform: 'capitalize'
                                        }}
                                    >
                                        {t === 'decimate' ? 'Simplify' : t}
                                    </button>
                                ))}
                            </div>

                            <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
                                {textureType === 'knurling' && (
                                    <>
                                        <div>
                                            <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                <label style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>Pitch (mm)</label>
                                                <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '500' }}>{pitch.toFixed(2)}</span>
                                            </div>
                                            <input
                                                type="range"
                                                min="0.1"
                                                max="10"
                                                step="0.01"
                                                value={pitch}
                                                onChange={(e) => setPitch(parseFloat(e.target.value))}
                                                style={{ width: '100%', accentColor: 'var(--accent-primary)', marginBottom: '12px' }}
                                            />
                                        </div>
                                        <div>
                                            <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                <label style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>Knurling Angle</label>
                                                <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '500' }}>{angle}°</span>
                                            </div>
                                            <div style={{ display: 'flex', gap: '4px', marginBottom: '16px' }}>
                                                {[0, 30, 45].map((a) => (
                                                    <button
                                                        key={a}
                                                        onClick={() => setAngle(a)}
                                                        style={{
                                                            flex: 1,
                                                            padding: '8px',
                                                            fontSize: '11px',
                                                            borderRadius: '4px',
                                                            backgroundColor: angle === a ? (isHistoryReEdit ? '#00d2ff' : 'var(--accent-primary)') : 'var(--bg-input)',
                                                            color: 'white',
                                                            border: '1px solid var(--border-color)',
                                                            cursor: 'pointer'
                                                        }}
                                                    >
                                                        {a === 0 ? 'Straight' : `${a}°`}
                                                    </button>
                                                ))}
                                            </div>
                                        </div>
                                        <div>
                                            <label style={{ display: 'block', fontSize: '11px', color: 'var(--text-muted)', marginBottom: '8px' }}>Knurl Pattern</label>
                                            <select
                                                value={knurlPattern}
                                                onChange={(e) => setKnurlPattern(e.target.value as KnurlPattern)}
                                                style={{
                                                    width: '100%',
                                                    backgroundColor: 'var(--bg-input)',
                                                    border: '1px solid var(--border-color)',
                                                    color: 'white',
                                                    padding: '8px',
                                                    borderRadius: '4px'
                                                }}
                                            >
                                                <option value="diamond">Diamond</option>
                                                <option value="straight">Straight</option>
                                                <option value="diagonal">Diagonal</option>
                                                <option value="square">Square</option>
                                            </select>
                                        </div>
                                    </>
                                )}

                                {textureType === 'honeycomb' && (
                                    <>
                                        <div>
                                            <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                <label style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>Cell Size</label>
                                                <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '500' }}>{cellSize}</span>
                                            </div>
                                            <input
                                                type="range"
                                                min="1"
                                                max="20"
                                                step="0.5"
                                                value={cellSize}
                                                onChange={(e) => setCellSize(parseFloat(e.target.value))}
                                                style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                            />
                                        </div>
                                    </>
                                )}

                                <div>
                                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                        <label style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>Depth (mm)</label>
                                        <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '500' }}>{depth}</span>
                                    </div>
                                    <input
                                        type="range"
                                        min="0.1"
                                        max="5"
                                        step="0.1"
                                        value={depth}
                                        onChange={(e) => setDepth(parseFloat(e.target.value))}
                                        style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                    />
                                </div>

                                <button
                                    onClick={handleApply}
                                    style={{
                                        width: '100%',
                                        padding: '12px',
                                        backgroundColor: isHistoryReEdit ? '#00d2ff' : 'var(--accent-primary)',
                                        color: 'white',
                                        border: 'none',
                                        borderRadius: 'var(--radius-sm)',
                                        cursor: 'pointer',
                                        fontWeight: 'bold',
                                        fontSize: '13px',
                                        marginTop: '16px'
                                    }}
                                >
                                    {isHistoryReEdit ? <RotateCcw size={14} style={{ marginRight: '8px' }} /> : null}
                                    {isHistoryReEdit ? 'Recalculate & Redo' : 'Apply Texture'}
                                </button>
                            </div>
                        </div>
                    )}
                </div>
            )}
        </div>
    );
};
