
import React, { useState } from 'react';
import { useStore } from '../../store/useStore';
import type { KnurlPattern } from '../../utils/texturizeUtils';
import { Eye, EyeOff, Box, Trash2, Scissors, Grid3X3, X, Info } from 'lucide-react';
import { estimateSelectionCircumference } from '../../utils/meshUtils';

interface ScenePanelProps {
    onClose?: () => void;
}

export const ScenePanel: React.FC<ScenePanelProps> = ({ onClose }) => {
    const models = useStore((state) => state.models);
    const selectedId = useStore((state) => state.selectedId);
    const selectModel = useStore((state) => state.selectModel);
    const updateModel = useStore((state) => state.updateModel);
    const removeModel = useStore((state) => state.removeModel);
    const subdivideSelection = useStore((state) => state.subdivideSelection);
    const applyTexturize = useStore((state) => state.applyTexturize);
    const selectAllFaces = useStore((state) => state.selectAllFaces);
    const smartSelection = useStore((state) => state.smartSelection);
    const transformMode = useStore((state) => state.transformMode);

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



    const selectedModel = models.find((m) => m.id === selectedId);
    const selection = selectedId ? smartSelection[selectedId] || [] : [];
    const hasSelection = selection.length > 0;
    const selectionCount = selection.length;
    const showSubdividePanel = transformMode === 'subdivide' && hasSelection;
    const showTexturizePanel = transformMode === 'texturize' && hasSelection;

    // Circumference estimation for seamless knurling
    const circumference = React.useMemo(() => {
        if (!selectedModel || selection.length === 0) return 0;
        return estimateSelectionCircumference(selectedModel.bufferGeometry, selection);
    }, [selectedModel, selection]);

    const suggestedPitches = React.useMemo(() => {
        if (circumference <= 0) return [];
        const suggestions = [];
        // Factor for diagonal pitch alignment
        // For a square grid rotated by A, the horizontal projection of one repeat unit
        // is P * (cos A + sin A)
        const angRad = (angle * Math.PI) / 180;
        const factor = Math.cos(angRad) + Math.sin(angRad);

        // Show pitches that result in an integer number of divisions around the circumference
        for (let n = 20; n <= 100; n += 10) {
            const hPitch = circumference / n;
            suggestions.push(hPitch / factor);
        }
        return suggestions;
    }, [circumference, angle]);

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
                        <button className="mobile-only" onClick={onClose} style={{ color: 'var(--text-muted)' }}>
                            <X size={16} />
                        </button>
                        <button style={{ color: 'var(--accent-primary)', fontSize: '18px' }}>+</button>
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
                    <h3 style={{ fontSize: '12px', color: 'var(--text-secondary)', textTransform: 'uppercase', marginBottom: '16px' }}>Properties</h3>

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

                    {/* Subdivide Tool Section */}
                    {showSubdividePanel && (
                        <div style={{
                            marginTop: '24px',
                            padding: '16px',
                            backgroundColor: 'rgba(255, 107, 0, 0.05)',
                            borderRadius: 'var(--radius-md)',
                            border: '1px solid rgba(255, 107, 0, 0.2)'
                        }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '16px' }}>
                                <div style={{ color: 'var(--accent-primary)' }}>
                                    <Scissors size={18} />
                                </div>
                                <h3 style={{ fontSize: '14px', fontWeight: '600', color: 'var(--text-primary)' }}>Subdivide Tool</h3>
                            </div>

                            <p style={{ fontSize: '12px', color: 'var(--text-secondary)', marginBottom: '16px' }}>
                                Refining the selection increases geometry density for better texture application.
                            </p>

                            <div style={{ marginBottom: '16px' }}>
                                <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                    <label style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>Steps</label>
                                    <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '500' }}>{subdivideSteps}</span>
                                </div>
                                <input
                                    type="range"
                                    min="1"
                                    max="3"
                                    step="1"
                                    value={subdivideSteps}
                                    onChange={(e) => setSubdivideSteps(parseInt(e.target.value))}
                                    style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                />
                            </div>

                            <button
                                onClick={() => subdivideSelection(selectedModel.id, subdivideSteps)}
                                style={{
                                    width: '100%',
                                    padding: '10px',
                                    backgroundColor: 'var(--accent-primary)',
                                    color: 'white',
                                    borderRadius: 'var(--radius-sm)',
                                    display: 'flex',
                                    alignItems: 'center',
                                    justifyContent: 'center',
                                    gap: '8px',
                                    fontSize: '13px',
                                    border: 'none',
                                    cursor: 'pointer'
                                }}
                            >
                                <Scissors size={14} />
                                Apply Subdivision
                            </button>
                        </div>
                    )}

                    {/* Texturize Tool Section */}
                    {showTexturizePanel && (
                        <div style={{
                            marginTop: '24px',
                            padding: '16px',
                            backgroundColor: 'rgba(0, 210, 255, 0.05)',
                            borderRadius: 'var(--radius-md)',
                            border: '1px solid rgba(0, 210, 255, 0.2)'
                        }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '16px' }}>
                                <div style={{ color: 'var(--accent-primary)' }}>
                                    <Grid3X3 size={18} />
                                </div>
                                <h3 style={{ fontSize: '14px', fontWeight: '600', color: 'var(--text-primary)' }}>Texturize Tool</h3>
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
                                            backgroundColor: textureType === t ? 'var(--accent-primary)' : 'transparent',
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

                            <p style={{ fontSize: '12px', color: 'var(--text-secondary)', marginBottom: '16px' }}>
                                Applying <strong>{textureType}</strong> to {selectionCount} faces.
                            </p>

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

                                            {circumference > 0 && (
                                                <div style={{
                                                    padding: '10px',
                                                    backgroundColor: 'rgba(255,107,0,0.1)',
                                                    borderRadius: '4px',
                                                    border: '1px solid rgba(255,107,0,0.2)',
                                                    marginBottom: '16px'
                                                }}>
                                                    <div style={{ display: 'flex', alignItems: 'center', gap: '6px', marginBottom: '8px' }}>
                                                        <Info size={12} color="var(--accent-primary)" />
                                                        <span style={{ fontSize: '11px', fontWeight: 'bold' }}>Seamless Suggestions</span>
                                                    </div>
                                                    <div style={{ display: 'flex', flexWrap: 'wrap', gap: '4px' }}>
                                                        {suggestedPitches.map((p, idx) => (
                                                            <button
                                                                key={idx}
                                                                onClick={() => setPitch(p)}
                                                                style={{
                                                                    fontSize: '10px',
                                                                    padding: '2px 6px',
                                                                    backgroundColor: 'var(--bg-input)',
                                                                    border: Math.abs(pitch - p) < 0.05 ? '1px solid var(--accent-primary)' : '1px solid transparent',
                                                                    borderRadius: '2px',
                                                                    color: 'var(--text-secondary)'
                                                                }}
                                                            >
                                                                {p.toFixed(2)}
                                                            </button>
                                                        ))}
                                                    </div>
                                                    <p style={{ fontSize: '9px', color: 'var(--text-muted)', marginTop: '8px' }}>
                                                        Est. Perimeter: {circumference.toFixed(1)}mm.
                                                    </p>
                                                </div>
                                            )}
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
                                                            backgroundColor: angle === a ? 'var(--accent-primary)' : 'var(--bg-input)',
                                                            color: 'white',
                                                            border: '1px solid var(--border-color)',
                                                            cursor: 'pointer',
                                                            transition: 'all 0.2s'
                                                        }}
                                                    >
                                                        {a === 0 ? 'Straight (0°)' : `${a}°`}
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
                                                    color: 'var(--text-primary)',
                                                    padding: '8px',
                                                    borderRadius: 'var(--radius-sm)',
                                                    fontSize: '13px',
                                                    outline: 'none',
                                                    cursor: 'pointer'
                                                }}
                                            >
                                                <option value="diamond">Diamond (Cross-Cut)</option>
                                                <option value="straight">Straight (Parallel)</option>
                                                <option value="diagonal">Diagonal (Single)</option>
                                                <option value="square">Square (Grid)</option>
                                            </select>
                                        </div>
                                    </>
                                )}

                                {textureType === 'honeycomb' && (
                                    <>
                                        <div>
                                            <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                <label style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>Cell Size (mm)</label>
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
                                        <div>
                                            <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                <label style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>Wall Thickness (mm)</label>
                                                <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '500' }}>{wallThickness}</span>
                                            </div>
                                            <input
                                                type="range"
                                                min="0.1"
                                                max="5"
                                                step="0.1"
                                                value={wallThickness}
                                                onChange={(e) => setWallThickness(parseFloat(e.target.value))}
                                                style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                            />
                                        </div>
                                    </>
                                )}

                                {textureType === 'decimate' && (
                                    <>
                                        <div>
                                            <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                                <label style={{ fontSize: '12px', color: 'var(--text-secondary)' }}>Reduction Intensity</label>
                                                <span style={{ fontSize: '12px', color: 'var(--accent-primary)', fontWeight: '500' }}>{Math.round(reductionRatio * 100)}%</span>
                                            </div>
                                            <input
                                                type="range"
                                                min="0.1"
                                                max="0.95"
                                                step="0.05"
                                                value={reductionRatio}
                                                onChange={(e) => setReductionRatio(parseFloat(e.target.value))}
                                                style={{ width: '100%', accentColor: 'var(--accent-primary)' }}
                                            />
                                            <p style={{ fontSize: '11px', color: 'var(--text-muted)', marginTop: '8px' }}>
                                                Higher intensity collapses more vertices into single points.
                                            </p>
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
                                    onClick={() => {
                                        if (textureType === 'knurling') {
                                            applyTexturize(selectedModel.id, { type: 'knurling', pitch, depth, angle, pattern: knurlPattern });
                                        } else if (textureType === 'honeycomb') {
                                            applyTexturize(selectedModel.id, { type: 'honeycomb', cellSize, wallThickness, depth });
                                        } else {
                                            applyTexturize(selectedModel.id, { type: 'decimate', reduction: reductionRatio });
                                        }
                                    }}
                                    style={{
                                        width: '100%',
                                        padding: '12px',
                                        backgroundColor: 'var(--accent-primary)',
                                        color: 'white',
                                        border: 'none',
                                        borderRadius: 'var(--radius-sm)',
                                        cursor: 'pointer',
                                        fontWeight: '600',
                                        fontSize: '13px',
                                        transition: 'all 0.2s ease',
                                        marginTop: '16px'
                                    }}
                                    onMouseEnter={(e) => e.currentTarget.style.filter = 'brightness(1.1)'}
                                    onMouseLeave={(e) => e.currentTarget.style.filter = 'none'}
                                >
                                    Apply Texture
                                </button>
                            </div>
                        </div>
                    )}
                </div>
            )}
        </div>
    );
};
