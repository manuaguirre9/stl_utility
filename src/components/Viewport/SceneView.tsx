import React, { Suspense, useState } from 'react';
import { Tooltip } from '../UI/Tooltip';
import { Layers, Box } from 'lucide-react';
import { Canvas } from '@react-three/fiber';
import { OrbitControls, Grid, GizmoHelper, GizmoViewport, Bounds, Line } from '@react-three/drei';
import { useStore } from '../../store/useStore';
import * as THREE from 'three';
import { generateSelectionGeometry, segmentMesh, generateClassificationGeometry } from '../../utils/selectionUtils';
import { useIsMobile } from '../../utils/useIsMobile';
import { ChevronDown, ChevronUp, Info } from 'lucide-react';

const EMPTY_ARRAY: number[] = [];

const Model = ({ model }: { model: any }) => {
    const selectModel = useStore((state) => state.selectModel);
    const showMesh = useStore((state) => state.showMesh);
    const showEdges = useStore((state) => state.showEdges);
    const setSmartSelection = useStore((state) => state.setSmartSelection);
    const showClassification = useStore((state) => state.showClassification);
    const classificationAngle = useStore((state) => state.classificationAngle);
    const addSelection = useStore((state) => state.addToSmartSelection);
    const removeSelection = useStore((state) => state.removeFromSmartSelection);
    const selection = useStore((state) => state.smartSelection[model.id] || EMPTY_ARRAY);
    const boundaryEdges = useStore((state) => state.boundaryEdges[model.id]);
    const nonManifoldEdges = useStore((state) => state.nonManifoldEdges[model.id]);
    const fillableEdges = useStore((state) => state.fillableEdges[model.id]);
    const unfillableEdges = useStore((state) => state.unfillableEdges[model.id]);

    const historyPreviewId = useStore((state) => state.historyPreviewId);
    const selectedHistoryIds = useStore((state) => state.selectedHistoryIds);
    const history = useStore((state) => state.history);

    // Derived geometry for highlight
    const highlightGeo = React.useMemo(() => {
        if (selection.length === 0) return null;
        return generateSelectionGeometry(model.bufferGeometry, selection);
    }, [model.bufferGeometry, selection]);

    const previewGeo = React.useMemo(() => {
        // Hover takes priority, fallback to last selected history item
        const activeId = historyPreviewId || (selectedHistoryIds.length > 0 ? selectedHistoryIds[selectedHistoryIds.length - 1] : null);
        if (!activeId) return null;

        const entryIdx = history.findIndex(h => h.id === activeId);
        if (entryIdx <= 0) return null;

        const entry = history[entryIdx];
        if (!('selection' in entry.action) || (entry.action as any).modelId !== model.id) return null;

        // Use the state BEFORE the action was applied to show what was selected correctly
        const prevState = history[entryIdx - 1];
        const prevModel = prevState.models.find(m => m.id === model.id);
        if (!prevModel) return null;

        try {
            return {
                geometry: generateSelectionGeometry(prevModel.bufferGeometry, (entry.action as any).selection),
                model: prevModel
            };
        } catch (e) {
            return null;
        }
    }, [historyPreviewId, selectedHistoryIds, history, model.id]);

    // Derived segmentation data for both visualization and smart selection
    const segmentation = React.useMemo(() => {
        const segments = segmentMesh(model.bufferGeometry, classificationAngle);
        // Build a quick lookup map: faceIndex -> segment
        const faceToSegmentMap = new Map<number, number[]>();
        segments.forEach((seg: { indices: number[] }) => {
            seg.indices.forEach((idx: number) => {
                faceToSegmentMap.set(idx, seg.indices);
            });
        });
        return { segments, faceToSegmentMap };
    }, [model.bufferGeometry, model.meshVersion, classificationAngle]);

    const classificationGeo = React.useMemo(() => {
        if (!showClassification) return null;
        return generateClassificationGeometry(model.bufferGeometry, segmentation.segments, selection);
    }, [model.bufferGeometry, showClassification, segmentation.segments, selection]);

    const handlePointerDown = (e: any) => {
        if (e.button !== 0) return; // Only process selection on left click
        e.stopPropagation();
        selectModel(model.id);

        if (e.faceIndex !== undefined) {
            const faceIndex = e.faceIndex;
            if (e.shiftKey) {
                if (selection.includes(faceIndex)) {
                    removeSelection(model.id, [faceIndex]);
                } else {
                    addSelection(model.id, [faceIndex]);
                }
            } else {
                setSmartSelection(model.id, [faceIndex]);
            }
        }
    };

    const handleDoubleClick = (e: any) => {
        e.stopPropagation();
        if (e.faceIndex !== undefined) {
            // Smart select the entire feature (patch) the clicked face belongs to
            const connected = segmentation.faceToSegmentMap.get(e.faceIndex) || [e.faceIndex];

            if (e.shiftKey) {
                if (selection.includes(e.faceIndex)) {
                    removeSelection(model.id, connected);
                } else {
                    addSelection(model.id, connected);
                }
            } else {
                setSmartSelection(model.id, connected);
            }
        }
    };

    // Ensure bounding box/sphere are precomputed and valid (prevents Bounds crash)
    React.useEffect(() => {
        const geo = model.bufferGeometry;
        if (geo.attributes.position) {
            const arr = geo.attributes.position.array as Float32Array;
            let nan = 0;
            for (let i = 0; i < arr.length; i++) {
                if (!Number.isFinite(arr[i])) { arr[i] = 0; nan++; }
            }
            if (nan > 0) {
                console.warn(`[Model] Sanitized ${nan} NaN position values in model ${model.id}`);
                geo.attributes.position.needsUpdate = true;
            }
        }
        if (geo.attributes.normal) {
            const arr = geo.attributes.normal.array as Float32Array;
            for (let i = 0; i < arr.length; i++) {
                if (!Number.isFinite(arr[i])) arr[i] = 0;
            }
        }
        geo.computeBoundingBox();
        geo.computeBoundingSphere();
    }, [model.bufferGeometry, model.meshVersion, model.id]);

    // Safe edges geometry: EdgesGeometry can produce NaN on degenerate meshes
    const safeEdgesGeo = React.useMemo(() => {
        if (!showEdges) return null;
        try {
            const edgesGeo = new THREE.EdgesGeometry(model.bufferGeometry, 15);
            // Sanitize edge positions
            if (edgesGeo.attributes.position) {
                const arr = edgesGeo.attributes.position.array as Float32Array;
                for (let i = 0; i < arr.length; i++) {
                    if (!Number.isFinite(arr[i])) arr[i] = 0;
                }
            }
            edgesGeo.computeBoundingBox();
            edgesGeo.computeBoundingSphere();
            return edgesGeo;
        } catch (e) {
            console.warn('[Model] EdgesGeometry failed:', e);
            return null;
        }
    }, [model.bufferGeometry, model.meshVersion, showEdges]);

    // Boundary edges geometry for highlighting non-manifold/open edges
    const boundaryGeo = React.useMemo(() => {
        if (!boundaryEdges) return null;
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.BufferAttribute(boundaryEdges, 3));
        return geo;
    }, [boundaryEdges]);

    const nonManifoldGeo = React.useMemo(() => {
        if (!nonManifoldEdges) return null;
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.BufferAttribute(nonManifoldEdges, 3));
        return geo;
    }, [nonManifoldEdges]);

    // Live preview for mesh repair tools (watertight fill)
    const fillableGeo = React.useMemo(() => {
        if (!fillableEdges) return null;
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.BufferAttribute(fillableEdges, 3));
        return geo;
    }, [fillableEdges]);

    const unfillableGeo = React.useMemo(() => {
        if (!unfillableEdges) return null;
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.BufferAttribute(unfillableEdges, 3));
        return geo;
    }, [unfillableEdges]);

    return (
        <mesh
            geometry={model.bufferGeometry}
            position={model.position}
            rotation={model.rotation}
            scale={model.scale}
            visible={model.visible}
            onPointerDown={handlePointerDown}
            onDoubleClick={handleDoubleClick}
        >
            <meshStandardMaterial
                color={model.color}
                roughness={0.5}
                metalness={0.1}
                side={THREE.DoubleSide}
                polygonOffset
                polygonOffsetFactor={1}
                polygonOffsetUnits={1}
                visible={!showClassification}
            />
            {classificationGeo && (
                <mesh geometry={classificationGeo} renderOrder={0.5}>
                    <meshBasicMaterial
                        vertexColors
                        side={THREE.DoubleSide}
                        polygonOffset
                        polygonOffsetFactor={0.1}
                        polygonOffsetUnits={0.1}
                    />
                </mesh>
            )}
            {previewGeo && (
                <mesh geometry={previewGeo.geometry} renderOrder={0.8}>
                    <meshBasicMaterial
                        color="#ffcc00"
                        side={THREE.DoubleSide}
                        polygonOffset
                        polygonOffsetFactor={0.2}
                        polygonOffsetUnits={0.2}
                        transparent
                        opacity={0.8}
                    />
                </mesh>
            )}
            {highlightGeo && (
                <mesh geometry={highlightGeo} renderOrder={1}>
                    <meshBasicMaterial
                        color="#ffcc00"
                        side={THREE.DoubleSide}
                        polygonOffset
                        polygonOffsetFactor={0.5}
                        polygonOffsetUnits={0.5}
                    />
                </mesh>
            )}
            {showMesh && (
                <mesh
                    geometry={model.bufferGeometry}
                    key={`wire-${model.id}-${model.meshVersion}`}
                    renderOrder={2}
                >
                    <meshBasicMaterial
                        color="#000000"
                        wireframe
                        transparent
                        opacity={0.3}
                        depthTest={true}
                        polygonOffset
                        polygonOffsetFactor={-1}
                        polygonOffsetUnits={-1}
                    />
                </mesh>
            )}
            {safeEdgesGeo && (
                <lineSegments
                    key={`edges-${model.id}-${model.meshVersion}`}
                    geometry={safeEdgesGeo}
                    renderOrder={3}
                >
                    <lineBasicMaterial color="#000000" />
                </lineSegments>
            )}
            {boundaryGeo && !fillableGeo && !unfillableGeo && boundaryEdges && boundaryEdges.length >= 6 && (
                <Line
                    points={Array.from(boundaryEdges)}
                    color="#ff6b00"
                    lineWidth={5}
                    depthTest={false}
                    transparent
                />
            )}
            {nonManifoldGeo && nonManifoldEdges && nonManifoldEdges.length >= 6 && (
                <Line
                    points={Array.from(nonManifoldEdges)}
                    color="#cc00ff"
                    lineWidth={8}
                    depthTest={false}
                    transparent
                />
            )}
            {fillableGeo && fillableEdges && fillableEdges.length >= 6 && (
                <Line
                    points={Array.from(fillableEdges)}
                    color="#ddff00"
                    lineWidth={6}
                    depthTest={false}
                    transparent
                />
            )}
            {unfillableGeo && unfillableEdges && unfillableEdges.length >= 6 && (
                <Line
                    points={Array.from(unfillableEdges)}
                    color="#ff3300"
                    lineWidth={8}
                    depthTest={false}
                    transparent
                />
            )}
        </mesh>
    );
};

export const SceneView: React.FC = () => {
    const isMobile = useIsMobile();
    const models = useStore((state) => state.models);
    const selectModel = useStore((state) => state.selectModel);
    const showMesh = useStore((state) => state.showMesh);
    const toggleShowMesh = useStore((state) => state.toggleShowMesh);
    const showEdges = useStore((state) => state.showEdges);
    const toggleShowEdges = useStore((state) => state.toggleShowEdges);
    const showClassification = useStore((state) => state.showClassification);
    const toggleClassification = useStore((state) => state.toggleClassification);
    const classificationAngle = useStore((state) => state.classificationAngle);
    const setClassificationAngle = useStore((state) => state.setClassificationAngle);
    const transformMode = useStore((state) => state.transformMode);
    const clearSmartSelection = useStore((state) => state.clearSmartSelection);
    const setSelectedHistoryIds = useStore((state) => state.setSelectedHistoryIds);
    const orbitControlsRef = React.useRef<any>(null);

    const [isLegendExpanded, setIsLegendExpanded] = useState(!isMobile);

    // ── Global NaN guard: auto-fix any geometry with NaN positions ──
    React.useEffect(() => {
        const original = THREE.BufferGeometry.prototype.computeBoundingBox;
        THREE.BufferGeometry.prototype.computeBoundingBox = function (this: THREE.BufferGeometry) {
            original.call(this);
            if (this.boundingBox) {
                const { min, max } = this.boundingBox;
                if ([min.x, min.y, min.z, max.x, max.y, max.z].some(v => !Number.isFinite(v))) {
                    const pos = this.attributes.position;
                    if (pos) {
                        const arr = pos.array as Float32Array;
                        let fixed = 0;
                        for (let i = 0; i < arr.length; i++) {
                            if (!Number.isFinite(arr[i])) { arr[i] = 0; fixed++; }
                        }
                        if (fixed > 0) {
                            console.warn(`[NaN Guard] Fixed ${fixed} NaN values in geometry ${this.uuid} (posCount=${pos.count}, arrayLen=${arr.length})`);
                            pos.needsUpdate = true;
                            // Recompute with clean data
                            original.call(this);
                        }
                    }
                }
            }
        };
        return () => { THREE.BufferGeometry.prototype.computeBoundingBox = original; };
    }, []);

    React.useEffect(() => {
        const handleKeyDown = (e: KeyboardEvent) => {
            if (e.key === 'Escape') {
                clearSmartSelection();
                selectModel(null);
            }
        };
        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [transformMode, clearSmartSelection, selectModel]);

    const setView = (view: string) => {
        if (!orbitControlsRef.current || models.length === 0) return;

        // Calculate bounding box of all visible models
        const box = new THREE.Box3();
        let hasVisible = false;

        models.forEach(model => {
            if (model.visible) {
                hasVisible = true;
                if (!model.bufferGeometry.boundingBox) model.bufferGeometry.computeBoundingBox();

                // Construct basic matrix from store props
                const matrix = new THREE.Matrix4();
                const euler = new THREE.Euler(...model.rotation);
                const quaternion = new THREE.Quaternion().setFromEuler(euler);
                matrix.compose(
                    new THREE.Vector3(...model.position),
                    quaternion,
                    new THREE.Vector3(...model.scale)
                );

                // transform local bounds to world
                const modelBox = model.bufferGeometry.boundingBox!.clone();
                modelBox.applyMatrix4(matrix);
                box.union(modelBox);
            }
        });

        if (!hasVisible) return;

        const center = new THREE.Vector3();
        box.getCenter(center);
        const size = new THREE.Vector3();
        box.getSize(size);

        const maxDim = Math.max(size.x, size.y, size.z);
        const fov = 50;
        const distance = maxDim / (2 * Math.tan(fov * Math.PI / 360));
        const offset = distance * 1.2;

        const camera = orbitControlsRef.current.object;

        switch (view) {
            case 'Front':
                camera.position.set(center.x, center.y, center.z + offset);
                break;
            case 'Top':
                camera.position.set(center.x, center.y + offset, center.z);
                break;
            case 'Right':
                camera.position.set(center.x + offset, center.y, center.z);
                break;
            case 'Perspective':
                camera.position.set(center.x + offset, center.y + offset, center.z + offset);
                break;
        }

        orbitControlsRef.current.target.copy(center);
        camera.lookAt(center);
        orbitControlsRef.current.update();
    };

    return (
        <div style={{ flex: 1, backgroundColor: '#121212', position: 'relative', overflow: 'hidden' }}>
            <Canvas
                camera={{ position: [50, 50, 50], fov: 50 }}
                onPointerMissed={(e) => {
                    if (e.button !== 0) return; // Only clear on left click
                    selectModel(null);
                    clearSmartSelection();
                    setSelectedHistoryIds(['initial']);
                }}
            >
                <ambientLight intensity={0.5} />
                <pointLight position={[10, 10, 10]} />
                <directionalLight position={[-10, 20, 10]} intensity={1} castShadow />

                <Grid
                    infiniteGrid
                    cellSize={10}
                    sectionSize={50}
                    fadeDistance={200}
                    sectionColor="#505050"
                    cellColor="#303030"
                />

                <Bounds fit clip observe margin={1.2}>
                    <Suspense fallback={null}>
                        {models.map((model) => (
                            <Model key={model.id} model={model} />
                        ))}
                    </Suspense>
                </Bounds>

                <OrbitControls
                    ref={orbitControlsRef}
                    makeDefault
                    screenSpacePanning
                    zoomSpeed={1.5}
                    panSpeed={1.5}
                />
                {!isMobile && (
                    <GizmoHelper alignment="bottom-right" margin={[80, 80]}>
                        <GizmoViewport axisColors={['#ff3b30', '#4cd964', '#007aff']} labelColor="black" />
                    </GizmoHelper>
                )}
            </Canvas>

            {/* Classification Legend */}
            {showClassification && (
                <div style={{
                    position: 'absolute',
                    top: isMobile ? '8px' : '20px',
                    left: isMobile ? '8px' : '20px',
                    backgroundColor: 'rgba(25, 25, 25, 0.9)',
                    backdropFilter: 'blur(12px)',
                    padding: isLegendExpanded ? (isMobile ? '10px' : '12px 16px') : '8px',
                    borderRadius: 'var(--radius-md)',
                    border: '1px solid var(--border-color)',
                    color: '#fff',
                    display: 'flex',
                    flexDirection: 'column',
                    gap: isLegendExpanded ? '10px' : '0',
                    zIndex: 200,
                    boxShadow: '0 8px 32px rgba(0,0,0,0.4)',
                    width: isLegendExpanded ? (isMobile ? '160px' : '200px') : 'auto',
                    transition: 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)'
                }}>
                    <div
                        onClick={() => setIsLegendExpanded(!isLegendExpanded)}
                        style={{
                            display: 'flex',
                            justifyContent: 'space-between',
                            alignItems: 'center',
                            cursor: 'pointer',
                            pointerEvents: 'auto'
                        }}
                    >
                        <div style={{
                            display: 'flex',
                            alignItems: 'center',
                            gap: '6px',
                            fontSize: '10px',
                            fontWeight: 'bold',
                            textTransform: 'uppercase',
                            color: 'var(--text-secondary)',
                            letterSpacing: '0.05em'
                        }}>
                            <Info size={12} />
                            {isLegendExpanded ? 'Analysis' : ''}
                        </div>
                        {isLegendExpanded ? <ChevronUp size={14} /> : <ChevronDown size={14} />}
                    </div>

                    {isLegendExpanded && (
                        <>
                            {[
                                { label: 'Plane', color: '#3498db' },
                                { label: 'Cylinder', color: '#2ecc71' },
                                { label: 'Cone', color: '#9b59b6' },
                                { label: 'Generic', color: '#95a5a6' }
                            ].map(item => (
                                <div key={item.label} style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                                    <div style={{ width: '10px', height: '10px', borderRadius: '3px', backgroundColor: item.color }} />
                                    <span style={{ fontSize: '12px', color: 'var(--text-primary)' }}>{item.label}</span>
                                </div>
                            ))}

                            <div style={{
                                marginTop: '4px',
                                paddingTop: '8px',
                                borderTop: '1px solid rgba(255,255,255,0.1)',
                                display: 'flex',
                                flexDirection: 'column',
                                gap: '6px',
                                pointerEvents: 'auto'
                            }}>
                                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                                    <span style={{ fontSize: '10px', color: 'var(--text-secondary)' }}>Angle</span>
                                    <span style={{ fontSize: '10px', color: 'var(--accent-primary)', fontWeight: 'bold' }}>{classificationAngle.toFixed(1)}°</span>
                                </div>
                                <input
                                    type="range"
                                    min="0.2"
                                    max="45"
                                    step="0.2"
                                    value={classificationAngle}
                                    onChange={(e) => setClassificationAngle(parseFloat(e.target.value))}
                                    style={{
                                        width: '100%',
                                        cursor: 'pointer',
                                        accentColor: 'var(--accent-primary)',
                                        height: '4px'
                                    }}
                                />
                            </div>
                        </>
                    )}
                </div>
            )}

            {transformMode === 'stitching' && (
                <div style={{
                    position: 'absolute',
                    bottom: isMobile ? '60px' : '80px',
                    right: '20px',
                    backgroundColor: 'rgba(20, 20, 20, 0.85)',
                    backdropFilter: 'blur(8px)',
                    padding: '12px',
                    borderRadius: '8px',
                    border: '1px solid var(--border-color)',
                    zIndex: 100,
                    display: 'flex',
                    flexDirection: 'column',
                    gap: '8px',
                    fontSize: '11px',
                    color: 'var(--text-secondary)'
                }}>
                    <strong style={{ color: '#fff', marginBottom: '4px' }}>Repair Diagnostics</strong>
                    <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                        <div style={{ width: '12px', height: '4px', backgroundColor: '#00ff00', borderRadius: '2px' }} />
                        <span>Ready to Repair (Safe)</span>
                    </div>
                    <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                        <div style={{ width: '12px', height: '4px', backgroundColor: '#ff6b00', borderRadius: '2px' }} />
                        <span>Raw Open Boundary</span>
                    </div>
                    <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                        <div style={{ width: '12px', height: '4px', backgroundColor: '#cc00ff', borderRadius: '2px' }} />
                        <span>Non-Manifold Geometry</span>
                    </div>
                    <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                        <div style={{ width: '12px', height: '4px', backgroundColor: '#ff3300', borderRadius: '2px' }} />
                        <span>Too Large / Complex</span>
                    </div>
                </div>
            )}

            {/* Bottom View Controls Overlay */}
            <div style={{
                position: 'absolute',
                bottom: isMobile ? '8px' : '20px',
                left: '50%',
                transform: 'translateX(-50%)',
                display: 'flex',
                gap: isMobile ? '4px' : '8px',
                backgroundColor: 'rgba(30, 30, 30, 0.85)',
                backdropFilter: 'blur(8px)',
                padding: isMobile ? '4px' : '8px',
                borderRadius: 'var(--radius-md)',
                border: '1px solid var(--border-color)',
                zIndex: 100,
                maxWidth: '95vw',
                overflowX: 'auto',
                scrollbarWidth: 'none'
            }}>

                {(isMobile ? ['Front', 'Top', 'Right'] : ['Perspective', 'Front', 'Top', 'Right']).map(view => (
                    <Tooltip key={view} content={`${view} View`}>
                        <button
                            onClick={() => setView(view)}
                            style={{
                                padding: isMobile ? '4px 8px' : '4px 12px',
                                borderRadius: 'var(--radius-sm)',
                                fontSize: isMobile ? '10px' : '12px',
                                backgroundColor: 'transparent',
                                color: 'var(--text-primary)',
                                border: '1px solid var(--border-color)',
                                whiteSpace: 'nowrap'
                            }}
                        >
                            {view}
                        </button>
                    </Tooltip>
                ))}
                <div style={{ width: '1px', backgroundColor: 'var(--border-color)', margin: '0 4px' }} />

                <Tooltip content="Toggle Wireframe Overlay" icon={<Box size={12} />}>
                    <button
                        onClick={toggleShowMesh}
                        style={{
                            padding: isMobile ? '4px 8px' : '4px 12px',
                            borderRadius: 'var(--radius-sm)',
                            fontSize: isMobile ? '10px' : '12px',
                            border: '1px solid var(--border-color)',
                            color: showMesh ? '#fff' : 'var(--text-primary)',
                            backgroundColor: showMesh ? 'var(--accent-primary)' : 'transparent',
                        }}
                    >
                        Mesh
                    </button>
                </Tooltip>

                <Tooltip content="Toggle Feature Edges" icon={<Layers size={12} />}>
                    <button
                        onClick={toggleShowEdges}
                        style={{
                            padding: isMobile ? '4px 8px' : '4px 12px',
                            borderRadius: 'var(--radius-sm)',
                            fontSize: isMobile ? '10px' : '12px',
                            border: '1px solid var(--border-color)',
                            color: showEdges ? '#fff' : 'var(--text-primary)',
                            backgroundColor: showEdges ? 'var(--accent-primary)' : 'transparent',
                        }}
                    >
                        Edges
                    </button>
                </Tooltip>

                <Tooltip content="Toggle Surface Analysis (Type Detection)" icon={<Info size={12} />}>
                    <button
                        onClick={toggleClassification}
                        style={{
                            padding: isMobile ? '4px 8px' : '4px 12px',
                            borderRadius: 'var(--radius-sm)',
                            fontSize: isMobile ? '10px' : '12px',
                            border: '1px solid var(--border-color)',
                            color: showClassification ? '#fff' : 'var(--text-primary)',
                            backgroundColor: showClassification ? 'var(--accent-primary)' : 'transparent',
                        }}
                    >
                        Analysis
                    </button>
                </Tooltip>
            </div>
        </div >
    );
};

