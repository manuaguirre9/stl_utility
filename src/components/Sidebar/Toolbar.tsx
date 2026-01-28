import React, { useRef } from 'react';
import { useStore } from '../../store/useStore';
import { Trash2, Upload, Download, MousePointerClick, Scissors, Grid3X3 } from 'lucide-react';
import { exportSceneToSTL } from '../../utils/exportUtils';
import { v4 as uuidv4 } from 'uuid';
import { STLLoader } from 'three/examples/jsm/loaders/STLLoader.js';
import { useIsMobile } from '../../utils/useIsMobile';

export const Toolbar: React.FC = () => {
    const isMobile = useIsMobile();
    const fileInputRef = useRef<HTMLInputElement>(null);
    const addModel = useStore((state) => state.addModel);
    const transformMode = useStore((state) => state.transformMode);
    const setTransformMode = useStore((state) => state.setTransformMode);
    const selectedId = useStore((state) => state.selectedId);
    const models = useStore((state) => state.models);
    const removeModel = useStore((state) => state.removeModel);

    const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
        const file = event.target.files?.[0];
        if (!file) return;

        const reader = new FileReader();
        reader.onload = (e) => {
            const contents = e.target?.result as ArrayBuffer;
            const loader = new STLLoader();
            const geometry = loader.parse(contents);

            // Validate geometry
            geometry.computeBoundingSphere();
            geometry.computeBoundingBox();

            const sphere = geometry.boundingSphere;
            const box = geometry.boundingBox;

            if (!sphere || !box || isNaN(sphere.radius) || !isFinite(sphere.radius) || sphere.radius <= 0) {
                alert("Invalid or corrupted STL file. Geometry bounds are infinite or zero.");
                return;
            }

            // Center geometry specifically
            geometry.center();

            addModel({
                id: uuidv4(),
                name: file.name,
                bufferGeometry: geometry,
                color: '#ff6b00', // Default orange
                visible: true,
                position: [0, 0, 0],
                rotation: [0, 0, 0],
                scale: [1, 1, 1],
                meshVersion: 0,
            });
        };
        reader.readAsArrayBuffer(file);
        // Reset input
        event.target.value = '';
    };

    const ToolButton = ({
        active,
        onClick,
        icon: Icon,
        primary = false
    }: {
        active?: boolean;
        onClick: () => void;
        icon: React.ElementType<{ size?: number | string }>;
        primary?: boolean
    }) => (
        <button
            onClick={onClick}
            style={{
                width: isMobile ? '36px' : '40px',
                height: isMobile ? '36px' : '40px',
                borderRadius: 'var(--radius-md)',
                backgroundColor: primary ? 'var(--accent-primary)' : (active ? 'var(--bg-hover)' : 'transparent'),
                color: primary ? '#fff' : (active ? '#fff' : 'var(--text-secondary)'),
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                marginBottom: isMobile ? '0' : '8px',
                marginRight: isMobile ? '8px' : '0',
                transition: 'all 0.2s ease',
            }}
        >
            <Icon size={isMobile ? 18 : 20} />
        </button>
    );

    return (
        <div style={{
            width: isMobile ? 'auto' : 'var(--toolbar-width)',
            height: isMobile ? 'auto' : '100%',
            backgroundColor: 'var(--bg-panel)',
            borderRight: isMobile ? 'none' : '1px solid var(--border-color)',
            borderBottom: isMobile ? '1px solid var(--border-color)' : 'none',
            display: 'flex',
            flexDirection: isMobile ? 'row' : 'column',
            alignItems: 'center',
            padding: isMobile ? 'var(--spacing-xs) var(--spacing-sm)' : 'var(--spacing-sm) 0',
            zIndex: 10,
            overflowX: isMobile ? 'auto' : 'visible'
        }}>
            <input
                type="file"
                ref={fileInputRef}
                onChange={handleFileUpload}
                accept=".stl"
                style={{ display: 'none' }}
            />

            <ToolButton
                primary
                icon={Upload}
                onClick={() => fileInputRef.current?.click()}
            />

            <ToolButton
                primary
                icon={Download}
                onClick={() => exportSceneToSTL(models)}
            />

            {!isMobile && <div style={{ height: '10px' }} />}
            {isMobile && <div style={{ width: '10px' }} />}

            <ToolButton
                active={transformMode === 'smart'}
                icon={MousePointerClick}
                onClick={() => setTransformMode('smart')}
            />
            <ToolButton
                active={transformMode === 'subdivide'}
                icon={Scissors}
                onClick={() => setTransformMode('subdivide')}
            />
            <ToolButton
                active={transformMode === 'texturize'}
                icon={Grid3X3}
                onClick={() => setTransformMode('texturize')}
            />

            {!isMobile && <div style={{ height: '10px', marginTop: 'auto' }} />}

            <ToolButton
                icon={Trash2}
                onClick={() => selectedId && removeModel(selectedId)}
            />
        </div>
    );
};

