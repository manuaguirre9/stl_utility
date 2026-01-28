import React from 'react';

import { Menu } from 'lucide-react';

interface TopBarProps {
    onToggleSidebar: () => void;
}

export const TopBar: React.FC<TopBarProps> = ({ onToggleSidebar }) => {
    return (
        <div style={{
            height: 'var(--header-height)',
            backgroundColor: 'var(--bg-panel-header)',
            display: 'flex',
            alignItems: 'center',
            padding: '0 var(--spacing-md)',
            borderBottom: '1px solid var(--border-color)',
            color: 'var(--text-primary)',
            fontWeight: 'bold',
            fontSize: '14px',
            justifyContent: 'space-between'
        }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
                <button
                    className="mobile-only"
                    onClick={onToggleSidebar}
                    style={{ padding: '4px' }}
                >
                    <Menu size={20} />
                </button>
                <span style={{
                    whiteSpace: 'nowrap',
                    overflow: 'hidden',
                    textOverflow: 'ellipsis',
                    maxWidth: '180px'
                }}>
                    3D Design Editor
                </span>
            </div>

            <div style={{ fontSize: '11px', color: 'var(--text-secondary)' }} className="desktop-only">
                V2.0.0 Stable
            </div>
        </div>
    );
};
