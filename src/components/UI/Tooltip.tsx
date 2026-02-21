import React, { useState } from 'react';
import type { ReactNode } from 'react';

interface TooltipProps {
    content: string;
    children: ReactNode;
    icon?: ReactNode;
    delay?: number;
    align?: 'left' | 'center' | 'right';
}

export const Tooltip: React.FC<TooltipProps> = ({ content, children, icon, align = 'center' }) => {
    const [isVisible, setIsVisible] = useState(false);
    const [pos, setPos] = useState({ x: 0, y: 0 });

    const handleMouseMove = (e: React.MouseEvent) => {
        setPos({ x: e.clientX, y: e.clientY });
    };

    const getTransform = () => {
        if (align === 'right') return 'translateX(15px)';
        if (align === 'left') return 'translateX(calc(-100% - 15px))';
        return 'translateX(-50%)';
    };

    return (
        <div
            onMouseEnter={(e) => {
                setPos({ x: e.clientX, y: e.clientY });
                setIsVisible(true);
            }}
            onMouseMove={handleMouseMove}
            onMouseLeave={() => setIsVisible(false)}
            style={{ display: 'inline-flex', alignItems: 'center', justifyContent: 'center' }}
        >
            {children}
            {isVisible && (
                <div style={{
                    position: 'fixed',
                    left: `${Math.max(10, Math.min(window.innerWidth - 10, pos.x))}px`,
                    top: `${pos.y - 40}px`,
                    transform: getTransform(),
                    backgroundColor: 'rgba(0, 0, 0, 0.85)',
                    color: 'white',
                    padding: '6px 12px',
                    borderRadius: '6px',
                    fontSize: '12px',
                    fontWeight: '500',
                    pointerEvents: 'none',
                    zIndex: 10000,
                    whiteSpace: 'nowrap',
                    boxShadow: '0 4px 12px rgba(0,0,0,0.5)',
                    backdropFilter: 'blur(4px)',
                    border: '1px solid rgba(255,255,255,0.1)',
                    display: 'flex',
                    alignItems: 'center',
                    gap: '8px',
                    animation: 'fadeInUp 0.15s ease-out'
                }}>
                    {icon && <div style={{ color: '#ffcc00', display: 'flex' }}>{icon}</div>}
                    {content}
                </div>
            )}
        </div>
    );
};
