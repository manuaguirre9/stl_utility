import React from 'react';
import { RotateCcw } from 'lucide-react';
import { useStore } from '../../store/useStore';

export const ProcessingOverlay: React.FC = () => {
    const isProcessing = useStore((state) => state.isProcessing);
    const processingStatus = useStore((state) => state.processingStatus);

    if (!isProcessing) return null;

    return (
        <div style={{
            position: 'fixed',
            inset: 0,
            backgroundColor: 'rgba(0, 0, 0, 0.4)',
            backdropFilter: 'blur(4px)',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            zIndex: 9999,
            transition: 'all 0.3s ease'
        }}>
            <div style={{
                backgroundColor: 'var(--bg-panel)',
                padding: '24px 32px',
                borderRadius: '16px',
                border: '1px solid var(--border-color)',
                display: 'flex',
                flexDirection: 'column',
                alignItems: 'center',
                gap: '16px',
                boxShadow: '0 12px 48px rgba(0,0,0,0.5)',
                animation: 'pulse 2s infinite ease-in-out'
            }}>
                <div style={{
                    width: '48px',
                    height: '48px',
                    borderRadius: '50%',
                    backgroundColor: 'rgba(255, 107, 0, 0.1)',
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    color: 'var(--accent-primary)'
                }}>
                    <RotateCcw size={24} className="spin" />
                </div>

                <div style={{ textAlign: 'center' }}>
                    <h3 style={{
                        margin: 0,
                        fontSize: '16px',
                        fontWeight: '600',
                        color: 'var(--text-primary)',
                        letterSpacing: '0.02em'
                    }}>
                        {processingStatus || 'Computing Geometry'}
                    </h3>
                    <p style={{
                        margin: '4px 0 0 0',
                        fontSize: '12px',
                        color: 'var(--text-secondary)',
                        opacity: 0.8
                    }}>
                        Applying changes to your model...
                    </p>
                </div>

                {/* Small indicator bar */}
                <div style={{
                    width: '120px',
                    height: '3px',
                    backgroundColor: 'var(--bg-tertiary)',
                    borderRadius: '2px',
                    overflow: 'hidden',
                    position: 'relative'
                }}>
                    <div style={{
                        position: 'absolute',
                        height: '100%',
                        width: '40%',
                        backgroundColor: 'var(--accent-primary)',
                        borderRadius: '2px',
                        animation: 'loading-slide 1.5s infinite ease-in-out'
                    }} />
                </div>
            </div>

            <style>{`
                @keyframes loading-slide {
                    0% { left: -40%; }
                    100% { left: 100%; }
                }
                @keyframes pulse {
                    0% { transform: scale(1); }
                    50% { transform: scale(1.02); }
                    100% { transform: scale(1); }
                }
                .spin {
                    animation: spin 2s linear infinite;
                }
                @keyframes spin {
                    from { transform: rotate(0deg); }
                    to { transform: rotate(360deg); }
                }
            `}</style>
        </div>
    );
};
