import React, { useState, useEffect, useRef } from 'react';
import { useStore } from '../../store/useStore';
import { useIsMobile } from '../../utils/useIsMobile';
import {
    Play,
    Pause,
    SkipBack,
    SkipForward,
    ChevronLeft,
    ChevronRight,
    Scissors,
    Grid3X3,
    FilePlus,
    FileMinus,
    Settings2,
    Home,
    Undo2,
    Redo2,
    Trash2
} from 'lucide-react';

const getIconForAction = (label: string) => {
    const l = label.toLowerCase();
    if (l.includes('initial') || l.includes('start')) return <Home size={14} />;
    if (l.includes('subdivide')) return <Scissors size={14} />;
    if (l.includes('texturize')) return <Grid3X3 size={14} />;
    if (l.includes('import')) return <FilePlus size={14} />;
    if (l.includes('remove')) return <FileMinus size={14} />;
    if (l.includes('edit') || l.includes('update')) return <Settings2 size={14} />;
    return <Settings2 size={14} />;
};

export const HistoryTimeline: React.FC = () => {
    const isMobile = useIsMobile();
    const history = useStore((state) => state.history);
    const historyIndex = useStore((state) => state.historyIndex);
    const selectedHistoryIds = useStore((state) => state.selectedHistoryIds);
    const setSelectedHistoryIds = useStore((state) => state.setSelectedHistoryIds);
    const deleteHistoryEntry = useStore((state) => state.deleteHistoryEntry);
    const jumpToHistory = useStore((state) => state.jumpToHistory);
    const undo = useStore((state) => state.undo);
    const redo = useStore((state) => state.redo);
    const historyPreviewId = useStore((state) => state.historyPreviewId);
    const setHistoryPreviewId = useStore((state) => state.setHistoryPreviewId);
    const reEditHistoryItem = useStore((state) => state.reEditHistoryItem);

    const [isPlaying, setIsPlaying] = useState(false);
    const playTimerRef = useRef<any>(null);

    useEffect(() => {
        if (isPlaying) {
            playTimerRef.current = setInterval(() => {
                const nextIndex = historyIndex + 1;
                if (nextIndex < history.length) {
                    jumpToHistory(nextIndex);
                } else {
                    setIsPlaying(false);
                }
            }, 800);
        } else {
            if (playTimerRef.current) clearInterval(playTimerRef.current);
        }
        return () => {
            if (playTimerRef.current) clearInterval(playTimerRef.current);
        };
    }, [isPlaying, historyIndex, history.length, jumpToHistory]);

    // Handle Delete key
    useEffect(() => {
        const handleKeyDown = (e: KeyboardEvent) => {
            if (e.key === 'Delete' || e.key === 'Backspace') {
                if (selectedHistoryIds.length > 0) {
                    deleteHistoryEntry(selectedHistoryIds);
                }
            }
        };
        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [selectedHistoryIds, deleteHistoryEntry]);

    const handleItemClick = (index: number, entryId: string, event: React.MouseEvent) => {
        stopPlayback();
        jumpToHistory(index);

        if (event.shiftKey && selectedHistoryIds.length > 0) {
            // Find current selection indices
            const selectedIndices = history
                .map((h, i) => selectedHistoryIds.includes(h.id) ? i : -1)
                .filter(i => i !== -1);

            const minSelected = Math.min(...selectedIndices);
            const maxSelected = Math.max(...selectedIndices);

            const start = Math.min(minSelected, index);
            const end = Math.max(maxSelected, index);

            const newSelection = history.slice(start, end + 1).map(h => h.id);
            setSelectedHistoryIds(newSelection);
        } else {
            setSelectedHistoryIds([entryId]);
        }
    };

    const stopPlayback = () => setIsPlaying(false);

    return (
        <div style={{
            height: 'var(--timeline-height)',
            backgroundColor: 'var(--bg-panel)',
            borderTop: '2px solid var(--border-color)',
            display: 'flex',
            alignItems: 'center',
            padding: isMobile ? '0 5px' : '0 10px',
            zIndex: 10,
            gap: isMobile ? '8px' : '15px'
        }}>
            {/* Playback Controls */}
            <div style={{
                display: 'flex',
                alignItems: 'center',
                gap: '2px',
                paddingRight: isMobile ? '8px' : '15px',
                borderRight: '1px solid var(--border-color)',
                height: '70%'
            }}>
                {!isMobile && (
                    <button
                        onClick={() => { stopPlayback(); jumpToHistory(0); }}
                        style={controlButtonStyle}
                        title="Rewind to beginning"
                    >
                        <SkipBack size={16} />
                    </button>
                )}
                <button
                    onClick={() => { stopPlayback(); jumpToHistory(Math.max(0, historyIndex - 1)); }}
                    disabled={historyIndex === 0}
                    style={controlButtonStyle}
                    title="Step backward"
                >
                    <ChevronLeft size={isMobile ? 22 : 18} />
                </button>
                <button
                    onClick={() => setIsPlaying(!isPlaying)}
                    style={{
                        ...controlButtonStyle,
                        backgroundColor: isPlaying ? 'var(--accent-primary)' : 'var(--bg-input)',
                        color: isPlaying ? 'white' : 'var(--text-primary)',
                        padding: isMobile ? '8px' : '4px',
                        minWidth: isMobile ? '36px' : '28px'
                    }}
                    title={isPlaying ? "Pause" : "Play"}
                >
                    {isPlaying ? <Pause size={isMobile ? 20 : 18} /> : <Play size={isMobile ? 20 : 18} />}
                </button>
                <button
                    onClick={() => { stopPlayback(); redo(); }}
                    disabled={historyIndex === history.length - 1}
                    style={controlButtonStyle}
                    title="Step forward"
                >
                    <ChevronRight size={isMobile ? 22 : 18} />
                </button>
                {!isMobile && (
                    <button
                        onClick={() => { stopPlayback(); jumpToHistory(history.length - 1); }}
                        style={controlButtonStyle}
                        title="Fast forward to end"
                    >
                        <SkipForward size={16} />
                    </button>
                )}
            </div>

            {/* Timeline Icons */}
            <div style={{
                flex: 1,
                overflowX: 'auto',
                display: 'flex',
                alignItems: 'center',
                gap: '6px',
                height: '100%',
                scrollbarWidth: 'none',
                WebkitOverflowScrolling: 'touch'
            }}>
                {history.map((entry, index) => {
                    const isActive = index === historyIndex;
                    const isSelected = selectedHistoryIds.includes(entry.id);
                    const isFuture = index > historyIndex;
                    const isPreviewed = historyPreviewId === entry.id;

                    return (
                        <div
                            key={entry.id}
                            onClick={(e) => handleItemClick(index, entry.id, e)}
                            onDoubleClick={() => reEditHistoryItem(index)}
                            onMouseEnter={() => setHistoryPreviewId(entry.id)}
                            onMouseLeave={() => setHistoryPreviewId(null)}
                            title={`${entry.label} (Double-click to re-edit)`}
                            style={{
                                minWidth: isMobile ? '36px' : '32px',
                                height: isMobile ? '36px' : '32px',
                                display: 'flex',
                                alignItems: 'center',
                                justifyContent: 'center',
                                backgroundColor: isSelected ? 'var(--accent-primary)' : (isPreviewed ? 'rgba(0, 210, 255, 0.3)' : 'var(--bg-input)'),
                                borderRadius: 'var(--radius-sm)',
                                cursor: 'pointer',
                                color: isSelected ? 'white' : (isFuture ? 'var(--text-muted)' : 'var(--text-primary)'),
                                border: isSelected ? (isActive ? '2px solid white' : '2px solid rgba(255,255,255,0.5)') : (isPreviewed ? '1px solid #00d2ff' : '1px solid var(--border-color)'),
                                opacity: isFuture ? 0.5 : 1,
                                transition: 'all 0.1s ease',
                                flexShrink: 0,
                                transform: (isSelected || isPreviewed) ? 'scale(1.05)' : 'none',
                                boxShadow: isSelected ? '0 0 12px rgba(255,107,0,0.4)' : (isPreviewed ? '0 0 12px rgba(0, 210, 255, 0.4)' : 'none'),
                                zIndex: (isSelected || isPreviewed) ? 5 : 1
                            }}
                        >
                            {getIconForAction(entry.label)}
                        </div>
                    );
                })}
            </div>

            {/* Bulk Actions & Navigation */}
            <div style={{
                display: 'flex',
                gap: '8px',
                paddingLeft: isMobile ? '8px' : '15px',
                borderLeft: '1px solid var(--border-color)',
                height: '70%',
                alignItems: 'center'
            }}>
                {selectedHistoryIds.length > 1 && (
                    <button
                        onClick={() => deleteHistoryEntry(selectedHistoryIds)}
                        style={{
                            ...iconOnlyBtnStyle,
                            backgroundColor: '#ff4444',
                            color: 'white',
                            borderColor: '#ff4444',
                            padding: '4px 8px',
                            fontSize: '11px',
                            fontWeight: 'bold',
                            gap: '4px'
                        }}
                    >
                        <Trash2 size={12} />
                        Delete {selectedHistoryIds.length}
                    </button>
                )}
                {!isMobile && (
                    <>
                        <button onClick={undo} disabled={historyIndex === 0} style={iconOnlyBtnStyle}>
                            <Undo2 size={14} />
                        </button>
                        <button onClick={redo} disabled={historyIndex === history.length - 1} style={iconOnlyBtnStyle}>
                            <Redo2 size={14} />
                        </button>
                    </>
                )}
            </div>
        </div>
    );
};

const controlButtonStyle: React.CSSProperties = {
    padding: '4px',
    backgroundColor: 'transparent',
    border: 'none',
    borderRadius: '4px',
    color: 'var(--text-primary)',
    cursor: 'pointer',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    transition: 'background-color 0.2s',
    minWidth: '28px'
};

const iconOnlyBtnStyle: React.CSSProperties = {
    padding: '4px',
    backgroundColor: 'transparent',
    border: '1px solid var(--border-color)',
    borderRadius: '3px',
    color: 'var(--text-muted)',
    cursor: 'pointer',
    display: 'flex',
    alignItems: 'center'
};

