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
import { Tooltip } from '../UI/Tooltip';

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

    const stopPlayback = () => setIsPlaying(false);

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

    const handleItemSelection = (index: number, entryId: string, event: React.MouseEvent) => {
        event.stopPropagation();
        if (event.shiftKey && selectedHistoryIds.length > 0) {
            const selectedIndices = history
                .map((h, i) => selectedHistoryIds.includes(h.id) ? i : -1)
                .filter(i => i !== -1);
            const start = Math.min(Math.min(...selectedIndices), index);
            const end = Math.max(Math.max(...selectedIndices), index);
            setSelectedHistoryIds(history.slice(start, end + 1).map(h => h.id));
        } else {
            setSelectedHistoryIds([entryId]);
        }
    };

    const [isDragging, setIsDragging] = useState(false);
    const timelineContainerRef = useRef<HTMLDivElement>(null);

    const itemWidth = isMobile ? 32 : 28;
    const gap = 8;
    const padding = 4;

    const handleScrub = (clientX: number) => {
        if (!timelineContainerRef.current) return;
        const rect = timelineContainerRef.current.getBoundingClientRect();
        const x = clientX - rect.left + timelineContainerRef.current.scrollLeft;

        let targetIndex = 0;
        let minDistance = Infinity;

        for (let i = 0; i < history.length; i++) {
            const midGap = padding + (i * (itemWidth + gap)) + itemWidth + gap / 2;
            const dist = Math.abs(x - midGap);
            if (dist < minDistance) {
                minDistance = dist;
                targetIndex = i;
            }
        }

        if (targetIndex !== historyIndex && targetIndex >= 0 && targetIndex < history.length) {
            stopPlayback();
            jumpToHistory(targetIndex);
        }
    };

    const handleContainerMouseDown = () => {
        setIsDragging(true);
    };

    useEffect(() => {
        if (!isDragging) return;
        const onMouseMove = (e: MouseEvent) => handleScrub(e.clientX);
        const onMouseUp = () => setIsDragging(false);
        window.addEventListener('mousemove', onMouseMove);
        window.addEventListener('mouseup', onMouseUp);
        return () => {
            window.removeEventListener('mousemove', onMouseMove);
            window.removeEventListener('mouseup', onMouseUp);
        };
    }, [isDragging, history.length, historyIndex]);

    const scrubberPos = padding + (historyIndex * (itemWidth + gap)) + itemWidth + gap / 2;

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
                height: '70%',
                flexShrink: 0
            }}>
                <Tooltip content="Rewind" icon={<SkipBack size={12} />}>
                    <button onClick={() => { stopPlayback(); jumpToHistory(0); }} style={controlButtonStyle}>
                        <SkipBack size={16} />
                    </button>
                </Tooltip>
                <Tooltip content="Previous" icon={<ChevronLeft size={12} />}>
                    <button onClick={() => { stopPlayback(); jumpToHistory(Math.max(0, historyIndex - 1)); }} disabled={historyIndex === 0} style={controlButtonStyle}>
                        <ChevronLeft size={18} />
                    </button>
                </Tooltip>
                <Tooltip content={isPlaying ? "Pause" : "Play"} icon={isPlaying ? <Pause size={12} /> : <Play size={12} />}>
                    <button onClick={() => setIsPlaying(!isPlaying)} style={{ ...controlButtonStyle, backgroundColor: isPlaying ? 'var(--accent-primary)' : 'var(--bg-input)', color: isPlaying ? 'white' : 'var(--text-primary)' }}>
                        {isPlaying ? <Pause size={18} /> : <Play size={18} />}
                    </button>
                </Tooltip>
                <Tooltip content="Next" icon={<ChevronRight size={12} />}>
                    <button onClick={() => { stopPlayback(); redo(); }} disabled={historyIndex === history.length - 1} style={controlButtonStyle}>
                        <ChevronRight size={18} />
                    </button>
                </Tooltip>
                <Tooltip content="Fast Forward" icon={<SkipForward size={12} />}>
                    <button onClick={() => { stopPlayback(); jumpToHistory(history.length - 1); }} style={controlButtonStyle}>
                        <SkipForward size={16} />
                    </button>
                </Tooltip>
            </div>

            {/* Timeline Area (Drag to Scrub) */}
            <div
                ref={timelineContainerRef}
                onMouseDown={handleContainerMouseDown}
                style={{
                    flex: 1,
                    height: '100%',
                    display: 'flex',
                    alignItems: 'center',
                    gap: `${gap}px`,
                    position: 'relative',
                    overflowX: 'auto',
                    scrollbarWidth: 'none',
                    padding: `0 ${padding}px`,
                    cursor: 'ew-resize',
                    userSelect: 'none'
                }}
            >
                {history.map((entry, index) => {
                    const isActive = index === historyIndex;
                    const isSelected = selectedHistoryIds.includes(entry.id);
                    const isFuture = index > historyIndex;
                    const isPreviewed = historyPreviewId === entry.id;

                    return (
                        <Tooltip key={entry.id} content={entry.label} icon={getIconForAction(entry.label)}>
                            <div
                                onMouseEnter={() => setHistoryPreviewId(entry.id)}
                                onMouseLeave={() => setHistoryPreviewId(null)}
                                onDoubleClick={(e) => { e.stopPropagation(); reEditHistoryItem(index); }}
                                onClick={(e) => handleItemSelection(index, entry.id, e)}
                                style={{
                                    minWidth: `${itemWidth}px`,
                                    height: `${itemWidth}px`,
                                    display: 'flex',
                                    alignItems: 'center',
                                    justifyContent: 'center',
                                    backgroundColor: isSelected ? 'var(--accent-primary)' : (isPreviewed ? 'rgba(255, 204, 0, 0.3)' : 'var(--bg-input)'),
                                    borderRadius: 'var(--radius-sm)',
                                    color: isSelected ? 'white' : (isFuture ? 'var(--text-muted)' : 'var(--text-primary)'),
                                    border: isSelected ? (isActive ? '2px solid white' : '2px solid rgba(255,255,255,0.5)') : (isPreviewed ? '1px solid #ffcc00' : '1px solid var(--border-color)'),
                                    opacity: isFuture ? 0.3 : 1,
                                    transition: 'all 0.15s ease',
                                    flexShrink: 0,
                                    zIndex: 1,
                                    cursor: 'pointer'
                                }}
                            >
                                {getIconForAction(entry.label)}
                            </div>
                        </Tooltip>
                    );
                })}

                {/* Vertical Orange Marker (Equidistant between blocks) */}
                <div style={{
                    position: 'absolute',
                    left: `${scrubberPos - 1}px`,
                    top: '10%',
                    bottom: '10%',
                    width: '2px',
                    backgroundColor: '#ff6b00',
                    zIndex: 10,
                    pointerEvents: 'none',
                    transition: isDragging ? 'none' : 'left 0.2s cubic-bezier(0.18, 0.89, 0.32, 1.28)'
                }}>
                    <div style={{ position: 'absolute', top: 0, left: '-4px', right: '-4px', height: '2px', backgroundColor: '#ff6b00' }} />
                    <div style={{ position: 'absolute', bottom: 0, left: '-4px', right: '-4px', height: '2px', backgroundColor: '#ff6b00' }} />
                </div>
            </div>

            {/* Bulk Actions & Undo/Redo */}
            <div style={{
                display: 'flex',
                gap: '8px',
                paddingLeft: isMobile ? '8px' : '15px',
                borderLeft: '1px solid var(--border-color)',
                height: '70%',
                alignItems: 'center',
                flexShrink: 0
            }}>
                {selectedHistoryIds.length > 1 && (
                    <Tooltip content="Delete Selected">
                        <button onClick={() => deleteHistoryEntry(selectedHistoryIds)} style={{ ...iconOnlyBtnStyle, backgroundColor: '#ff4444', color: 'white', borderColor: '#ff4444' }}>
                            <Trash2 size={12} />
                        </button>
                    </Tooltip>
                )}
                <Tooltip content="Undo (Ctrl+Z)" icon={<Undo2 size={12} />}>
                    <button onClick={undo} disabled={historyIndex === 0} style={iconOnlyBtnStyle}>
                        <Undo2 size={14} />
                    </button>
                </Tooltip>
                <Tooltip content="Redo (Ctrl+Y)" icon={<Redo2 size={12} />}>
                    <button onClick={redo} disabled={historyIndex === history.length - 1} style={iconOnlyBtnStyle}>
                        <Redo2 size={14} />
                    </button>
                </Tooltip>
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
