import React, { useEffect } from 'react';
import { TopBar } from './components/Layout/TopBar';
import { Toolbar } from './components/Sidebar/Toolbar';
import { ScenePanel } from './components/Sidebar/ScenePanel';
import { SceneView } from './components/Viewport/SceneView';
import { HistoryTimeline } from './components/Timeline/HistoryTimeline';
import { useStore } from './store/useStore';
import { useIsMobile } from './utils/useIsMobile';
import { Analytics } from '@vercel/analytics/react';
import './styles/global.css';

const App: React.FC = () => {
  const isMobile = useIsMobile();
  const undo = useStore((state) => state.undo);
  const redo = useStore((state) => state.redo);
  const selectedHistoryIds = useStore((state) => state.selectedHistoryIds);
  const deleteHistoryEntry = useStore((state) => state.deleteHistoryEntry);
  const [isSidebarOpen, setSidebarOpen] = React.useState(!isMobile);

  useEffect(() => {
    // Sync sidebar state with mobile status
    setSidebarOpen(!isMobile);
  }, [isMobile]);

  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      // Don't trigger if user is typing in an input
      if (document.activeElement instanceof HTMLInputElement || document.activeElement instanceof HTMLTextAreaElement) {
        return;
      }

      if (e.ctrlKey || e.metaKey) {
        if (e.key === 'z') {
          if (e.shiftKey) {
            redo();
          } else {
            undo();
          }
        } else if (e.key === 'y') {
          redo();
        }
      } else if (e.key === 'Delete') {
        if (selectedHistoryIds.length > 0) {
          deleteHistoryEntry(selectedHistoryIds);
        }
      }
    };

    window.addEventListener('keydown', handleKeyDown);
    return () => window.removeEventListener('keydown', handleKeyDown);
  }, [undo, redo, selectedHistoryIds, deleteHistoryEntry]);

  return (
    <div style={{ height: '100vh', display: 'flex', flexDirection: 'column', overflow: 'hidden', backgroundColor: 'var(--bg-dark)' }}>
      <Analytics />
      <TopBar onToggleSidebar={() => setSidebarOpen(!isSidebarOpen)} />

      <div style={{
        flex: 1,
        display: 'flex',
        flexDirection: isMobile ? 'column' : 'row',
        overflow: 'hidden',
        position: 'relative'
      }}>

        {/* Toolbar: Side on Desktop, Top on Mobile */}
        <div className="toolbar-drawer">
          <Toolbar />
        </div>

        {/* This container ensures Viewport and Timeline are correctly stacked vertically */}
        <div style={{
          flex: 1,
          display: 'flex',
          flexDirection: 'column',
          overflow: 'hidden',
          position: 'relative'
        }}>

          <div style={{ flex: 1, display: 'flex', overflow: 'hidden', position: 'relative' }}>
            <SceneView />

            <div
              className={`sidebar-drawer ${!isSidebarOpen ? 'sidebar-closed' : ''}`}
              style={{ pointerEvents: 'auto' }}
            >
              <ScenePanel onClose={() => setSidebarOpen(false)} />
            </div>

            {/* Mobile Overlay for Sidebar */}
            {isMobile && isSidebarOpen && (
              <div
                onClick={() => setSidebarOpen(false)}
                style={{
                  position: 'absolute',
                  inset: 0,
                  backgroundColor: 'rgba(0,0,0,0.5)',
                  zIndex: 999,
                  backdropFilter: 'blur(2px)'
                }}
              />
            )}
          </div>

          {/* Timeline is always at the bottom of the content area */}
          <HistoryTimeline />
        </div>
      </div>
    </div>
  );
};

export default App;

