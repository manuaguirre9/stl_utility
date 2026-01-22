import * as THREE from 'three';
import { STLLoader } from 'three/addons/loaders/STLLoader.js';
import { STLExporter } from 'three/addons/exporters/STLExporter.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

/**
 * STL Fuzzy Sculpt Pro - Core Application
 */

class App {
    constructor() {
        this.initScene();
        this.initListeners();
        this.history = [];
        this.historyIndex = -1;
        this.maxHistory = 20;

        this.brush = {
            size: 20,
            thickness: 0.3,
            density: 0.8,
            isPainting: false
        };

        this.modelState = {
            mesh: null,
            originalVertices: null,
            maxDim: 1,
            adjacency: null
        };

        // UI Refs
        this.overlay = document.getElementById('overlay');
        this.loadingText = document.getElementById('loading-text');
        this.progressContainer = document.getElementById('progress-container');
        this.progressBar = document.getElementById('progress-bar');

        this.hideOverlay();
        this.animate();
    }

    initScene() {
        this.canvas = document.getElementById('viewport');
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x0d0d0e);

        this.camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 5000);
        this.camera.position.set(100, 100, 100);

        this.renderer = new THREE.WebGLRenderer({
            canvas: this.canvas,
            antialias: true,
            alpha: true
        });
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        this.renderer.toneMapping = THREE.ACESFilmicToneMapping;

        // Lighting (Backup for non-matcap materials)
        const ambient = new THREE.AmbientLight(0xffffff, 0.4);
        this.scene.add(ambient);

        // Grid & Helpers
        const grid = new THREE.GridHelper(200, 40, 0x222222, 0x111111);
        grid.position.y = -0.1;
        this.scene.add(grid);

        // Controls
        this.controls = new OrbitControls(this.camera, this.canvas);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;
        this.controls.enableZoom = false; // Custom zoom
        this.controls.mouseButtons = {
            LEFT: null, // Painting
            MIDDLE: THREE.MOUSE.DOLLY,
            RIGHT: THREE.MOUSE.ROTATE
        };

        // Raycaster for interaction
        this.raycaster = new THREE.Raycaster();
        this.mouse = new THREE.Vector2();

        // 3D Brush Cursor
        const ringGeo = new THREE.RingGeometry(0.9, 1, 32);
        const ringMat = new THREE.MeshBasicMaterial({ color: 0x3b82f6, side: THREE.DoubleSide });
        this.brushCursor = new THREE.Mesh(ringGeo, ringMat);
        this.brushCursor.visible = false;
        this.scene.add(this.brushCursor);

        // Matcap Loader
        const loader = new THREE.TextureLoader();
        // Using a reliable high-quality matcap from a common CDN
        this.matcapTexture = loader.load('https://raw.githubusercontent.com/nidorx/matcaps/master/thumbnail/renders/666666_cccccc_999999_bbbbbb.jpg');
    }

    initListeners() {
        window.addEventListener('resize', () => this.onResize());

        // Custom Zoom
        this.canvas.addEventListener('wheel', (e) => {
            e.preventDefault();
            const factor = 1.1;
            if (e.deltaY > 0) this.controls.dollyOut(factor);
            else this.controls.dollyIn(factor);
            this.controls.update();
        }, { passive: false });

        // Drag & Drop
        const dropZone = document.getElementById('drop-zone');
        const fileInput = document.getElementById('file-input');

        fileInput.addEventListener('change', (e) => this.loadFile(e.target.files[0]));

        // Brush Controls
        this.setupSlider('input-brush-size', 'val-brush-size', (v) => this.brush.size = parseFloat(v));
        this.setupSlider('input-thickness', 'val-thickness', (v) => this.brush.thickness = parseFloat(v));
        this.setupSlider('input-density', 'val-density', (v) => this.brush.density = parseFloat(v));

        // Painting Events
        this.canvas.addEventListener('mousedown', (e) => {
            if (e.button === 0 && this.modelState.mesh) {
                this.brush.isPainting = true;
                this.saveHistory(); // Start of stroke
                this.paint(e);
            }
        });

        window.addEventListener('mousemove', (e) => {
            this.updateBrushCursor(e);
            if (this.brush.isPainting) this.paint(e);
        });

        window.addEventListener('mouseup', () => this.brush.isPainting = false);

        // Actions
        document.getElementById('btn-undo').onclick = () => this.undo();
        document.getElementById('btn-redo').onclick = () => this.redo();
        document.getElementById('btn-reset').onclick = () => this.resetGeometry();
        document.getElementById('btn-export').onclick = () => this.exportSTL();

        // View Snapping
        document.getElementById('view-top').onclick = () => this.snapView('top');
        document.getElementById('view-front').onclick = () => this.snapView('front');
        document.getElementById('view-side').onclick = () => this.snapView('side');
        document.getElementById('toggle-wireframe').onclick = () => this.toggleWireframe();

        // Keyboard
        window.addEventListener('keydown', (e) => {
            if (e.ctrlKey && e.key === 'z') this.undo();
            if (e.ctrlKey && e.key === 'y') this.redo();
        });

        this.canvas.addEventListener('dblclick', (e) => this.onDoubleClick(e));
    }

    setupSlider(id, valId, callback) {
        const input = document.getElementById(id);
        const display = document.getElementById(valId);
        input.addEventListener('input', (e) => {
            display.textContent = parseFloat(e.target.value).toFixed(2);
            callback(e.target.value);
        });
    }

    async loadFile(file) {
        if (!file || !file.name.toLowerCase().endsWith('.stl')) return;

        this.showOverlay(`LOADING ${file.name.toUpperCase()}`);

        const reader = new FileReader();
        reader.onload = async (e) => {
            try {
                const loader = new STLLoader();
                const geometry = loader.parse(e.target.result);
                this.setupModel(geometry);
            } catch (err) {
                alert("Error loading STL: " + err.message);
                this.hideOverlay();
            }
        };
        reader.readAsArrayBuffer(file);
    }

    setupModel(geometry) {
        if (this.modelState.mesh) {
            this.scene.remove(this.modelState.mesh);
            this.modelState.mesh.geometry.dispose();
            this.modelState.mesh.material.dispose();
        }

        geometry.center();
        geometry.computeBoundingBox();
        const size = new THREE.Vector3();
        geometry.boundingBox.getSize(size);
        this.modelState.maxDim = Math.max(size.x, size.y, size.z);

        const material = new THREE.MeshMatcapMaterial({
            matcap: this.matcapTexture,
            color: 0xffffff,
            side: THREE.DoubleSide
        });

        this.modelState.mesh = new THREE.Mesh(geometry, material);
        this.scene.add(this.modelState.mesh);

        // Initial Snapping
        this.camera.position.set(this.modelState.maxDim * 1.5, this.modelState.maxDim * 1.5, this.modelState.maxDim * 1.5);
        this.controls.target.set(0, 0, 0);
        this.controls.update();

        // Stats
        document.getElementById('stat-verts').textContent = `Vertices: ${geometry.attributes.position.count}`;
        document.getElementById('stat-dims').textContent = `Dim: ${size.x.toFixed(1)} x ${size.y.toFixed(1)} x ${size.z.toFixed(1)} mm`;

        // Save original for reset
        this.modelState.originalVertices = geometry.attributes.position.array.slice();
        this.history = [];
        this.historyIndex = -1;
        this.saveHistory(); // Initial state

        document.getElementById('controls-container').style.display = 'block';
        this.hideOverlay();

        // Async Adjacency
        this.buildAdjacency();
    }

    async buildAdjacency() {
        if (!this.modelState.mesh) return;
        this.showOverlay("ANALYZING MESH...");

        return new Promise((resolve) => {
            setTimeout(() => {
                const geometry = this.modelState.mesh.geometry;
                const pos = geometry.attributes.position;
                const triangleCount = pos.count / 3;

                // Map vertex position to triangle indices
                const vertexToTriangles = new Map();
                const precision = 5;

                for (let i = 0; i < triangleCount; i++) {
                    if (i % 1000 === 0) this.setProgress((i / triangleCount) * 50); // First 50%
                    for (let j = 0; j < 3; j++) {
                        const idx = i * 3 + j;
                        const x = pos.getX(idx).toFixed(precision);
                        const y = pos.getY(idx).toFixed(precision);
                        const z = pos.getZ(idx).toFixed(precision);
                        const key = `${x},${y},${z}`;

                        if (!vertexToTriangles.has(key)) vertexToTriangles.set(key, []);
                        vertexToTriangles.get(key).push(i);
                    }
                }

                this.modelState.adjacency = new Array(triangleCount).fill(0).map(() => new Set());

                for (let i = 0; i < triangleCount; i++) {
                    if (i % 1000 === 0) this.setProgress(50 + (i / triangleCount) * 50); // Last 50%
                    for (let j = 0; j < 3; j++) {
                        const idx = i * 3 + j;
                        const key = `${pos.getX(idx).toFixed(precision)},${pos.getY(idx).toFixed(precision)},${pos.getZ(idx).toFixed(precision)}`;
                        const sharingTriangles = vertexToTriangles.get(key);
                        for (const otherIdx of sharingTriangles) {
                            if (otherIdx !== i) this.modelState.adjacency[i].add(otherIdx);
                        }
                    }
                }

                this.hideOverlay();
                console.log("Adjacency Graph Built", triangleCount, "triangles");
                resolve();
            }, 10);
        });
    }

    onDoubleClick(event) {
        if (!this.modelState.mesh || !this.modelState.adjacency) return;

        this.mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
        this.mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

        this.raycaster.setFromCamera(this.mouse, this.camera);
        const intersects = this.raycaster.intersectObject(this.modelState.mesh);

        if (intersects.length > 0) {
            const startFace = intersects[0].faceIndex;
            this.smartSelectRegion(startFace);
        }
    }

    smartSelectRegion(startFace) {
        const geo = this.modelState.mesh.geometry;
        const normals = geo.attributes.normal;
        const adjacency = this.modelState.adjacency;

        const startNormal = new THREE.Vector3().fromBufferAttribute(normals, startFace * 3);
        const queue = [startFace];
        const visited = new Set([startFace]);
        const selection = new Set([startFace]);

        const threshold = Math.cos(30 * (Math.PI / 180)); // 30 degrees

        const currentNormal = new THREE.Vector3();
        const neighborNormal = new THREE.Vector3();

        while (queue.length > 0) {
            const current = queue.shift();
            currentNormal.fromBufferAttribute(normals, current * 3);

            for (const neighbor of adjacency[current]) {
                if (visited.has(neighbor)) continue;

                neighborNormal.fromBufferAttribute(normals, neighbor * 3);

                // Normal angle check
                if (currentNormal.dot(neighborNormal) > threshold) {
                    visited.add(neighbor);
                    selection.add(neighbor);
                    queue.push(neighbor);
                }
            }
        }

        this.applyFuzzyToSelection(selection);
    }

    applyFuzzyToSelection(selection) {
        this.saveHistory();
        const geo = this.modelState.mesh.geometry;
        const pos = geo.attributes.position;
        const normals = geo.attributes.normal;
        const strength = this.brush.thickness;

        // Find unique vertices in selection
        const vertexIndices = new Set();
        selection.forEach(faceIdx => {
            vertexIndices.add(faceIdx * 3);
            vertexIndices.add(faceIdx * 3 + 1);
            vertexIndices.add(faceIdx * 3 + 2);
        });

        const tempPoint = new THREE.Vector3();
        const vNormal = new THREE.Vector3();
        const strokeMap = new Map();

        vertexIndices.forEach(i => {
            tempPoint.fromBufferAttribute(pos, i);
            const key = `${tempPoint.x.toFixed(4)},${tempPoint.y.toFixed(4)},${tempPoint.z.toFixed(4)}`;

            let noise;
            if (strokeMap.has(key)) {
                noise = strokeMap.get(key);
            } else {
                noise = (Math.random() * 2 - 1) * strength;
                strokeMap.set(key, noise);
            }

            vNormal.fromBufferAttribute(normals, i).applyQuaternion(this.modelState.mesh.quaternion);
            tempPoint.addScaledVector(vNormal, noise);
            this.modelState.mesh.worldToLocal(tempPoint);
            pos.setXYZ(i, tempPoint.x, tempPoint.y, tempPoint.z);
        });

        pos.needsUpdate = true;
        geo.computeVertexNormals();
        console.log(`Applied fuzzy to ${selection.size} faces.`);
    }
    updateBrushCursor(event) {
        if (!this.modelState.mesh) return;

        this.mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
        this.mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

        this.raycaster.setFromCamera(this.mouse, this.camera);
        const intersects = this.raycaster.intersectObject(this.modelState.mesh);

        if (intersects.length > 0) {
            const hit = intersects[0];
            this.brushCursor.visible = true;
            this.brushCursor.position.copy(hit.point).addScaledVector(hit.face.normal, 0.1);

            // Align with normal
            const lookAtTarget = hit.point.clone().add(hit.face.normal);
            this.brushCursor.lookAt(lookAtTarget);

            // Scale cursor based on UI size
            const brushSizeWorld = (this.brush.size / 100) * this.modelState.maxDim * 0.5;
            this.brushCursor.scale.set(brushSizeWorld, brushSizeWorld, 1);
        } else {
            this.brushCursor.visible = false;
        }
    }

    paint(event) {
        if (!this.modelState.mesh) return;

        this.raycaster.setFromCamera(this.mouse, this.camera);
        const intersects = this.raycaster.intersectObject(this.modelState.mesh);

        if (intersects.length > 0) {
            const hit = intersects[0];
            this.applyFuzzy(hit.point, hit.face.normal);
        }
    }

    applyFuzzy(center, normal) {
        const geo = this.modelState.mesh.geometry;
        const pos = geo.attributes.position;
        const radius = (this.brush.size / 100) * this.modelState.maxDim * 0.5;
        const strength = this.brush.thickness;

        const tempPoint = new THREE.Vector3();
        const vNormal = new THREE.Vector3();
        const normals = geo.attributes.normal;

        // Stroke map for watertightness
        const strokeMap = new Map();

        for (let i = 0; i < pos.count; i++) {
            tempPoint.fromBufferAttribute(pos, i);
            this.modelState.mesh.localToWorld(tempPoint);

            const dist = tempPoint.distanceTo(center);
            if (dist < radius) {
                const key = `${tempPoint.x.toFixed(4)},${tempPoint.y.toFixed(4)},${tempPoint.z.toFixed(4)}`;

                let noise;
                if (strokeMap.has(key)) {
                    noise = strokeMap.get(key);
                } else {
                    const falloff = 1 - (dist / radius);
                    noise = (Math.random() * 2 - 1) * strength * falloff;
                    strokeMap.set(key, noise);
                }

                if (normals) vNormal.fromBufferAttribute(normals, i).applyQuaternion(this.modelState.mesh.quaternion);
                else vNormal.copy(normal);

                tempPoint.addScaledVector(vNormal, noise);
                this.modelState.mesh.worldToLocal(tempPoint);
                pos.setXYZ(i, tempPoint.x, tempPoint.y, tempPoint.z);
            }
        }

        pos.needsUpdate = true;
        geo.computeVertexNormals();
    }

    /**
     * History System
     */
    saveHistory() {
        if (!this.modelState.mesh) return;

        // Remove forward history if we are in the middle of a stack
        if (this.historyIndex < this.history.length - 1) {
            this.history = this.history.slice(0, this.historyIndex + 1);
        }

        const snapshot = this.modelState.mesh.geometry.attributes.position.array.slice();
        this.history.push(snapshot);

        if (this.history.length > this.maxHistory) {
            this.history.shift();
        } else {
            this.historyIndex++;
        }
        this.updateHistoryButtons();
    }

    undo() {
        if (this.historyIndex <= 0) return;
        this.historyIndex--;
        this.restoreHistory(this.history[this.historyIndex]);
        this.updateHistoryButtons();
    }

    redo() {
        if (this.historyIndex >= this.history.length - 1) return;
        this.historyIndex++;
        this.restoreHistory(this.history[this.historyIndex]);
        this.updateHistoryButtons();
    }

    restoreHistory(data) {
        const pos = this.modelState.mesh.geometry.attributes.position;
        pos.array.set(data);
        pos.needsUpdate = true;
        this.modelState.mesh.geometry.computeVertexNormals();
    }

    updateHistoryButtons() {
        document.getElementById('btn-undo').disabled = this.historyIndex <= 0;
        document.getElementById('btn-redo').disabled = this.historyIndex >= this.history.length - 1;
    }

    resetGeometry() {
        if (!this.modelState.originalVertices) return;
        this.restoreHistory(this.modelState.originalVertices);
        this.saveHistory();
    }

    /**
     * Viewport Actions
     */
    snapView(side) {
        const d = this.modelState.maxDim * 2;
        switch (side) {
            case 'top': this.camera.position.set(0, d, 0); break;
            case 'front': this.camera.position.set(0, 0, d); break;
            case 'side': this.camera.position.set(d, 0, 0); break;
        }
        this.controls.target.set(0, 0, 0);
        this.controls.update();
    }

    toggleWireframe() {
        if (!this.modelState.mesh) return;
        const mat = this.modelState.mesh.material;
        mat.wireframe = !mat.wireframe;
        document.getElementById('toggle-wireframe').classList.toggle('active', mat.wireframe);
    }

    exportSTL() {
        if (!this.modelState.mesh) return;
        const exporter = new STLExporter();
        const result = exporter.parse(this.modelState.mesh);
        const blob = new Blob([result], { type: 'application/octet-stream' });
        const link = document.createElement('a');
        link.href = URL.createObjectURL(blob);
        link.download = 'refined_fuzzy_model.stl';
        link.click();
    }

    onResize() {
        this.camera.aspect = window.innerWidth / window.innerHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(window.innerWidth, window.innerHeight);
    }

    showOverlay(text, showProgress = true) {
        this.loadingText.textContent = text;
        this.overlay.style.display = 'flex';
        this.overlay.style.opacity = '1';
        if (showProgress) {
            this.progressContainer.style.display = 'block';
            this.setProgress(0);
        } else {
            this.progressContainer.style.display = 'none';
        }
    }

    setProgress(percent) {
        if (this.progressBar) {
            this.progressBar.style.width = `${percent}%`;
        }
    }

    hideOverlay() {
        this.setProgress(100);
        this.overlay.style.opacity = '0';
        setTimeout(() => {
            this.overlay.style.display = 'none';
            this.setProgress(0);
        }, 500);
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }
}

// Start App
new App();
