import { STLExporter } from 'three-stdlib';
import type { ModelData } from '../store/useStore';
import * as THREE from 'three';

export const exportSceneToSTL = (models: ModelData[]) => {
    if (models.length === 0) {
        console.warn('No models to export');
        return;
    }

    const exporter = new STLExporter();

    // Create a temporary scene or group to hold visible meshes with their transforms applied
    // Or we can just process geometry directly if we want to merge them, but STLExporter handles object arrays.

    const exportObjects: THREE.Object3D[] = [];

    models.forEach(model => {
        if (!model.visible) return;

        // Reconstruct the mesh with transforms for export
        const mesh = new THREE.Mesh(model.bufferGeometry);
        mesh.position.set(...model.position);
        mesh.rotation.set(...model.rotation);
        mesh.scale.set(...model.scale);
        mesh.updateMatrixWorld();

        exportObjects.push(mesh);
    });

    if (exportObjects.length === 0) {
        alert("No visible objects to export.");
        return;
    }

    // Export scene (or individual objects if we want separate files, but usually one scene export)
    // If we pass an array, STLExporter might expect a single Object3D usually. 
    // It's safer to add them to a group.
    const group = new THREE.Group();
    exportObjects.forEach(obj => group.add(obj));

    const result = exporter.parse(group, { binary: true });
    const blob = new Blob([result as any], { type: 'application/octet-stream' });
    const url = URL.createObjectURL(blob);

    const link = document.createElement('a');
    link.href = url;
    link.download = 'scene_export.stl';
    link.click();

    URL.revokeObjectURL(url);
};
