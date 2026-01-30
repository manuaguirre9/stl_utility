
import * as THREE from 'three';

/**
 * Subdivides a set of triangles.
 * Each triangle is split into 4 triangles by connecting midpoints of edges.
 * @param positions Array of vertex positions (flat [x,y,z, x,y,z, ...])
 * @param steps Number of subdivision iterations
 */
export function subdivideTriangles(positions: Float32Array, steps: number): Float32Array {
    if (steps <= 0) return positions;

    let currentPositions = positions;

    for (let s = 0; s < steps; s++) {
        const triangleCount = currentPositions.length / 9;
        const nextPositions = new Float32Array(triangleCount * 4 * 9);

        for (let i = 0; i < triangleCount; i++) {
            const offset = i * 9;

            const v1 = new THREE.Vector3(currentPositions[offset + 0], currentPositions[offset + 1], currentPositions[offset + 2]);
            const v2 = new THREE.Vector3(currentPositions[offset + 3], currentPositions[offset + 4], currentPositions[offset + 5]);
            const v3 = new THREE.Vector3(currentPositions[offset + 6], currentPositions[offset + 7], currentPositions[offset + 8]);

            const m1 = new THREE.Vector3().addVectors(v1, v2).multiplyScalar(0.5);
            const m2 = new THREE.Vector3().addVectors(v2, v3).multiplyScalar(0.5);
            const m3 = new THREE.Vector3().addVectors(v3, v1).multiplyScalar(0.5);

            const newTriangles = [
                v1, m1, m3,
                m1, v2, m2,
                m3, m2, v3,
                m1, m2, m3
            ];

            const nextOffset = i * 4 * 9;
            for (let t = 0; t < 4; t++) {
                const tOffset = nextOffset + (t * 9);
                nextPositions[tOffset + 0] = newTriangles[t * 3 + 0].x;
                nextPositions[tOffset + 1] = newTriangles[t * 3 + 0].y;
                nextPositions[tOffset + 2] = newTriangles[t * 3 + 0].z;

                nextPositions[tOffset + 3] = newTriangles[t * 3 + 1].x;
                nextPositions[tOffset + 4] = newTriangles[t * 3 + 1].y;
                nextPositions[tOffset + 5] = newTriangles[t * 3 + 1].z;

                nextPositions[tOffset + 6] = newTriangles[t * 3 + 2].x;
                nextPositions[tOffset + 7] = newTriangles[t * 3 + 2].y;
                nextPositions[tOffset + 8] = newTriangles[t * 3 + 2].z;
            }
        }
        currentPositions = nextPositions;
    }

    return currentPositions;
}

/**
 * Specifically subdivides selected faces until they reach a target triangle size.
 */
export function subdivideSelectedFacesToSize(
    geometry: THREE.BufferGeometry,
    faceIndices: number[],
    targetSize: number
): THREE.BufferGeometry {
    if (faceIndices.length === 0 || targetSize <= 0) return geometry;

    const posAttr = geometry.attributes.position;
    const totalFaces = posAttr.count / 3;
    const selectedIndicesSet = new Set(faceIndices);

    let currentPositions: Float32Array = new Float32Array(faceIndices.length * 9);
    faceIndices.forEach((faceIdx, i) => {
        const offset = faceIdx * 3;
        for (let v = 0; v < 3; v++) {
            currentPositions[i * 9 + v * 3 + 0] = posAttr.getX(offset + v);
            currentPositions[i * 9 + v * 3 + 1] = posAttr.getY(offset + v);
            currentPositions[i * 9 + v * 3 + 2] = posAttr.getZ(offset + v);
        }
    });

    const MAX_STEPS = 6;
    for (let s = 0; s < MAX_STEPS; s++) {
        const triCount = currentPositions.length / 9;
        let needsMore = false;

        for (let i = 0; i < triCount; i++) {
            const o = i * 9;
            const d12 = distSq(currentPositions[o], currentPositions[o + 1], currentPositions[o + 2], currentPositions[o + 3], currentPositions[o + 4], currentPositions[o + 5]);
            const d23 = distSq(currentPositions[o + 3], currentPositions[o + 4], currentPositions[o + 5], currentPositions[o + 6], currentPositions[o + 7], currentPositions[o + 8]);
            const d31 = distSq(currentPositions[o + 6], currentPositions[o + 7], currentPositions[o + 8], currentPositions[o], currentPositions[o + 1], currentPositions[o + 2]);

            if (Math.max(d12, d23, d31) > targetSize * targetSize) {
                needsMore = true;
                break;
            }
        }

        if (!needsMore) break;
        currentPositions = subdivideTriangles(currentPositions, 1);
    }

    const unselectedFaces: number[] = [];
    for (let i = 0; i < totalFaces; i++) {
        if (!selectedIndicesSet.has(i)) unselectedFaces.push(i);
    }

    const finalTriangleCount = unselectedFaces.length + (currentPositions.length / 9);
    const finalPositions = new Float32Array(finalTriangleCount * 9);

    unselectedFaces.forEach((faceIdx, i) => {
        const srcOffset = faceIdx * 3;
        for (let v = 0; v < 3; v++) {
            finalPositions[i * 9 + v * 3 + 0] = posAttr.getX(srcOffset + v);
            finalPositions[i * 9 + v * 3 + 1] = posAttr.getY(srcOffset + v);
            finalPositions[i * 9 + v * 3 + 2] = posAttr.getZ(srcOffset + v);
        }
    });

    finalPositions.set(currentPositions, unselectedFaces.length * 9);

    const newGeometry = new THREE.BufferGeometry();
    newGeometry.setAttribute('position', new THREE.BufferAttribute(finalPositions, 3));
    newGeometry.computeVertexNormals();
    newGeometry.computeBoundingBox();
    newGeometry.computeBoundingSphere();

    return newGeometry;
}

function distSq(x1: number, y1: number, z1: number, x2: number, y2: number, z2: number): number {
    return (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2;
}

/**
 * Specifically subdivides selected faces within a geometry.
 */
export function subdivideSelectedFaces(
    geometry: THREE.BufferGeometry,
    faceIndices: number[],
    steps: number
): THREE.BufferGeometry {
    if (steps <= 0 || faceIndices.length === 0) return geometry;

    const posAttr = geometry.attributes.position;
    const totalFaces = posAttr.count / 3;
    const selectedIndicesSet = new Set(faceIndices);

    const unselectedFaces: number[] = [];
    for (let i = 0; i < totalFaces; i++) {
        if (!selectedIndicesSet.has(i)) unselectedFaces.push(i);
    }

    const selectedPositions = new Float32Array(faceIndices.length * 9);
    faceIndices.forEach((faceIdx, i) => {
        const offset = faceIdx * 3;
        for (let v = 0; v < 3; v++) {
            selectedPositions[i * 9 + v * 3 + 0] = posAttr.getX(offset + v);
            selectedPositions[i * 9 + v * 3 + 1] = posAttr.getY(offset + v);
            selectedPositions[i * 9 + v * 3 + 2] = posAttr.getZ(offset + v);
        }
    });

    const subdividedPositions = subdivideTriangles(selectedPositions, steps);

    const finalTriangleCount = unselectedFaces.length + (subdividedPositions.length / 9);
    const finalPositions = new Float32Array(finalTriangleCount * 9);

    unselectedFaces.forEach((faceIdx, i) => {
        const srcOffset = faceIdx * 3;
        for (let v = 0; v < 3; v++) {
            finalPositions[i * 9 + v * 3 + 0] = posAttr.getX(srcOffset + v);
            finalPositions[i * 9 + v * 3 + 1] = posAttr.getY(srcOffset + v);
            finalPositions[i * 9 + v * 3 + 2] = posAttr.getZ(srcOffset + v);
        }
    });

    finalPositions.set(subdividedPositions, unselectedFaces.length * 9);

    const newGeometry = new THREE.BufferGeometry();
    newGeometry.setAttribute('position', new THREE.BufferAttribute(finalPositions, 3));
    newGeometry.computeVertexNormals();
    newGeometry.computeBoundingBox();
    newGeometry.computeBoundingSphere();

    return newGeometry;
}

/**
 * Estimates the circumference of the selected faces.
 * Useful for calculating seamless texture repeats.
 */
export function estimateSelectionCircumference(
    geometry: THREE.BufferGeometry,
    faceIndices: number[]
): number {
    if (faceIndices.length === 0) return 0;

    const posAttr = geometry.attributes.position;
    const box = new THREE.Box3();
    const centroid = new THREE.Vector3();
    let pointCount = 0;

    for (const faceIdx of faceIndices) {
        for (let v = 0; v < 3; v++) {
            const idx = faceIdx * 3 + v;
            const x = posAttr.getX(idx);
            const y = posAttr.getY(idx);
            const z = posAttr.getZ(idx);
            const p = new THREE.Vector3(x, y, z);
            box.expandByPoint(p);
            centroid.add(p);
            pointCount++;
        }
    }
    centroid.divideScalar(pointCount);

    const size = new THREE.Vector3();
    box.getSize(size);

    let radiusSum = 0;

    if (size.z > size.x && size.z > size.y) {
        for (const faceIdx of faceIndices) {
            for (let v = 0; v < 3; v++) {
                const idx = faceIdx * 3 + v;
                const dx = posAttr.getX(idx) - centroid.x;
                const dy = posAttr.getY(idx) - centroid.y;
                radiusSum += Math.sqrt(dx * dx + dy * dy);
            }
        }
    } else if (size.y > size.x && size.y > size.z) {
        for (const faceIdx of faceIndices) {
            for (let v = 0; v < 3; v++) {
                const idx = faceIdx * 3 + v;
                const dx = posAttr.getX(idx) - centroid.x;
                const dz = posAttr.getZ(idx) - centroid.z;
                radiusSum += Math.sqrt(dx * dx + dz * dz);
            }
        }
    } else {
        for (const faceIdx of faceIndices) {
            for (let v = 0; v < 3; v++) {
                const idx = faceIdx * 3 + v;
                const dy = posAttr.getY(idx) - centroid.y;
                const dz = posAttr.getZ(idx) - centroid.z;
                radiusSum += Math.sqrt(dy * dy + dz * dz);
            }
        }
    }

    const avgRadius = radiusSum / pointCount;
    return 2 * Math.PI * avgRadius;
}

/**
 * Groups face indices into connected islands/patches.
 * Two faces are in the same island if they share at least one vertex
 * AND (optional) the angle between their normals is below a threshold.
 */
export function getContiguousIslands(
    geometry: THREE.BufferGeometry,
    faceIndices: number[],
    angleThresholdDeg: number = 20 // Lower threshold to split cylinder/cone transitions
): number[][] {
    if (faceIndices.length === 0) return [];

    const posAttr = geometry.attributes.position;
    const vertexToFaces = new Map<string, number[]>();

    // Pre-calculate normals for all selected faces for efficient angle comparison
    const faceNormals = new Map<number, THREE.Vector3>();
    for (const fIdx of faceIndices) {
        const off = fIdx * 3;
        const vA = new THREE.Vector3().fromBufferAttribute(posAttr, off);
        const vB = new THREE.Vector3().fromBufferAttribute(posAttr, off + 1);
        const vC = new THREE.Vector3().fromBufferAttribute(posAttr, off + 2);
        const normal = new THREE.Vector3()
            .crossVectors(new THREE.Vector3().subVectors(vB, vA), new THREE.Vector3().subVectors(vC, vA))
            .normalize();
        faceNormals.set(fIdx, normal);
    }

    const dotThreshold = Math.cos((angleThresholdDeg * Math.PI) / 180);

    // Key function for vertices to handle slight floating point issues
    const getVKey = (fIdx: number, vIdxInFace: number) => {
        const off = fIdx * 3 + vIdxInFace;
        return `${Math.round(posAttr.getX(off) * 10000)},${Math.round(posAttr.getY(off) * 10000)},${Math.round(posAttr.getZ(off) * 10000)}`;
    };

    // Build vertex-to-face adjacency (only for the selected faces)
    for (const fIdx of faceIndices) {
        for (let v = 0; v < 3; v++) {
            const key = getVKey(fIdx, v);
            if (!vertexToFaces.has(key)) vertexToFaces.set(key, []);
            vertexToFaces.get(key)!.push(fIdx);
        }
    }

    const visited = new Set<number>();
    const islands: number[][] = [];

    for (const fIdx of faceIndices) {
        if (visited.has(fIdx)) continue;

        const island: number[] = [];
        const stack = [fIdx];
        visited.add(fIdx);

        while (stack.length > 0) {
            const current = stack.pop()!;
            island.push(current);
            const currentNormal = faceNormals.get(current)!;

            // Check all 3 vertices of current face for adjacent selected faces
            for (let v = 0; v < 3; v++) {
                const key = getVKey(current, v);
                const candidates = vertexToFaces.get(key) || [];
                for (const adj of candidates) {
                    if (!visited.has(adj)) {
                        const adjNormal = faceNormals.get(adj)!;
                        // Only join if normals are similar enough
                        if (currentNormal.dot(adjNormal) >= dotThreshold) {
                            visited.add(adj);
                            stack.push(adj);
                        }
                    }
                }
            }
        }
        islands.push(island);
    }

    return islands;
}
