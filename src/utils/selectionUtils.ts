
import * as THREE from 'three';

// Optimization: O(N) edge-based adjacency build.
// Instead of O(N^2) face pairs per vertex, we hash each of the 3 edges of a face.
// If two faces share an edge hash, they are neighbors.
export const buildAdjacencyGraph = (geometry: THREE.BufferGeometry, precision = 1e-3) => {
    const position = geometry.attributes.position;
    const index = geometry.index;
    const faceCount = index ? index.count / 3 : position.count / 3;

    const vKey = (idx: number) => {
        const x = Math.round(position.getX(idx) / precision);
        const y = Math.round(position.getY(idx) / precision);
        const z = Math.round(position.getZ(idx) / precision);
        return `${x},${y},${z}`;
    };

    // Map: edgeKey -> [faceIdx1, faceIdx2, ...]
    const edgeToFaces = new Map<string, number[]>();

    for (let f = 0; f < faceCount; f++) {
        const i0 = index ? index.getX(f * 3 + 0) : f * 3 + 0;
        const i1 = index ? index.getX(f * 3 + 1) : f * 3 + 1;
        const i2 = index ? index.getX(f * 3 + 2) : f * 3 + 2;

        const k0 = vKey(i0), k1 = vKey(i1), k2 = vKey(i2);

        // Face has 3 edges: (k0,k1), (k1,k2), (k2,k0)
        const edges = [
            k0 < k1 ? `${k0}|${k1}` : `${k1}|${k0}`,
            k1 < k2 ? `${k1}|${k2}` : `${k2}|${k1}`,
            k2 < k0 ? `${k2}|${k0}` : `${k0}|${k2}`
        ];

        for (const edge of edges) {
            let list = edgeToFaces.get(edge);
            if (!list) {
                list = [];
                edgeToFaces.set(edge, list);
            }
            list.push(f);
        }
    }

    const adjacency: number[][] = Array.from({ length: faceCount }, () => []);

    for (const neighbors of edgeToFaces.values()) {
        if (neighbors.length < 2) continue;
        // Typically length is 2 for a manifold edge.
        // If it's > 2, it's non-manifold, but we still link them for selection purposes.
        for (let i = 0; i < neighbors.length; i++) {
            for (let j = i + 1; j < neighbors.length; j++) {
                const a = neighbors[i];
                const b = neighbors[j];
                adjacency[a].push(b);
                adjacency[b].push(a);
            }
        }
    }

    return adjacency;
};

/**
 * Smart selection: uses region growing and dynamic model fitting 
 * to distinguish between planar and cylindrical features.
 */
export const getConnectedFaces = (
    geometry: THREE.BufferGeometry,
    startFace: number,
    thresholdAngleDeg: number,
    adjacency?: number[][]
) => {
    if (!geometry.attributes.normal) geometry.computeVertexNormals();

    const pos = geometry.attributes.position;
    const index = geometry.index;
    const faceCount = index ? index.count / 3 : pos.count / 3;
    const normals: THREE.Vector3[] = [];

    // Precalculate normals for each face (triangle)
    const pA = new THREE.Vector3(), pB = new THREE.Vector3(), pC = new THREE.Vector3();
    const cb = new THREE.Vector3(), ab = new THREE.Vector3();

    for (let f = 0; f < faceCount; f++) {
        const vA_idx = index ? index.getX(f * 3 + 0) : f * 3 + 0;
        const vB_idx = index ? index.getX(f * 3 + 1) : f * 3 + 1;
        const vC_idx = index ? index.getX(f * 3 + 2) : f * 3 + 2;

        pA.fromBufferAttribute(pos, vA_idx);
        pB.fromBufferAttribute(pos, vB_idx);
        pC.fromBufferAttribute(pos, vC_idx);
        cb.subVectors(pC, pB);
        ab.subVectors(pA, pB);
        cb.cross(ab).normalize();
        normals.push(cb.clone());
    }

    const graph = adjacency || buildAdjacencyGraph(geometry);
    const startNormal = normals[startFace];

    /**
     * STEP 1: INITIAL PATCH GROWING (FOR CLASSIFICATION)
     * We grow a small representative patch to determine the feature type reliably.
     */
    const patchLimit = 25;
    const patchIndices: number[] = [startFace];
    const patchVisited = new Set<number>([startFace]);
    const patchQueue = [startFace];

    while (patchQueue.length > 0 && patchIndices.length < patchLimit) {
        const curr = patchQueue.shift()!;
        for (const neighbor of graph[curr]) {
            if (!patchVisited.has(neighbor)) {
                if (startNormal.dot(normals[neighbor]) > 0.7) {
                    patchVisited.add(neighbor);
                    patchIndices.push(neighbor);
                    patchQueue.push(neighbor);
                }
            }
            if (patchIndices.length >= patchLimit) break;
        }
    }

    /**
     * STEP 2: FEATURE ANALYSIS (PLANE VS CYLINDER)
     */

    // Plane evaluation: Average Normal & Deviation
    const avgNorm = new THREE.Vector3(0, 0, 0);
    patchIndices.forEach(idx => avgNorm.add(normals[idx]));
    avgNorm.normalize();

    let planeResidual = 0;
    patchIndices.forEach(idx => {
        planeResidual += (1.0 - Math.abs(avgNorm.dot(normals[idx])));
    });
    planeResidual /= patchIndices.length;

    // Cylinder evaluation: Axis Estimation
    let bestAxis = new THREE.Vector3(0, 1, 0);
    let bestAxisPoint = new THREE.Vector3(0, 0, 0);
    let avgRadius = 0;
    let cylinderResidual = 1.0;

    if (patchIndices.length > 5) {
        let n1 = startFace, n2 = -1;
        let minDot = 1.0;
        for (let i = 0; i < patchIndices.length; i++) {
            for (let j = i + 1; j < patchIndices.length; j++) {
                const d = Math.abs(normals[patchIndices[i]].dot(normals[patchIndices[j]]));
                if (d < minDot) {
                    minDot = d;
                    n1 = patchIndices[i];
                    n2 = patchIndices[j];
                }
            }
        }

        if (minDot < 0.9999) {
            const d1 = new THREE.Vector3().subVectors(normals[n1], startNormal);
            const d2 = new THREE.Vector3().subVectors(normals[n2], startNormal);
            bestAxis.crossVectors(d1, d2).normalize();

            if (bestAxis.length() > 0.5) {
                const p0 = new THREE.Vector3().fromBufferAttribute(pos, startFace * 3).projectOnPlane(bestAxis);
                const v0 = startNormal.clone().projectOnPlane(bestAxis).normalize();
                const p1 = new THREE.Vector3().fromBufferAttribute(pos, n1 * 3).projectOnPlane(bestAxis);
                const v1 = normals[n1].clone().projectOnPlane(bestAxis).normalize();

                const w0 = new THREE.Vector3().subVectors(p0, p1);
                const a = v0.dot(v0), b = v0.dot(v1), c = v1.dot(v1), d = v0.dot(w0), e = v1.dot(w0);
                const denom = a * c - b * b;

                if (Math.abs(denom) > 1e-6) {
                    const t = (b * e - c * d) / denom;
                    bestAxisPoint.copy(p0).addScaledVector(v0, t);

                    let radiusSum = 0;
                    const axisLine = new THREE.Line3(bestAxisPoint, bestAxisPoint.clone().add(bestAxis));
                    const tempP = new THREE.Vector3(), closest = new THREE.Vector3();

                    patchIndices.forEach(idx => {
                        tempP.fromBufferAttribute(pos, idx * 3);
                        axisLine.closestPointToPoint(tempP, false, closest);
                        radiusSum += tempP.distanceTo(closest);
                    });
                    avgRadius = radiusSum / patchIndices.length;

                    let radDev = 0;
                    patchIndices.forEach(idx => {
                        tempP.fromBufferAttribute(pos, idx * 3);
                        axisLine.closestPointToPoint(tempP, false, closest);
                        radDev += Math.abs(tempP.distanceTo(closest) - avgRadius);
                    });
                    cylinderResidual = radDev / (avgRadius + 1e-6) / patchIndices.length;
                }
            }
        }
    }

    /**
     * STEP 3: FINAL REGION GROWING
     */
    // Decision logic to select the best primitive
    const isPlane = planeResidual < cylinderResidual || cylinderResidual > 0.1;
    const isCylinder = !isPlane && cylinderResidual < 0.05 && avgRadius > 0.1;

    const visited = new Set<number>();
    const queue = [startFace];
    visited.add(startFace);
    const selected: number[] = [];

    const coplanarTolerance = 0.9999;
    const thresholdDot = Math.cos(thresholdAngleDeg * Math.PI / 180);
    const radiusTol = Math.max(avgRadius * 0.1, 0.1);

    while (queue.length > 0) {
        const current = queue.shift()!;
        selected.push(current);

        const neighbors = graph[current];
        for (const neighbor of neighbors) {
            if (visited.has(neighbor)) continue;

            const neighborNormal = normals[neighbor];
            let shouldAdd = false;

            if (isPlane) {
                // Plane Expansion
                if (avgNorm.dot(neighborNormal) >= coplanarTolerance) {
                    shouldAdd = true;
                }
            } else if (isCylinder) {
                // Cylinder Expansion (Radius check)
                const vA_idx = index ? index.getX(neighbor * 3 + 0) : neighbor * 3 + 0;
                const pN = new THREE.Vector3().fromBufferAttribute(pos, vA_idx);
                const axisL = new THREE.Line3(bestAxisPoint, bestAxisPoint.clone().add(bestAxis));
                const closeP = new THREE.Vector3();
                axisL.closestPointToPoint(pN, false, closeP);
                const dist = pN.distanceTo(closeP);

                if (Math.abs(dist - avgRadius) < radiusTol) {
                    // Safety check to ensure we stay on the same smooth feature
                    if (normals[current].dot(neighborNormal) > 0.8) {
                        shouldAdd = true;
                    }
                }
            } else {
                // Fallback (Generic curvature)
                if (normals[current].dot(neighborNormal) >= thresholdDot) {
                    shouldAdd = true;
                }
            }

            if (shouldAdd) {
                visited.add(neighbor);
                queue.push(neighbor);
            }
        }
    }

    return selected;
};

export type FeatureType = 'plane' | 'cylinder' | 'cone' | 'generic';

/**
 * Computes the characteristic of the normal distribution using PCA.
 * Returns eigenvalues and eigenvectors of the covariance matrix.
 */
const computeNormalPCA = (normals: THREE.Vector3[]) => {
    const n = normals.length;
    let meanX = 0, meanY = 0, meanZ = 0;
    normals.forEach(v => { meanX += v.x; meanY += v.y; meanZ += v.z; });
    meanX /= n; meanY /= n; meanZ /= n;

    let cxx = 0, cxy = 0, cxz = 0, cyy = 0, cyz = 0, czz = 0;
    normals.forEach(v => {
        const dx = v.x - meanX, dy = v.y - meanY, dz = v.z - meanZ;
        cxx += dx * dx; cxy += dx * dy; cxz += dx * dz;
        cyy += dy * dy; cyz += dy * dz; czz += dz * dz;
    });

    const cov = [
        [cxx / n, cxy / n, cxz / n],
        [cxy / n, cyy / n, cyz / n],
        [cxz / n, cyz / n, czz / n]
    ];

    const eigVals: number[] = [0, 0, 0];
    const eigVecs: THREE.Vector3[] = [new THREE.Vector3(), new THREE.Vector3(), new THREE.Vector3()];

    // Jacobi Rotation for 3x3 Symmetric Matrix
    const A = [cov[0].slice(), cov[1].slice(), cov[2].slice()];
    const V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

    for (let iter = 0; iter < 10; iter++) {
        let p = 0, q = 1;
        if (Math.abs(A[0][2]) > Math.abs(A[p][q])) { p = 0; q = 2; }
        if (Math.abs(A[1][2]) > Math.abs(A[p][q])) { p = 1; q = 2; }

        if (Math.abs(A[p][q]) < 1e-9) break;

        const theta = 0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]);
        const c = Math.cos(theta), s = Math.sin(theta);

        const app = c * c * A[p][p] - 2 * s * c * A[p][q] + s * s * A[q][q];
        const aqq = s * s * A[p][p] + 2 * s * c * A[p][q] + c * c * A[q][q];
        const apq = (c * c - s * s) * A[p][q] + s * c * (A[p][p] - A[q][q]);

        [0, 1, 2].forEach(k => {
            if (k !== p && k !== q) {
                const akp = c * A[k][p] - s * A[k][q];
                const akq = s * A[k][p] + c * A[k][q];
                A[k][p] = A[p][k] = akp;
                A[k][q] = A[q][k] = akq;
            }
        });

        A[p][p] = app; A[q][q] = aqq; A[p][q] = A[q][p] = apq;

        for (let i = 0; i < 3; i++) {
            const vip = c * V[i][p] - s * V[i][q];
            const viq = s * V[i][p] + c * V[i][q];
            V[i][p] = vip; V[i][q] = viq;
        }
    }

    eigVals[0] = A[0][0]; eigVals[1] = A[1][1]; eigVals[2] = A[2][2];
    eigVecs[0].set(V[0][0], V[1][0], V[2][0]);
    eigVecs[1].set(V[0][1], V[1][1], V[2][1]);
    eigVecs[2].set(V[0][2], V[1][2], V[2][2]);

    const indices = [0, 1, 2].sort((a, b) => eigVals[a] - eigVals[b]);
    return {
        values: indices.map(i => eigVals[i]),
        vectors: indices.map(i => eigVecs[i]),
        meanNormal: new THREE.Vector3(meanX, meanY, meanZ)
    };
};

/**
 * Classifies a set of faces based on their collective geometric properties.
 */
export const getPatchClassification = (
    geometry: THREE.BufferGeometry,
    faceIndices: number[],
    normals: THREE.Vector3[]
): FeatureType => {
    if (faceIndices.length < 3) return 'generic';

    const pos = geometry.attributes.position;
    const index = geometry.index;
    const patchNormals: THREE.Vector3[] = faceIndices.map(idx => normals[idx]);

    // Perform PCA on Normals
    const pca = computeNormalPCA(patchNormals);
    const [e1, e2, e3] = pca.values; // Sorted ascending
    const axis = pca.vectors[0]; // eigenvector of smallest eigenvalue is the rotation axis

    // 1. PLANE CHECK - Stricter to avoid large cylinders being marked as planes
    if (e3 < 0.005) return 'plane';

    // 2. REVOLUTION SURFACE CHECK
    // 2. REVOLUTION SURFACE CHECK
    // e1 < 0.20 and e2 > 0.01 captures large/noisy industrial cylinders (40mm-65mm).
    if (e1 < 0.20 && e2 > 0.01) {
        // ... Axis finding logic remains the same (already optimized with centroids) ...
        let sXX = 0, sXY = 0, sXZ = 0, sYY = 0, sYZ = 0, sZZ = 0;
        let bX = 0, bY = 0, bZ = 0;

        faceIndices.forEach(idx => {
            const nProj = normals[idx].clone().projectOnPlane(axis).normalize();
            const vA_idx = index ? index.getX(idx * 3 + 0) : idx * 3 + 0;
            const vB_idx = index ? index.getX(idx * 3 + 1) : idx * 3 + 1;
            const vC_idx = index ? index.getX(idx * 3 + 2) : idx * 3 + 2;
            const vA = new THREE.Vector3().fromBufferAttribute(pos, vA_idx);
            const vB = new THREE.Vector3().fromBufferAttribute(pos, vB_idx);
            const vC = new THREE.Vector3().fromBufferAttribute(pos, vC_idx);
            const q = new THREE.Vector3().add(vA).add(vB).add(vC).divideScalar(3);

            const mat = [1 - nProj.x * nProj.x, -nProj.x * nProj.y, -nProj.x * nProj.z, -nProj.y * nProj.x, 1 - nProj.y * nProj.y, -nProj.y * nProj.z, -nProj.z * nProj.x, -nProj.z * nProj.y, 1 - nProj.z * nProj.z];
            sXX += mat[0]; sXY += mat[1]; sXZ += mat[2]; sYY += mat[4]; sYZ += mat[5]; sZZ += mat[8];
            bX += mat[0] * q.x + mat[1] * q.y + mat[2] * q.z; bY += mat[3] * q.x + mat[4] * q.y + mat[5] * q.z; bZ += mat[6] * q.x + mat[7] * q.y + mat[8] * q.z;
        });

        const det = sXX * (sYY * sZZ - sYZ * sYZ) - sXY * (sXY * sZZ - sYZ * sXZ) + sXZ * (sXY * sYZ - sYY * sXZ);
        let bestAxisPoint = new THREE.Vector3(0, 0, 0);
        if (Math.abs(det) > 1e-8) {
            bestAxisPoint.x = ((sYY * sZZ - sYZ * sYZ) * bX - (sXY * sZZ - sYZ * sXZ) * bY + (sXY * sYZ - sYY * sXZ) * bZ) / det;
            bestAxisPoint.y = (-(sXY * sZZ - sYZ * sXZ) * bX + (sXX * sZZ - sXZ * sXZ) * bY - (sXX * sYZ - sXY * sXZ) * bZ) / det;
            bestAxisPoint.z = ((sXY * sYZ - sYY * sXZ) * bX - (sXX * sYZ - sXY * sXZ) * bY + (sXX * sYY - sXY * sXY) * bZ) / det;
        } else {
            faceIndices.forEach(idx => {
                const vA_idx = index ? index.getX(idx * 3 + 0) : idx * 3 + 0;
                bestAxisPoint.add(new THREE.Vector3().fromBufferAttribute(pos, vA_idx));
            });
            bestAxisPoint.divideScalar(faceIndices.length);
        }

        const axisLine = new THREE.Line3(bestAxisPoint, bestAxisPoint.clone().add(axis));
        const heights: number[] = [], radii: number[] = [];
        const p = new THREE.Vector3(), proj = new THREE.Vector3(), vA = new THREE.Vector3(), vB = new THREE.Vector3(), vC = new THREE.Vector3();

        faceIndices.forEach(idx => {
            const vA_idx = index ? index.getX(idx * 3 + 0) : idx * 3 + 0;
            const vB_idx = index ? index.getX(idx * 3 + 1) : idx * 3 + 1;
            const vC_idx = index ? index.getX(idx * 3 + 2) : idx * 3 + 2;
            vA.fromBufferAttribute(pos, vA_idx); vB.fromBufferAttribute(pos, vB_idx); vC.fromBufferAttribute(pos, vC_idx);
            p.set(0, 0, 0).add(vA).add(vB).add(vC).divideScalar(3);
            axisLine.closestPointToPoint(p, false, proj);
            radii.push(p.distanceTo(proj)); heights.push(p.dot(axis));
        });

        const n = radii.length;
        let sH = 0, sR = 0, sHH = 0, sHR = 0;
        for (let i = 0; i < n; i++) {
            sH += heights[i]; sR += radii[i]; sHH += heights[i] * heights[i]; sHR += heights[i] * radii[i];
        }
        const denom = (n * sHH - sH * sH);
        if (Math.abs(denom) > 1e-9) {
            const slope = (n * sHR - sH * sR) / denom;
            const intercept = (sR - slope * sH) / n;
            let totalError = 0;
            const avgR = sR / n;
            for (let i = 0; i < n; i++) totalError += Math.abs(radii[i] - (slope * heights[i] + intercept));
            const normalizedError = totalError / (avgR * n + 1e-6);

            // Slope < 0.01 is a cylinder. Error < 0.3 for noisy meshes.
            if (normalizedError < 0.30 && avgR > 0.01) {
                if (Math.abs(slope) < 0.02) return 'cylinder'; // Slightly more tolerant cylinder
                return 'cone';
            }
        }
    }

    return 'generic';
};

/**
 * Segments the entire mesh into feature patches.
 * Uses a strict visitation pattern to ensure 100% coverage without overlaps.
 */
export const segmentMesh = (geometry: THREE.BufferGeometry, angle: number = 20, adjacency?: number[][]) => {
    const index = geometry.index;
    const faceCount = index ? index.count / 3 : geometry.attributes.position.count / 3;
    const graph = adjacency || buildAdjacencyGraph(geometry);

    // Pre-calculate normals
    const normals: THREE.Vector3[] = [];
    const pos = geometry.attributes.position;
    const pA = new THREE.Vector3(), pB = new THREE.Vector3(), pC = new THREE.Vector3();
    const cb = new THREE.Vector3(), ab = new THREE.Vector3();

    for (let f = 0; f < faceCount; f++) {
        const vA_idx = index ? index.getX(f * 3 + 0) : f * 3 + 0;
        const vB_idx = index ? index.getX(f * 3 + 1) : f * 3 + 1;
        const vC_idx = index ? index.getX(f * 3 + 2) : f * 3 + 2;

        pA.fromBufferAttribute(pos, vA_idx);
        pB.fromBufferAttribute(pos, vB_idx);
        pC.fromBufferAttribute(pos, vC_idx);
        cb.subVectors(pC, pB);
        ab.subVectors(pA, pB);
        cb.cross(ab).normalize();
        normals.push(cb.clone());
    }

    const visited = new Uint8Array(faceCount);
    const initialSegments: { type: FeatureType, indices: number[], axis?: THREE.Vector3, center?: THREE.Vector3, radius?: number }[] = [];
    const thresholdDot = Math.cos(angle * Math.PI / 180);

    for (let i = 0; i < faceCount; i++) {
        if (visited[i]) continue;

        const patch: number[] = [];
        const queue = [i];
        visited[i] = 1;

        while (queue.length > 0) {
            const curr = queue.shift()!;
            patch.push(curr);

            for (const neighbor of graph[curr]) {
                if (!visited[neighbor]) {
                    if (normals[curr].dot(normals[neighbor]) > thresholdDot) {
                        visited[neighbor] = 1;
                        queue.push(neighbor);
                    }
                }
            }
        }

        if (patch.length > 0) {
            const stats = getPatchDetailedStats(geometry, patch, normals);
            initialSegments.push({
                type: stats.type,
                indices: patch,
                axis: stats.axis,
                center: stats.center,
                radius: stats.radius
            });
        }
    }

    // Pass 2: Connected Merge for Compatible Features
    // Only merge segments that are same type, same axis AND PHYSICALLY ADJACENT
    const finalSegments: { type: FeatureType, indices: number[] }[] = [];
    const merged = new Uint8Array(initialSegments.length);

    // Build segment adjacency map
    const segAdjacency = Array.from({ length: initialSegments.length }, () => new Set<number>());
    const faceToSegIdx = new Int32Array(faceCount).fill(-1);
    initialSegments.forEach((seg, sIdx) => {
        seg.indices.forEach(fIdx => faceToSegIdx[fIdx] = sIdx);
    });

    for (let f = 0; f < faceCount; f++) {
        const s1 = faceToSegIdx[f];
        if (s1 === -1) continue;
        for (const neighbor of graph[f]) {
            const s2 = faceToSegIdx[neighbor];
            if (s2 !== -1 && s2 !== s1) {
                segAdjacency[s1].add(s2);
                segAdjacency[s2].add(s1);
            }
        }
    }

    for (let i = 0; i < initialSegments.length; i++) {
        if (merged[i]) continue;

        const clusterIndices = [...initialSegments[i].indices];
        const type = initialSegments[i].type;
        const queue = [i];
        merged[i] = 1;

        while (queue.length > 0) {
            const currIdx = queue.shift()!;
            const segA = initialSegments[currIdx];

            for (const neighborIdx of segAdjacency[currIdx]) {
                if (merged[neighborIdx]) continue;
                const segB = initialSegments[neighborIdx];

                if (type === segB.type && (type === 'cylinder' || type === 'cone')) {
                    if (segA.axis && segB.axis && segA.center && segB.center) {
                        const axisMatch = Math.abs(segA.axis.dot(segB.axis)) > 0.99;
                        const distToAxis = new THREE.Vector3().subVectors(segB.center, segA.center).projectOnPlane(segA.axis).length();
                        const posMatch = distToAxis < 0.8;
                        const radMatch = Math.abs((segA.radius || 0) - (segB.radius || 0)) < 0.5;

                        if (axisMatch && posMatch && (type === 'cone' || radMatch)) {
                            // Merge!
                            clusterIndices.push(...segB.indices);
                            merged[neighborIdx] = 1;
                            queue.push(neighborIdx);
                        }
                    }
                } else if (type === segB.type && type === 'plane') {
                    // Planes are already quite connected in Pass 1, but we can merge adjacent coplanar ones
                    // Actually Pass 1 handled normals, but let's leave planes for now as they are usually fine.
                }
            }
        }
        finalSegments.push({ type, indices: clusterIndices });
    }

    return finalSegments;
};

const getPatchDetailedStats = (geometry: THREE.BufferGeometry, faceIndices: number[], normals: THREE.Vector3[]) => {
    const type = getPatchClassification(geometry, faceIndices, normals);
    if (type === 'generic' || type === 'plane') return { type };

    const pos = geometry.attributes.position;
    const index = geometry.index;
    const patchNormals = faceIndices.map(idx => normals[idx]);
    const pca = computeNormalPCA(patchNormals);
    const axis = pca.vectors[0];

    let sXX = 0, sXY = 0, sXZ = 0, sYY = 0, sYZ = 0, sZZ = 0, bX = 0, bY = 0, bZ = 0;
    faceIndices.forEach(idx => {
        const nProj = normals[idx].clone().projectOnPlane(axis).normalize();
        const vA_idx = index ? index.getX(idx * 3 + 0) : idx * 3 + 0;
        const q = new THREE.Vector3().fromBufferAttribute(pos, vA_idx);
        const mat = [1 - nProj.x * nProj.x, -nProj.x * nProj.y, -nProj.x * nProj.z, -nProj.y * nProj.x, 1 - nProj.y * nProj.y, -nProj.y * nProj.z, -nProj.z * nProj.x, -nProj.z * nProj.y, 1 - nProj.z * nProj.z];
        sXX += mat[0]; sXY += mat[1]; sXZ += mat[2]; sYY += mat[4]; sYZ += mat[5]; sZZ += mat[8];
        bX += mat[0] * q.x + mat[1] * q.y + mat[2] * q.z; bY += mat[3] * q.x + mat[4] * q.y + mat[5] * q.z; bZ += mat[6] * q.x + mat[7] * q.y + mat[8] * q.z;
    });
    const det = sXX * (sYY * sZZ - sYZ * sYZ) - sXY * (sXY * sZZ - sYZ * sXZ) + sXZ * (sXY * sYZ - sYY * sXZ);
    let center = new THREE.Vector3();
    if (Math.abs(det) > 1e-7) {
        center.x = ((sYY * sZZ - sYZ * sYZ) * bX - (sXY * sZZ - sYZ * sXZ) * bY + (sXY * sYZ - sYY * sXZ) * bZ) / det;
        center.y = (-(sXY * sZZ - sYZ * sXZ) * bX + (sXX * sZZ - sXZ * sXZ) * bY - (sXX * sYZ - sXY * sXZ) * bZ) / det;
        center.z = ((sXY * sYZ - sYY * sXZ) * bX - (sXX * sYZ - sXY * sXZ) * bY + (sXX * sYY - sXY * sXY) * bZ) / det;
    } else {
        faceIndices.forEach(idx => {
            const vA_idx = index ? index.getX(idx * 3 + 0) : idx * 3 + 0;
            center.add(new THREE.Vector3().fromBufferAttribute(pos, vA_idx));
        });
        center.divideScalar(faceIndices.length);
    }
    let radiusSum = 0;
    const axisLine = new THREE.Line3(center, center.clone().add(axis));
    const vA = new THREE.Vector3(), vB = new THREE.Vector3(), vC = new THREE.Vector3(), p = new THREE.Vector3();

    faceIndices.forEach(idx => {
        const vA_idx = index ? index.getX(idx * 3 + 0) : idx * 3 + 0;
        const vB_idx = index ? index.getX(idx * 3 + 1) : idx * 3 + 1;
        const vC_idx = index ? index.getX(idx * 3 + 2) : idx * 3 + 2;
        vA.fromBufferAttribute(pos, vA_idx);
        vB.fromBufferAttribute(pos, vB_idx);
        vC.fromBufferAttribute(pos, vC_idx);
        p.set(0, 0, 0).add(vA).add(vB).add(vC).divideScalar(3);

        const proj = new THREE.Vector3();
        axisLine.closestPointToPoint(p, false, proj);
        radiusSum += p.distanceTo(proj);
    });
    return { type, axis, center, radius: radiusSum / faceIndices.length };
};

/**
 * Generates a colorful BufferGeometry representing the classification.
 */
export const generateClassificationGeometry = (
    geometry: THREE.BufferGeometry,
    segments: { type: FeatureType, indices: number[] }[],
    selection: number[] = []
) => {
    const classificationGeometry = new THREE.BufferGeometry();
    const pos = geometry.attributes.position;
    const index = geometry.index;
    const faceCount = index ? index.count / 3 : pos.count / 3;
    const selectionSet = new Set(selection);

    const newPos = new Float32Array(faceCount * 9);
    const newColors = new Float32Array(faceCount * 9);

    const colors: Record<FeatureType, THREE.Color> = {
        'plane': new THREE.Color(0x3498db),    // Blue
        'cylinder': new THREE.Color(0x2ecc71), // Green
        'cone': new THREE.Color(0x9b59b6),     // Purple
        'generic': new THREE.Color(0x95a5a6)   // Gray
    };
    const highlightColor = new THREE.Color(0xffcc00); // Yellow-ish orange selection

    let vertexOffset = 0;
    segments.forEach(seg => {
        const baseColor = colors[seg.type];
        seg.indices.forEach(faceIdx => {
            const isSelected = selectionSet.has(faceIdx);
            const color = isSelected ? highlightColor : baseColor;

            for (let v = 0; v < 3; v++) {
                const vIdx = index ? index.getX(faceIdx * 3 + v) : faceIdx * 3 + v;
                const outIdx = vertexOffset * 3;
                newPos[outIdx + 0] = pos.getX(vIdx);
                newPos[outIdx + 1] = pos.getY(vIdx);
                newPos[outIdx + 2] = pos.getZ(vIdx);

                newColors[outIdx + 0] = color.r;
                newColors[outIdx + 1] = color.g;
                newColors[outIdx + 2] = color.b;
                vertexOffset++;
            }
        });
    });

    classificationGeometry.setAttribute('position', new THREE.BufferAttribute(newPos, 3));
    classificationGeometry.setAttribute('color', new THREE.BufferAttribute(newColors, 3));
    classificationGeometry.computeVertexNormals();
    return classificationGeometry;
};

// Generate a BufferGeometry for the selected triangles to render highlight
export const generateSelectionGeometry = (originalGeometry: THREE.BufferGeometry, faceIndices: number[]) => {
    const selectedGeometry = new THREE.BufferGeometry();
    const pos = originalGeometry.attributes.position;
    const index = originalGeometry.index;

    const newPos = new Float32Array(faceIndices.length * 9); // 3 vertices * 3 coords

    faceIndices.forEach((faceIdx, i) => {
        const vA_idx = index ? index.getX(faceIdx * 3 + 0) : faceIdx * 3 + 0;
        const vB_idx = index ? index.getX(faceIdx * 3 + 1) : faceIdx * 3 + 1;
        const vC_idx = index ? index.getX(faceIdx * 3 + 2) : faceIdx * 3 + 2;

        newPos[i * 9 + 0] = pos.getX(vA_idx);
        newPos[i * 9 + 1] = pos.getY(vA_idx);
        newPos[i * 9 + 2] = pos.getZ(vA_idx);

        newPos[i * 9 + 3] = pos.getX(vB_idx);
        newPos[i * 9 + 4] = pos.getY(vB_idx);
        newPos[i * 9 + 5] = pos.getZ(vB_idx);

        newPos[i * 9 + 6] = pos.getX(vC_idx);
        newPos[i * 9 + 7] = pos.getY(vC_idx);
        newPos[i * 9 + 8] = pos.getZ(vC_idx);
    });

    selectedGeometry.setAttribute('position', new THREE.BufferAttribute(newPos, 3));
    return selectedGeometry;
};
