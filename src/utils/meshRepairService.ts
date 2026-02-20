import * as THREE from 'three';

/**
 * Pure JavaScript mesh repair service.
 * Refined to avoid filling large functional holes while still fixing small tears.
 */

interface IndexedMesh {
    vertices: number[];   // flat xyz
    faces: number[];      // flat v0,v1,v2
}

function edgeKey(a: number, b: number): string {
    return a < b ? `${a}_${b}` : `${b}_${a}`;
}

function buildIndexedMesh(positions: ArrayLike<number>): IndexedMesh {
    const numFaces = (positions.length / 3) / 3;
    const PREC = 100000;
    const vertexMap = new Map<string, number>();
    const uniqueVerts: number[] = [];
    const faceIndices: number[] = new Array(numFaces * 3);

    let numUniqueVerts = 0;
    for (let fi = 0; fi < numFaces; fi++) {
        for (let vi = 0; vi < 3; vi++) {
            const srcIdx = (fi * 3 + vi) * 3;
            const x = positions[srcIdx], y = positions[srcIdx + 1], z = positions[srcIdx + 2];
            const key = `${Math.round(x * PREC)},${Math.round(y * PREC)},${Math.round(z * PREC)}`;
            let idx = vertexMap.get(key);
            if (idx === undefined) {
                idx = numUniqueVerts++;
                vertexMap.set(key, idx);
                uniqueVerts.push(x, y, z);
            }
            faceIndices[fi * 3 + vi] = idx;
        }
    }

    return {
        vertices: uniqueVerts,
        faces: faceIndices
    };
}

function snapNearbyVertices(mesh: IndexedMesh, tolerance: number): number {
    const numFaces = mesh.faces.length / 3;
    const edgeCount = new Map<string, number>();
    for (let fi = 0; fi < numFaces; fi++) {
        const o = fi * 3;
        for (const [a, b] of [[mesh.faces[o], mesh.faces[o + 1]], [mesh.faces[o + 1], mesh.faces[o + 2]], [mesh.faces[o + 2], mesh.faces[o]]]) {
            if (a === b) continue;
            const key = edgeKey(a, b);
            edgeCount.set(key, (edgeCount.get(key) || 0) + 1);
        }
    }

    const boundaryVerts = new Set<number>();
    for (const [key, count] of edgeCount) {
        if (count === 1) {
            const [a, b] = key.split('_').map(Number);
            boundaryVerts.add(a); boundaryVerts.add(b);
        }
    }

    if (boundaryVerts.size === 0) return 0;

    const cellSize = tolerance * 2;
    const grid = new Map<string, number[]>();
    const verts = mesh.vertices;

    for (const vi of boundaryVerts) {
        const cx = Math.floor(verts[vi * 3] / cellSize), cy = Math.floor(verts[vi * 3 + 1] / cellSize), cz = Math.floor(verts[vi * 3 + 2] / cellSize);
        const cellKey = `${cx},${cy},${cz}`;
        let cell = grid.get(cellKey);
        if (!cell) { cell = []; grid.set(cellKey, cell); }
        cell.push(vi);
    }

    const numUnique = verts.length / 3;
    const mergeMap = new Int32Array(numUnique);
    for (let i = 0; i < numUnique; i++) mergeMap[i] = i;

    function root(v: number): number {
        while (mergeMap[v] !== v) v = mergeMap[v];
        return v;
    }

    const tolSq = tolerance * tolerance;
    let mergeCount = 0;

    for (const vi of boundaryVerts) {
        const x = verts[vi * 3], y = verts[vi * 3 + 1], z = verts[vi * 3 + 2];
        const cx = Math.floor(x / cellSize), cy = Math.floor(y / cellSize), cz = Math.floor(z / cellSize);

        for (let dx = -1; dx <= 1; dx++) {
            for (let dy = -1; dy <= 1; dy++) {
                for (let dz = -1; dz <= 1; dz++) {
                    const cell = grid.get(`${cx + dx},${cy + dy},${cz + dz}`);
                    if (!cell) continue;
                    for (const vj of cell) {
                        if (vj <= vi) continue;
                        const rv = root(vi), rj = root(vj);
                        if (rv === rj) continue;
                        const distSq = (verts[vj * 3] - x) ** 2 + (verts[vj * 3 + 1] - y) ** 2 + (verts[vj * 3 + 2] - z) ** 2;
                        if (distSq <= tolSq) {
                            mergeMap[rj] = rv;
                            verts[vj * 3] = x; verts[vj * 3 + 1] = y; verts[vj * 3 + 2] = z;
                            mergeCount++;
                        }
                    }
                }
            }
        }
    }

    if (mergeCount === 0) return 0;
    for (let i = 0; i < mesh.faces.length; i++) mesh.faces[i] = root(mesh.faces[i]);
    return mergeCount;
}

function fillHoles(mesh: IndexedMesh, maxHoleEdges: number, maxInradius: number = Infinity): number {
    const numFaces = mesh.faces.length / 3;
    const edgeFaceCount = new Map<string, number>();
    for (let fi = 0; fi < numFaces; fi++) {
        const o = fi * 3;
        for (const [a, b] of [[mesh.faces[o], mesh.faces[o + 1]], [mesh.faces[o + 1], mesh.faces[o + 2]], [mesh.faces[o + 2], mesh.faces[o]]]) {
            const key = edgeKey(a, b);
            edgeFaceCount.set(key, (edgeFaceCount.get(key) || 0) + 1);
        }
    }

    const adj = new Map<number, number[]>();
    const presentDirectedEdges = new Set<string>();
    for (let fi = 0; fi < numFaces; fi++) {
        const o = fi * 3;
        presentDirectedEdges.add(`${mesh.faces[o]}_${mesh.faces[o + 1]}`);
        presentDirectedEdges.add(`${mesh.faces[o + 1]}_${mesh.faces[o + 2]}`);
        presentDirectedEdges.add(`${mesh.faces[o + 2]}_${mesh.faces[o]}`);
    }

    for (const [key, count] of edgeFaceCount) {
        if (count === 1) {
            const [u, v] = key.split('_').map(Number);
            if (presentDirectedEdges.has(`${u}_${v}`)) {
                if (!adj.has(v)) adj.set(v, []); adj.get(v)!.push(u);
            } else {
                if (!adj.has(u)) adj.set(u, []); adj.get(u)!.push(v);
            }
        }
    }

    const visitedEdges = new Set<string>();
    let filledCount = 0;

    for (const startNode of adj.keys()) {
        const neighbors = adj.get(startNode);
        if (!neighbors) continue;
        for (const nextNode of neighbors) {
            const edgeId = `${startNode}_${nextNode}`;
            if (visitedEdges.has(edgeId)) continue;

            const loop: number[] = [startNode];
            let curr = nextNode;
            visitedEdges.add(edgeId);
            let stuck = false;
            while (curr !== startNode) {
                loop.push(curr);
                if (loop.length > maxHoleEdges) { stuck = true; break; }
                const nexts = adj.get(curr);
                if (!nexts) { stuck = true; break; }
                let foundNext = -1;
                for (const candidate of nexts) {
                    if (!visitedEdges.has(`${curr}_${candidate}`)) { foundNext = candidate; break; }
                }
                if (foundNext === -1) { stuck = true; break; }
                visitedEdges.add(`${curr}_${foundNext}`);
                curr = foundNext;
            }

            if (!stuck && loop.length >= 3) {
                // HEURISTIC: Skip if this is a large functional opening (like a pipe end)
                if (maxInradius !== Infinity) {
                    let cx = 0, cy = 0, cz = 0;
                    for (const vi of loop) { cx += mesh.vertices[vi * 3]; cy += mesh.vertices[vi * 3 + 1]; cz += mesh.vertices[vi * 3 + 2]; }
                    cx /= loop.length; cy /= loop.length; cz /= loop.length;

                    let minDistSq = Infinity;
                    for (const vi of loop) {
                        const d2 = (mesh.vertices[vi * 3] - cx) ** 2 + (mesh.vertices[vi * 3 + 1] - cy) ** 2 + (mesh.vertices[vi * 3 + 2] - cz) ** 2;
                        if (d2 < minDistSq) minDistSq = d2;
                    }
                    // If the 'inradius' is much larger than our tolerance, it's a design hole, not a tear.
                    if (Math.sqrt(minDistSq) > maxInradius * 1.5) {
                        continue;
                    }
                }

                triangulatePolygonRobust(mesh, loop);
                filledCount++;
            }
        }
    }
    return filledCount;
}

function triangulatePolygonRobust(mesh: IndexedMesh, loopIndices: number[]) {
    const polygon = [...loopIndices];
    const verts = mesh.vertices;
    let safety = 0;
    while (polygon.length > 3 && safety < 2000) {
        safety++;
        let bestEarIdx = -1;
        let maxMinAngle = -Infinity;
        for (let i = 0; i < polygon.length; i++) {
            const iPrev = (i - 1 + polygon.length) % polygon.length;
            const iNext = (i + 1) % polygon.length;
            const p = polygon[i], a = polygon[iPrev], b = polygon[iNext];
            const va = new THREE.Vector3(verts[a * 3] - verts[p * 3], verts[a * 3 + 1] - verts[p * 3 + 1], verts[a * 3 + 2] - verts[p * 3 + 2]);
            const vb = new THREE.Vector3(verts[b * 3] - verts[p * 3], verts[b * 3 + 1] - verts[p * 3 + 1], verts[b * 3 + 2] - verts[p * 3 + 2]);
            const vab = new THREE.Vector3(verts[b * 3] - verts[a * 3], verts[b * 3 + 1] - verts[a * 3 + 1], verts[b * 3 + 2] - verts[a * 3 + 2]);
            const la = va.length(), lb = vb.length(), lab = vab.length();
            if (la < 1e-9 || lb < 1e-9 || lab < 1e-9) { bestEarIdx = i; break; }
            const cosP = va.dot(vb) / (la * lb);
            const minAngleCos = Math.max(cosP, va.clone().negate().dot(vab) / (la * lab), vb.clone().negate().dot(vab.clone().negate()) / (lb * lab));
            let isEar = true;
            for (let j = 0; j < polygon.length; j++) {
                if (j === i || j === iPrev || j === iNext) continue;
                if (isPointInTriangle(polygon[j], a, p, b, verts)) { isEar = false; break; }
            }
            if (isEar) {
                const quality = -minAngleCos;
                if (quality > maxMinAngle) { maxMinAngle = quality; bestEarIdx = i; }
            }
        }
        if (bestEarIdx !== -1) {
            const iPrev = (bestEarIdx - 1 + polygon.length) % polygon.length;
            const iNext = (bestEarIdx + 1) % polygon.length;
            mesh.faces.push(polygon[iPrev], polygon[bestEarIdx], polygon[iNext]);
            polygon.splice(bestEarIdx, 1);
        } else {
            mesh.faces.push(polygon[0], polygon[1], polygon[2]);
            polygon.splice(1, 1);
        }
    }
    if (polygon.length === 3) mesh.faces.push(polygon[0], polygon[1], polygon[2]);
}

function isPointInTriangle(pIdx: number, aIdx: number, bIdx: number, cIdx: number, verts: number[]): boolean {
    const P = new THREE.Vector3(verts[pIdx * 3], verts[pIdx * 3 + 1], verts[pIdx * 3 + 2]);
    const A = new THREE.Vector3(verts[aIdx * 3], verts[aIdx * 3 + 1], verts[aIdx * 3 + 2]);
    const B = new THREE.Vector3(verts[bIdx * 3], verts[bIdx * 3 + 1], verts[bIdx * 3 + 2]);
    const C = new THREE.Vector3(verts[cIdx * 3], verts[cIdx * 3 + 1], verts[cIdx * 3 + 2]);
    const v0 = new THREE.Vector3().subVectors(C, A), v1 = new THREE.Vector3().subVectors(B, A), v2 = new THREE.Vector3().subVectors(P, A);
    const dot00 = v0.dot(v0), dot01 = v0.dot(v1), dot02 = v0.dot(v2), dot11 = v1.dot(v1), dot12 = v1.dot(v2);
    const invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    const u = (dot11 * dot02 - dot01 * dot12) * invDenom, v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    return (u >= 0) && (v >= 0) && (u + v < 1);
}

function removeDegenerateFaces(mesh: IndexedMesh): IndexedMesh {
    const valid: number[] = [];
    for (let i = 0; i < mesh.faces.length / 3; i++) {
        const a = mesh.faces[i * 3], b = mesh.faces[i * 3 + 1], c = mesh.faces[i * 3 + 2];
        if (a !== b && b !== c && c !== a) valid.push(a, b, c);
    }
    mesh.faces = valid; return mesh;
}

function removeUnconnectedFacets(mesh: IndexedMesh): IndexedMesh {
    const edgeCount = new Map<string, number>();
    for (let i = 0; i < mesh.faces.length / 3; i++) {
        const v0 = mesh.faces[i * 3], v1 = mesh.faces[i * 3 + 1], v2 = mesh.faces[i * 3 + 2];
        for (const [a, b] of [[v0, v1], [v1, v2], [v2, v0]]) {
            const key = edgeKey(a, b);
            edgeCount.set(key, (edgeCount.get(key) || 0) + 1);
        }
    }
    const keep: number[] = [];
    for (let i = 0; i < mesh.faces.length / 3; i++) {
        const v0 = mesh.faces[i * 3], v1 = mesh.faces[i * 3 + 1], v2 = mesh.faces[i * 3 + 2];
        let conn = 0;
        for (const [a, b] of [[v0, v1], [v1, v2], [v2, v0]]) { if ((edgeCount.get(edgeKey(a, b)) || 0) >= 2) conn++; }
        if (conn >= 1) keep.push(v0, v1, v2);
    }
    mesh.faces = keep; return mesh;
}

function fixNormalDirections(mesh: IndexedMesh): void {
    const numFaces = mesh.faces.length / 3;
    const edgeToFaces = new Map<string, number[]>();
    for (let i = 0; i < numFaces; i++) {
        const v0 = mesh.faces[i * 3], v1 = mesh.faces[i * 3 + 1], v2 = mesh.faces[i * 3 + 2];
        for (const [a, b] of [[v0, v1], [v1, v2], [v2, v0]]) {
            const key = edgeKey(a, b);
            let list = edgeToFaces.get(key); if (!list) { list = []; edgeToFaces.set(key, list); }
            list.push(i);
        }
    }
    const visited = new Uint8Array(numFaces), queue: number[] = [];
    function getDir(f: number, a: number, b: number) {
        const v0 = mesh.faces[f * 3], v1 = mesh.faces[f * 3 + 1], v2 = mesh.faces[f * 3 + 2];
        return (v0 === a && v1 === b) || (v1 === a && v2 === b) || (v2 === a && v0 === b);
    }
    for (let s = 0; s < numFaces; s++) {
        if (visited[s]) continue; queue.push(s); visited[s] = 1;
        while (queue.length > 0) {
            const f = queue.shift()!;
            for (const [a, b] of [[mesh.faces[f * 3], mesh.faces[f * 3 + 1]], [mesh.faces[f * 3 + 1], mesh.faces[f * 3 + 2]], [mesh.faces[f * 3 + 2], mesh.faces[f * 3]]]) {
                const neighbors = edgeToFaces.get(edgeKey(a, b)); if (!neighbors) continue;
                for (const nf of neighbors) {
                    if (visited[nf]) continue; visited[nf] = 1;
                    if (getDir(f, a, b) === getDir(nf, a, b)) {
                        const tmp = mesh.faces[nf * 3 + 1]; mesh.faces[nf * 3 + 1] = mesh.faces[nf * 3 + 2]; mesh.faces[nf * 3 + 2] = tmp;
                    }
                    queue.push(nf);
                }
            }
        }
    }
}

function fixVolumeSign(mesh: IndexedMesh): void {
    const v = mesh.vertices; let vol = 0;
    for (let i = 0; i < mesh.faces.length / 3; i++) {
        const a = mesh.faces[i * 3] * 3, b = mesh.faces[i * 3 + 1] * 3, c = mesh.faces[i * 3 + 2] * 3;
        vol += v[a] * (v[b + 1] * v[c + 2] - v[b + 2] * v[c + 1]) + v[a + 1] * (v[b + 2] * v[c] - v[b] * v[c + 2]) + v[a + 2] * (v[b] * v[c + 1] - v[b + 1] * v[c]);
    }
    if (vol < 0) {
        for (let i = 0; i < mesh.faces.length / 3; i++) {
            const o = i * 3, tmp = mesh.faces[o + 1];
            mesh.faces[o + 1] = mesh.faces[o + 2]; mesh.faces[o + 2] = tmp;
        }
    }
}

export async function repairGeometry(
    geometry: THREE.BufferGeometry,
    maxTolerance?: number
): Promise<THREE.BufferGeometry> {
    try {
        let mesh = buildIndexedMesh(geometry.attributes.position.array);
        const tolerances = [0.001, 0.01, 0.1, 0.5];
        if (maxTolerance && maxTolerance > 0.5) tolerances.push(maxTolerance);
        for (const t of tolerances) if (snapNearbyVertices(mesh, t) > 0) removeDegenerateFaces(mesh);
        removeUnconnectedFacets(mesh);

        // Use threshold to intelligently skip large openings
        fillHoles(mesh, 2000, maxTolerance || Infinity);

        fixNormalDirections(mesh);
        fixVolumeSign(mesh);
        const result = new THREE.BufferGeometry();
        const pos = new Float32Array(mesh.faces.length * 3);
        const verts = mesh.vertices;
        for (let i = 0; i < mesh.faces.length; i++) {
            const vi = mesh.faces[i] * 3;
            pos[i * 3] = verts[vi]; pos[i * 3 + 1] = verts[vi + 1]; pos[i * 3 + 2] = verts[vi + 2];
        }
        result.setAttribute('position', new THREE.BufferAttribute(pos, 3));
        result.computeVertexNormals();
        return result;
    } catch (err) {
        console.error('[MeshRepair] Failed:', err); return geometry;
    }
}
