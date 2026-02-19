import * as THREE from 'three';

/**
 * Pure JavaScript mesh repair service.
 * Pipeline:
 * 1. Build indexed mesh (Exact deduplication)
 * 2. Remove degenerate faces
 * 3. Stitch: Snap nearby vertices (Tolerance-based)
 * 4. Remove unconnected facets (Cleanup orphans)
 * 5. Fill Holes: Close boundary loops (Advanced Ear Clipping)
 * 6. Fix normal directions & Volume sign
 */

interface IndexedMesh {
    vertices: Float32Array;   // flat xyz
    faces: number[];          // flat v0,v1,v2
    numVertices: number;
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
        vertices: new Float32Array(uniqueVerts),
        faces: faceIndices,
        numVertices: numUniqueVerts
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

    const mergeMap = new Int32Array(mesh.numVertices);
    for (let i = 0; i < mesh.numVertices; i++) mergeMap[i] = i;

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

function fillHoles(mesh: IndexedMesh, maxHoleEdges: number): number {
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

    const visited = new Set<string>();
    let filledCount = 0;

    for (const startVal of adj.keys()) {
        const neighbors = adj.get(startVal);
        if (!neighbors) continue;
        for (const nextVal of neighbors) {
            const loopKey = `${startVal}_${nextVal}`;
            if (visited.has(loopKey)) continue;

            const loop: number[] = [startVal];
            let curr = nextVal;
            visited.add(loopKey);
            let stuck = false;
            while (curr !== startVal) {
                loop.push(curr);
                if (loop.length > maxHoleEdges) { stuck = true; break; }
                const nexts = adj.get(curr);
                if (!nexts) { stuck = true; break; }
                let foundNext = -1;
                for (const n of nexts) { if (!visited.has(`${curr}_${n}`)) { foundNext = n; break; } }
                if (foundNext === -1) { stuck = true; break; }
                visited.add(`${curr}_${foundNext}`);
                curr = foundNext;
            }

            if (!stuck && loop.length >= 3) {
                triangulatePolygon(mesh, loop, mesh.vertices);
                filledCount++;
            }
        }
    }
    return filledCount;
}

function triangulatePolygon(mesh: IndexedMesh, loopIndices: number[], verts: Float32Array) {
    const polygon = loopIndices.slice();
    let safety = 0;
    while (polygon.length > 3 && safety < 1000) {
        safety++;
        let bestIdx = -1, minScore = Infinity;
        for (let i = 0; i < polygon.length; i++) {
            const prev = polygon[(i - 1 + polygon.length) % polygon.length];
            const next = polygon[(i + 1) % polygon.length];
            const distSq = (verts[prev * 3] - verts[next * 3]) ** 2 + (verts[prev * 3 + 1] - verts[next * 3 + 1]) ** 2 + (verts[prev * 3 + 2] - verts[next * 3 + 2]) ** 2;
            if (distSq < minScore) { minScore = distSq; bestIdx = i; }
        }
        if (bestIdx !== -1) {
            mesh.faces.push(polygon[(bestIdx - 1 + polygon.length) % polygon.length], polygon[bestIdx], polygon[(bestIdx + 1) % polygon.length]);
            polygon.splice(bestIdx, 1);
        } else break;
    }
    if (polygon.length === 3) mesh.faces.push(polygon[0], polygon[1], polygon[2]);
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

export async function repairGeometry(geometry: THREE.BufferGeometry): Promise<THREE.BufferGeometry> {
    try {
        let mesh = buildIndexedMesh(geometry.attributes.position.array);
        const tolerances = [0.001, 0.01, 0.1, 0.5];
        for (const t of tolerances) if (snapNearbyVertices(mesh, t) > 0) removeDegenerateFaces(mesh);
        removeUnconnectedFacets(mesh);
        fillHoles(mesh, 2000);
        fixNormalDirections(mesh);
        fixVolumeSign(mesh);
        const result = new THREE.BufferGeometry();
        const pos = new Float32Array(mesh.faces.length * 3);
        for (let i = 0; i < mesh.faces.length; i++) {
            const vi = mesh.faces[i] * 3;
            pos[i * 3] = mesh.vertices[vi]; pos[i * 3 + 1] = mesh.vertices[vi + 1]; pos[i * 3 + 2] = mesh.vertices[vi + 2];
        }
        result.setAttribute('position', new THREE.BufferAttribute(pos, 3));
        result.computeVertexNormals();
        return result;
    } catch (err) {
        console.error('[MeshRepair] Failed:', err); return geometry;
    }
}
