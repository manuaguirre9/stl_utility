import * as THREE from 'three';
import initManifold from 'manifold-3d';
import type { ManifoldToplevel } from 'manifold-3d';

/**
 * Conservative Mesh Repair Service
 * 
 * Philosophy: ONLY fill small gaps between surfaces to make the mesh watertight.
 * - Never remove existing faces or components
 * - Never reorient existing faces
 * - Only add very small triangles to fill boundary holes
 * - Skip large holes (they might be intentional or too complex)
 * - The result should be visually identical to the original
 */

let wasm: ManifoldToplevel | null = null;

async function getWasm(): Promise<ManifoldToplevel> {
    if (!wasm) {
        wasm = await initManifold({ locateFile: () => '/manifold.wasm' });
        if (wasm.setup) wasm.setup();
    }
    return wasm;
}

// ─── Types ───────────────────────────────────────────────────────────

interface IndexedMesh {
    vertices: number[];   // flat xyz
    faces: number[];      // flat v0, v1, v2
}

// ─── Utilities ───────────────────────────────────────────────────────

function edgeKey(a: number, b: number): string {
    return a < b ? `${a}_${b}` : `${b}_${a}`;
}

function vertexDist(mesh: IndexedMesh, a: number, b: number): number {
    const ax = mesh.vertices[a * 3], ay = mesh.vertices[a * 3 + 1], az = mesh.vertices[a * 3 + 2];
    const bx = mesh.vertices[b * 3], by = mesh.vertices[b * 3 + 1], bz = mesh.vertices[b * 3 + 2];
    const dx = bx - ax, dy = by - ay, dz = bz - az;
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/** Compute average edge length of the mesh (for size thresholds). */
function computeAvgEdgeLength(mesh: IndexedMesh): number {
    let totalLen = 0;
    let count = 0;
    const numFaces = mesh.faces.length / 3;
    for (let fi = 0; fi < numFaces; fi++) {
        const a = mesh.faces[fi * 3], b = mesh.faces[fi * 3 + 1], c = mesh.faces[fi * 3 + 2];
        totalLen += vertexDist(mesh, a, b);
        totalLen += vertexDist(mesh, b, c);
        totalLen += vertexDist(mesh, c, a);
        count += 3;
    }
    return count > 0 ? totalLen / count : 1;
}

// ─── Pipeline stages ─────────────────────────────────────────────────

/**
 * Build indexed mesh from flat position array, welding coincident vertices.
 */
function buildIndexedMesh(positions: ArrayLike<number>): IndexedMesh {
    const PREC = 1e5;
    const vertexMap = new Map<string, number>();
    const vertices: number[] = [];
    const faces: number[] = [];
    let vertCount = 0;

    const numVerts = positions.length / 3;
    for (let i = 0; i < numVerts; i++) {
        const x = positions[i * 3], y = positions[i * 3 + 1], z = positions[i * 3 + 2];
        const key = `${Math.round(x * PREC)},${Math.round(y * PREC)},${Math.round(z * PREC)}`;
        let idx = vertexMap.get(key);
        if (idx === undefined) {
            idx = vertCount++;
            vertexMap.set(key, idx);
            vertices.push(x, y, z);
        }
        faces.push(idx);
    }
    return { vertices, faces };
}

/**
 * Remove only truly degenerate faces (where 2+ vertex INDICES are the same).
 * Does NOT remove any face that has valid geometry.
 */
function removeDegenerateFaces(mesh: IndexedMesh): number {
    const clean: number[] = [];
    let removed = 0;
    for (let i = 0; i < mesh.faces.length; i += 3) {
        const a = mesh.faces[i], b = mesh.faces[i + 1], c = mesh.faces[i + 2];
        if (a !== b && b !== c && c !== a) {
            clean.push(a, b, c);
        } else {
            removed++;
        }
    }
    mesh.faces = clean;
    return removed;
}

/**
 * Find all boundary edges (edges shared by exactly 1 face).
 * Returns a map of directed boundary edges for hole tracing.
 */
function findBoundaryEdges(mesh: IndexedMesh): {
    boundaryCount: number;
    /** Map from vertex A → list of vertex B where A→B is a boundary direction */
    adj: Map<number, number[]>;
} {
    const numFaces = mesh.faces.length / 3;
    const edgeFaceCount = new Map<string, number>();
    const directedEdges = new Set<string>();

    for (let fi = 0; fi < numFaces; fi++) {
        const v0 = mesh.faces[fi * 3], v1 = mesh.faces[fi * 3 + 1], v2 = mesh.faces[fi * 3 + 2];
        directedEdges.add(`${v0}_${v1}`);
        directedEdges.add(`${v1}_${v2}`);
        directedEdges.add(`${v2}_${v0}`);
        for (const [a, b] of [[v0, v1], [v1, v2], [v2, v0]]) {
            const k = edgeKey(a, b);
            edgeFaceCount.set(k, (edgeFaceCount.get(k) || 0) + 1);
        }
    }

    // Build adjacency for boundary edges only, in the REVERSE direction
    // (so traversing follows the hole boundary consistently)
    const adj = new Map<number, number[]>();
    let boundaryCount = 0;

    for (const [key, count] of edgeFaceCount) {
        if (count !== 1) continue;
        boundaryCount++;
        const [uStr, vStr] = key.split('_');
        const u = Number(uStr), v = Number(vStr);

        // Find which direction this edge exists as a directed edge
        if (directedEdges.has(`${u}_${v}`)) {
            // The face has u→v, so the hole boundary goes v→u
            let list = adj.get(v);
            if (!list) { list = []; adj.set(v, list); }
            list.push(u);
        } else {
            // The face has v→u, so the hole boundary goes u→v
            let list = adj.get(u);
            if (!list) { list = []; adj.set(u, list); }
            list.push(v);
        }
    }

    return { boundaryCount, adj };
}

/**
 * Trace boundary loops from the adjacency map.
 * Each loop is a list of vertex indices forming a closed polygon (the hole boundary).
 */
function traceBoundaryLoops(adj: Map<number, number[]>, maxEdges: number = 2000): number[][] {
    const visitedEdges = new Set<string>();
    const loops: number[][] = [];

    for (const startNode of adj.keys()) {
        const neighbors = adj.get(startNode);
        if (!neighbors) continue;

        for (const nextNode of neighbors) {
            const startEdgeKey = `${startNode}_${nextNode}`;
            if (visitedEdges.has(startEdgeKey)) continue;

            // Trace a loop
            const loop: number[] = [startNode];
            let curr = nextNode;
            visitedEdges.add(startEdgeKey);
            let stuck = false;

            while (curr !== startNode) {
                loop.push(curr);
                if (loop.length > maxEdges) { stuck = true; break; }
                const nexts = adj.get(curr);
                if (!nexts) { stuck = true; break; }

                let foundNext = -1;
                for (const candidate of nexts) {
                    if (!visitedEdges.has(`${curr}_${candidate}`)) {
                        foundNext = candidate;
                        break;
                    }
                }
                if (foundNext === -1) { stuck = true; break; }
                visitedEdges.add(`${curr}_${foundNext}`);
                curr = foundNext;
            }

            if (!stuck && loop.length >= 3) {
                loops.push(loop);
            }
        }
    }

    return loops;
}

/**
 * Compute the maximum edge length in a boundary loop.
 */
function loopMaxEdgeLength(mesh: IndexedMesh, loop: number[]): number {
    let maxLen = 0;
    for (let i = 0; i < loop.length; i++) {
        const next = (i + 1) % loop.length;
        const d = vertexDist(mesh, loop[i], loop[next]);
        if (d > maxLen) maxLen = d;
    }
    return maxLen;
}

/**
 * Compute the "diameter" of a boundary loop (maximum distance between any two vertices).
 * We approximate by checking start vs all others (O(n) instead of O(n²)).
 */
function loopDiameter(mesh: IndexedMesh, loop: number[]): number {
    if (loop.length < 2) return 0;
    let maxDist = 0;
    // Check from first vertex and from the midpoint vertex
    const checkFrom = [0, Math.floor(loop.length / 2)];
    for (const fromIdx of checkFrom) {
        const from = loop[fromIdx];
        for (let i = 0; i < loop.length; i++) {
            const d = vertexDist(mesh, from, loop[i]);
            if (d > maxDist) maxDist = d;
        }
    }
    return maxDist;
}

/**
 * Compute best-fit plane normal via Newell's method, then
 * return RMS deviation of boundary vertices from that plane.
 * Returns a ratio = RMS_deviation / hole_diameter.
 */
function holePlanarityRatio(mesh: IndexedMesh, loop: number[]): number {
    if (loop.length < 3) return 1;
    // Centroid
    let cx = 0, cy = 0, cz = 0;
    for (const vi of loop) {
        cx += mesh.vertices[vi * 3];
        cy += mesh.vertices[vi * 3 + 1];
        cz += mesh.vertices[vi * 3 + 2];
    }
    cx /= loop.length; cy /= loop.length; cz /= loop.length;

    // Newell normal
    let nx = 0, ny = 0, nz = 0;
    for (let i = 0; i < loop.length; i++) {
        const a = loop[i], b = loop[(i + 1) % loop.length];
        const ax = mesh.vertices[a * 3] - cx, ay = mesh.vertices[a * 3 + 1] - cy, az = mesh.vertices[a * 3 + 2] - cz;
        const bx = mesh.vertices[b * 3] - cx, by = mesh.vertices[b * 3 + 1] - cy, bz = mesh.vertices[b * 3 + 2] - cz;
        nx += ay * bz - az * by;
        ny += az * bx - ax * bz;
        nz += ax * by - ay * bx;
    }
    const nl = Math.sqrt(nx * nx + ny * ny + nz * nz);
    if (nl < 1e-10) return 1;
    nx /= nl; ny /= nl; nz /= nl;

    // RMS deviation from plane
    let sumSq = 0;
    let maxDist = 0;
    for (const vi of loop) {
        const dx = mesh.vertices[vi * 3] - cx;
        const dy = mesh.vertices[vi * 3 + 1] - cy;
        const dz = mesh.vertices[vi * 3 + 2] - cz;
        const d = Math.abs(dx * nx + dy * ny + dz * nz);
        sumSq += d * d;
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (r > maxDist) maxDist = r;
    }
    const rms = Math.sqrt(sumSq / loop.length);
    return maxDist > 0 ? rms / maxDist : 0;
}

/**
 * Fan-fill: add a centroid vertex, then create N triangles connecting
 * each boundary edge to the centroid. Works perfectly for flat holes
 * like missing cylinder caps — creates uniform, small triangles.
 */
function fillHoleWithFan(mesh: IndexedMesh, loop: number[]): number {
    if (loop.length < 3) return 0;

    // Compute centroid
    let cx = 0, cy = 0, cz = 0;
    for (const vi of loop) {
        cx += mesh.vertices[vi * 3];
        cy += mesh.vertices[vi * 3 + 1];
        cz += mesh.vertices[vi * 3 + 2];
    }
    cx /= loop.length; cy /= loop.length; cz /= loop.length;

    // Add centroid as a new vertex
    const centroidIdx = mesh.vertices.length / 3;
    mesh.vertices.push(cx, cy, cz);

    // Create a triangle for each boundary edge
    let added = 0;
    for (let i = 0; i < loop.length; i++) {
        const a = loop[i];
        const b = loop[(i + 1) % loop.length];
        mesh.faces.push(a, b, centroidIdx);
        added++;
    }
    return added;
}

/**
 * Ear-clipping triangulation for non-planar or concave holes.
 * Greedily picks the ear with the smallest longest edge at each step.
 */
function fillHoleEarClip(mesh: IndexedMesh, loop: number[]): number {
    if (loop.length < 3) return 0;
    const polygon = [...loop];
    let addedFaces = 0;
    let safety = 0;

    while (polygon.length > 3 && safety < 5000) {
        safety++;
        let bestEarIdx = -1;
        let bestScore = Infinity;

        for (let i = 0; i < polygon.length; i++) {
            const p = polygon[(i - 1 + polygon.length) % polygon.length];
            const ear = polygon[i];
            const n = polygon[(i + 1) % polygon.length];
            const maxEdge = Math.max(
                vertexDist(mesh, p, ear),
                vertexDist(mesh, ear, n),
                vertexDist(mesh, p, n)
            );
            if (maxEdge < bestScore) {
                bestScore = maxEdge;
                bestEarIdx = i;
            }
        }

        if (bestEarIdx === -1) break;
        const p = polygon[(bestEarIdx - 1 + polygon.length) % polygon.length];
        const n = polygon[(bestEarIdx + 1) % polygon.length];
        mesh.faces.push(p, polygon[bestEarIdx], n);
        polygon.splice(bestEarIdx, 1);
        addedFaces++;
    }

    if (polygon.length === 3) {
        mesh.faces.push(polygon[0], polygon[1], polygon[2]);
        addedFaces++;
    }
    return addedFaces;
}

/**
 * Density Refinement: Subdivide triangles in the patch whose edges exceed maxLen.
 * This ensures the patch has a similar polygon density to the surrounding mesh.
 */
function refinePatch(mesh: IndexedMesh, patchStartFaceIdx: number, maxLen: number, loop: number[]) {
    // Original hole boundary edges should not be split
    const loopEdges = new Set<string>();
    for (let i = 0; i < loop.length; i++) {
        const a = loop[i];
        const b = loop[(i + 1) % loop.length];
        loopEdges.add(`${a}_${b}`);
        loopEdges.add(`${b}_${a}`);
    }

    let i = patchStartFaceIdx;
    let safety = 0;

    // Process faces dynamically as they are added
    while (i < mesh.faces.length / 3 && safety < 50000) {
        safety++;
        const a = mesh.faces[i * 3];
        const b = mesh.faces[i * 3 + 1];
        const c = mesh.faces[i * 3 + 2];

        const edges = [
            { u: a, v: b, opp: c, d: vertexDist(mesh, a, b) },
            { u: b, v: c, opp: a, d: vertexDist(mesh, b, c) },
            { u: c, v: a, opp: b, d: vertexDist(mesh, c, a) }
        ];

        edges.sort((e1, e2) => e2.d - e1.d);

        let splitEdge = null;
        for (const edge of edges) {
            if (edge.d > maxLen && !loopEdges.has(`${edge.u}_${edge.v}`)) {
                splitEdge = edge;
                break;
            }
        }

        if (splitEdge) {
            const { u, v, opp } = splitEdge;
            const midX = (mesh.vertices[u * 3] + mesh.vertices[v * 3]) / 2;
            const midY = (mesh.vertices[u * 3 + 1] + mesh.vertices[v * 3 + 1]) / 2;
            const midZ = (mesh.vertices[u * 3 + 2] + mesh.vertices[v * 3 + 2]) / 2;
            const midIdx = mesh.vertices.length / 3;
            mesh.vertices.push(midX, midY, midZ);

            // Replace current face i with (u, mid, opp) - preserves winding
            mesh.faces[i * 3] = u;
            mesh.faces[i * 3 + 1] = midIdx;
            mesh.faces[i * 3 + 2] = opp;
            // Add new face (mid, v, opp)
            mesh.faces.push(midIdx, v, opp);

            // Find neighbor face sharing this internal edge (v, u)
            let neighborFaceIdx = -1;
            let n_u = -1, n_v = -1, n_opp = -1;

            // Search from patchStartFaceIdx (internal edges are only shared with other patch faces, 
            // since we don't split boundary loop edges)
            for (let j = patchStartFaceIdx; j < mesh.faces.length / 3; j++) {
                if (j === i) continue;
                const na = mesh.faces[j * 3];
                const nb = mesh.faces[j * 3 + 1];
                const nc = mesh.faces[j * 3 + 2];

                // Directed edge in neighbor should be (v, u) to maintain manifoldness
                if (na === v && nb === u) { n_u = na; n_v = nb; n_opp = nc; neighborFaceIdx = j; break; }
                if (nb === v && nc === u) { n_u = nb; n_v = nc; n_opp = na; neighborFaceIdx = j; break; }
                if (nc === v && na === u) { n_u = nc; n_v = na; n_opp = nb; neighborFaceIdx = j; break; }
                // Also check opposite direction if orientation was flipped somehow
                if (na === u && nb === v) { n_u = na; n_v = nb; n_opp = nc; neighborFaceIdx = j; break; }
                if (nb === u && nc === v) { n_u = nb; n_v = nc; n_opp = na; neighborFaceIdx = j; break; }
                if (nc === u && na === v) { n_u = nc; n_v = na; n_opp = nb; neighborFaceIdx = j; break; }
            }

            if (neighborFaceIdx !== -1) {
                mesh.faces[neighborFaceIdx * 3] = n_u;
                mesh.faces[neighborFaceIdx * 3 + 1] = midIdx;
                mesh.faces[neighborFaceIdx * 3 + 2] = n_opp;
                mesh.faces.push(midIdx, n_v, n_opp);
            }
            // Do not increment i, re-evaluate the two new halves of face i
        } else {
            i++;
        }
    }
}

/**
 * Laplacian Fairing: Smooths the interior vertices of the patch to naturally match the curvature
 * of the surrounding hole boundary by using Bi-Laplacian (curvature) forces.
 */
function fairPatch(mesh: IndexedMesh, patchStartFaceIdx: number, loop: number[], iterations: number = 30) {
    const boundarySet = new Set(loop);

    const patchVertices = new Set<number>();
    for (let i = patchStartFaceIdx; i < mesh.faces.length / 3; i++) {
        patchVertices.add(mesh.faces[i * 3]);
        patchVertices.add(mesh.faces[i * 3 + 1]);
        patchVertices.add(mesh.faces[i * 3 + 2]);
    }

    const interiorVerts: number[] = [];
    patchVertices.forEach(v => {
        if (!boundarySet.has(v)) interiorVerts.push(v);
    });

    if (interiorVerts.length === 0) return;

    // Use a 1-ring expansion around the boundary to anchor the curvature.
    const activeVertices = new Set<number>(patchVertices);
    const adj = new Map<number, Set<number>>();

    for (let i = 0; i < mesh.faces.length / 3; i++) {
        const a = mesh.faces[i * 3], b = mesh.faces[i * 3 + 1], c = mesh.faces[i * 3 + 2];
        if (activeVertices.has(a) || activeVertices.has(b) || activeVertices.has(c)) {
            if (!adj.has(a)) adj.set(a, new Set());
            if (!adj.has(b)) adj.set(b, new Set());
            if (!adj.has(c)) adj.set(c, new Set());
            adj.get(a)!.add(b); adj.get(a)!.add(c);
            adj.get(b)!.add(a); adj.get(b)!.add(c);
            adj.get(c)!.add(a); adj.get(c)!.add(b);
            activeVertices.add(a);
            activeVertices.add(b);
            activeVertices.add(c);
        }
    }

    const laplaciansX = new Float32Array(mesh.vertices.length);
    const laplaciansY = new Float32Array(mesh.vertices.length);
    const laplaciansZ = new Float32Array(mesh.vertices.length);
    const activeArray = Array.from(activeVertices);

    for (let iter = 0; iter < iterations; iter++) {
        // Step 1: Umbrella Laplacian L(v) for all active vertices
        for (let i = 0; i < activeArray.length; i++) {
            const v = activeArray[i];
            const neighbors = adj.get(v)!;
            let sumX = 0, sumY = 0, sumZ = 0;
            let count = 0;
            for (const n of neighbors) {
                sumX += mesh.vertices[n * 3];
                sumY += mesh.vertices[n * 3 + 1];
                sumZ += mesh.vertices[n * 3 + 2];
                count++;
            }
            if (count > 0) {
                laplaciansX[v] = (sumX / count) - mesh.vertices[v * 3];
                laplaciansY[v] = (sumY / count) - mesh.vertices[v * 3 + 1];
                laplaciansZ[v] = (sumZ / count) - mesh.vertices[v * 3 + 2];
            } else {
                laplaciansX[v] = 0; laplaciansY[v] = 0; laplaciansZ[v] = 0;
            }
        }

        // Step 2: Bi-Laplacian explicitly L^2(v) for interior vertices
        for (let i = 0; i < interiorVerts.length; i++) {
            const v = interiorVerts[i];
            const neighbors = adj.get(v)!;
            let sumLX = 0, sumLY = 0, sumLZ = 0;
            let count = 0;
            for (const n of neighbors) {
                sumLX += laplaciansX[n];
                sumLY += laplaciansY[n];
                sumLZ += laplaciansZ[n];
                count++;
            }

            if (count > 0) {
                const biLapX = (sumLX / count) - laplaciansX[v];
                const biLapY = (sumLY / count) - laplaciansY[v];
                const biLapZ = (sumLZ / count) - laplaciansZ[v];

                // Update: v_new = v - k * L^2(v)
                const k = 0.25;
                mesh.vertices[v * 3] -= biLapX * k;
                mesh.vertices[v * 3 + 1] -= biLapY * k;
                mesh.vertices[v * 3 + 2] -= biLapZ * k;
            }
        }
    }
}

/**
 * Main fill dispatcher: choose fan-fill for planar holes (flat caps)
 * or ear-clipping for complex non-planar holes.
 */
function fillHoleWithSmallTriangles(mesh: IndexedMesh, loop: number[]): number {
    const planarityRatio = holePlanarityRatio(mesh, loop);
    // If boundary vertices lie close to a plane (< 15% deviation), use fan fill
    if (planarityRatio < 0.15) {
        return fillHoleWithFan(mesh, loop);
    }
    return fillHoleEarClip(mesh, loop);
}


// ─── Annular Gap (Ring) Fastening ─────────────────────────────────────

/**
 * Zips two boundary loops together using a shortest-diagonal advancing front.
 * Assumes loopA and loopB are roughly parallel and facing each other.
 * Returns the number of triangles added.
 */
function bridgeAnnularGap(mesh: IndexedMesh, loopA: number[], loopB: number[]): number {
    if (loopA.length === 0 || loopB.length === 0) return 0;

    // Find the pair of vertices with the absolute minimum distance between the two loops
    let minD = Infinity;
    let startA = 0;
    let startB = 0;

    for (let i = 0; i < loopA.length; i++) {
        for (let j = 0; j < loopB.length; j++) {
            const d = vertexDist(mesh, loopA[i], loopB[j]);
            if (d < minD) {
                minD = d;
                startA = i;
                startB = j;
            }
        }
    }

    // Since the loops might be oriented in the same or opposite winding relative to the "tube"
    // we are forming, we must determine the correct stepping direction for loopB so the triangles
    // face outward (manifoldness). We test both directions for the first step and pick the
    // one that minimizes the crossed diagonal length.
    const a1 = loopA[(startA + 1) % loopA.length];

    const bPlus = loopB[(startB + 1) % loopB.length];
    const bMinus = loopB[(startB - 1 + loopB.length) % loopB.length];

    const dPlus = vertexDist(mesh, a1, bPlus);
    const dMinus = vertexDist(mesh, a1, bMinus);

    const stepB = dPlus < dMinus ? 1 : -1;

    let ptrA = startA;
    let ptrB = startB;
    let countA = 0;
    let countB = 0;
    let addedFaces = 0;

    // Advance the front until we have traversed both loops entirely
    while (countA < loopA.length || countB < loopB.length) {
        const currA = loopA[ptrA];
        const currB = loopB[ptrB];

        const nextAIdx = (ptrA + 1) % loopA.length;
        const nextBIdx = (ptrB + stepB + loopB.length) % loopB.length;

        const nextA = loopA[nextAIdx];
        const nextB = loopB[nextBIdx];

        // If one loop is exhausted, we must step the other
        if (countA === loopA.length) {
            mesh.faces.push(currA, currB, nextB);
            ptrB = nextBIdx;
            countB++;
            addedFaces++;
        } else if (countB === loopB.length) {
            mesh.faces.push(currA, currB, nextA);
            ptrA = nextAIdx;
            countA++;
            addedFaces++;
        } else {
            // Compare the two possible diagonals to advance the front
            // Diagonal 1: currA to nextB
            const diag1 = vertexDist(mesh, currA, nextB);
            // Diagonal 2: nextA to currB
            const diag2 = vertexDist(mesh, nextA, currB);

            if (diag2 <= diag1) {
                // Advance A: Triangle is (currA, currB, nextA)
                // Winding handles outward face if the loops originated from solid hull
                mesh.faces.push(currA, currB, nextA);
                ptrA = nextAIdx;
                countA++;
            } else {
                // Advance B: Triangle is (currA, currB, nextB)
                mesh.faces.push(currA, currB, nextB);
                ptrB = nextBIdx;
                countB++;
            }
            addedFaces++;
        }
    }

    return addedFaces;
}

/**
 * Sweeps all discovered loops. If two loops are found to be physically adjacent (annular gap),
 * bridges them and removes them from the pending list.
 * Returns the unbridged remaining loops.
 */
function detectAndBridgeAnnularGaps(mesh: IndexedMesh, loops: number[][], searchDiameter: number): { remainingLoops: number[][], totalFacesAdded: number } {
    const remaining = [...loops];
    let totalFacesAdded = 0;

    let i = 0;
    while (i < remaining.length) {
        let bridged = false;
        const loopA = remaining[i];

        // Find best candidate for Loop A
        let bestCandidateIdx = -1;
        let bestDist = searchDiameter;

        for (let j = i + 1; j < remaining.length; j++) {
            const loopB = remaining[j];

            // Fast check: are the bounding boxes or centroids close?
            // To be robust, we just find the closest point pair between A and B
            let minD = Infinity;
            for (let a = 0; a < Math.min(20, loopA.length); a++) {
                const va = loopA[Math.floor((a / 20) * loopA.length)];
                for (let b = 0; b < Math.min(20, loopB.length); b++) {
                    const vb = loopB[Math.floor((b / 20) * loopB.length)];
                    const d = vertexDist(mesh, va, vb);
                    if (d < minD) minD = d;
                }
            }

            if (minD < bestDist) {
                bestDist = minD;
                bestCandidateIdx = j;
            }
        }

        if (bestCandidateIdx !== -1) {
            const loopB = remaining[bestCandidateIdx];
            const added = bridgeAnnularGap(mesh, loopA, loopB);
            console.log(`[MeshRepair] Bridged annular gap between loop of ${loopA.length} edges and ${loopB.length} edges. Added ${added} triangles. Gap dist: ${bestDist.toFixed(4)}`);
            totalFacesAdded += added;

            // Remove both loops from pending list
            remaining.splice(bestCandidateIdx, 1);
            remaining.splice(i, 1);
            bridged = true;
        }

        if (!bridged) {
            i++;
        }
    }

    return { remainingLoops: remaining, totalFacesAdded };
}


// ─── Output conversion ──────────────────────────────────────────────

function indexedMeshToBufferGeometry(mesh: IndexedMesh): THREE.BufferGeometry {
    const numFaces = mesh.faces.length / 3;
    const numVerts = mesh.vertices.length / 3;
    const posArray = new Float32Array(numFaces * 9);
    const v = mesh.vertices;

    for (let fi = 0; fi < numFaces; fi++) {
        for (let vi = 0; vi < 3; vi++) {
            const vertIdx = mesh.faces[fi * 3 + vi];
            if (vertIdx >= 0 && vertIdx < numVerts) {
                posArray[fi * 9 + vi * 3] = v[vertIdx * 3];
                posArray[fi * 9 + vi * 3 + 1] = v[vertIdx * 3 + 1];
                posArray[fi * 9 + vi * 3 + 2] = v[vertIdx * 3 + 2];
            }
        }
    }

    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(posArray, 3));
    geo.computeVertexNormals();
    return geo;
}

// ─── Manifold validation (optional) ──────────────────────────────────

function validateWithManifold(Manifold: any, Mesh: any, mesh: IndexedMesh): boolean {
    try {
        const vertProperties = new Float32Array(mesh.vertices);
        const triVerts = new Uint32Array(mesh.faces);
        const meshObj = new Mesh({ numProp: 3, vertProperties, triVerts });
        meshObj.merge();
        const manifoldObj = new Manifold(meshObj);
        const empty = manifoldObj.isEmpty();
        manifoldObj.delete();
        return !empty;
    } catch {
        return false;
    }
}

export interface RepairOptions {
    /** Multiplier for max fill triangle edge size (relative to avg edge). Default: 3 */
    edgeMultiplier?: number;
    /** Multiplier for max hole diameter to fill (relative to avg edge). Default: 10 */
    holeDiameterMultiplier?: number;
}

export interface HoleAnalysis {
    diameter: number;
    maxEdge: number;
    edgeCount: number;
    positions: Float32Array;
}

export function analyzeRepairableHoles(geometry: THREE.BufferGeometry): {
    avgEdgeLength: number;
    holes: HoleAnalysis[];
} {
    const positions = geometry.attributes.position.array;
    const mesh = buildIndexedMesh(positions);
    removeDegenerateFaces(mesh);

    const { adj } = findBoundaryEdges(mesh);
    const loops = traceBoundaryLoops(adj);
    const avgEdgeLength = computeAvgEdgeLength(mesh);

    const holes: HoleAnalysis[] = loops.map(loop => {
        const diameter = loopDiameter(mesh, loop);
        const maxEdge = loopMaxEdgeLength(mesh, loop);

        const linePos = new Float32Array(loop.length * 6);
        for (let i = 0; i < loop.length; i++) {
            const curIdx = loop[i];
            const nextIdx = loop[(i + 1) % loop.length];
            linePos[i * 6 + 0] = mesh.vertices[curIdx * 3 + 0];
            linePos[i * 6 + 1] = mesh.vertices[curIdx * 3 + 1];
            linePos[i * 6 + 2] = mesh.vertices[curIdx * 3 + 2];
            linePos[i * 6 + 3] = mesh.vertices[nextIdx * 3 + 0];
            linePos[i * 6 + 4] = mesh.vertices[nextIdx * 3 + 1];
            linePos[i * 6 + 5] = mesh.vertices[nextIdx * 3 + 2];
        }

        return {
            diameter,
            maxEdge,
            edgeCount: loop.length,
            positions: linePos
        };
    });

    return { avgEdgeLength, holes };
}

export async function repairWithManifold(
    geometry: THREE.BufferGeometry,
    options: RepairOptions = {}
): Promise<THREE.BufferGeometry> {
    const {
        edgeMultiplier = 3,
        holeDiameterMultiplier = 10,
    } = options;

    const positions = geometry.attributes.position.array;
    const inputTriCount = positions.length / 9;

    console.log(`[MeshRepair] === Starting conservative repair: ${inputTriCount} triangles ===`);
    console.log(`[MeshRepair] Options: edgeMultiplier=${edgeMultiplier}, holeDiameterMultiplier=${holeDiameterMultiplier}`);

    // Step 1: Build indexed mesh (weld coincident vertices)
    const mesh = buildIndexedMesh(positions);
    const origFaceCount = mesh.faces.length / 3;

    // Step 2: Remove only truly degenerate faces (same vertex indices)
    const degRemoved = removeDegenerateFaces(mesh);
    if (degRemoved > 0) {
        console.log(`[MeshRepair] Removed ${degRemoved} degenerate faces`);
    }

    // Step 3: Find boundary edges (gaps in the mesh)
    const { boundaryCount, adj } = findBoundaryEdges(mesh);
    console.log(`[MeshRepair] Found ${boundaryCount} boundary edges`);

    if (boundaryCount === 0) {
        console.log('[MeshRepair] Mesh is already watertight! No repair needed.');
        return indexedMeshToBufferGeometry(mesh);
    }

    // Step 4: Trace boundary loops (each loop = one hole)
    const loops = traceBoundaryLoops(adj);
    console.log(`[MeshRepair] Traced ${loops.length} boundary loops (holes)`);

    // Step 5: Compute size threshold — configurable via options
    const avgEdge = computeAvgEdgeLength(mesh);
    const maxTriSize = avgEdge * edgeMultiplier;
    const maxHoleDiameter = avgEdge * holeDiameterMultiplier;

    console.log(`[MeshRepair] Avg edge length: ${avgEdge.toFixed(4)}, max fill tri: ${maxTriSize.toFixed(4)}, max hole diameter: ${maxHoleDiameter.toFixed(4)}`);

    // Step 6: Process Annular Gaps (Rings) First
    const { remainingLoops, totalFacesAdded: bridgedFaces } = detectAndBridgeAnnularGaps(mesh, loops, maxHoleDiameter);
    let totalFilled = bridgedFaces;
    let holesSkipped = 0;

    // Step 7: Fill remaining small standalone holes
    for (const loop of remainingLoops) {
        const diameter = loopDiameter(mesh, loop);
        const maxEdge = loopMaxEdgeLength(mesh, loop);

        if (diameter > maxHoleDiameter) {
            console.log(`[MeshRepair] Skipping large hole: ${loop.length} edges, diameter=${diameter.toFixed(4)} (max=${maxHoleDiameter.toFixed(4)})`);
            holesSkipped++;
            continue;
        }

        const patchStartFaceIdx = mesh.faces.length / 3;
        const added = fillHoleWithSmallTriangles(mesh, loop);

        if (added > 0) {
            refinePatch(mesh, patchStartFaceIdx, maxTriSize, loop);
            fairPatch(mesh, patchStartFaceIdx, loop, 50); // 50 Laplacian iterations

            const totalFacesAdded = (mesh.faces.length / 3) - patchStartFaceIdx;
            console.log(`[MeshRepair] Filled & faired hole: ${loop.length} edges, ${totalFacesAdded} triangles added (diameter=${diameter.toFixed(4)}, maxEdge=${maxEdge.toFixed(4)})`);
            totalFilled += totalFacesAdded;
        } else {
            console.log(`[MeshRepair] Could not fill hole: ${loop.length} edges (triangles too large)`);
            holesSkipped++;
        }
    }

    // Step 7: Check remaining boundaries
    const after = findBoundaryEdges(mesh);
    const finalFaces = mesh.faces.length / 3;

    console.log(
        `[MeshRepair] Results: ` +
        `${origFaceCount} → ${finalFaces} faces (+${finalFaces - origFaceCount}), ` +
        `${totalFilled} fill triangles, ${holesSkipped} holes skipped, ` +
        `${after.boundaryCount} boundary edges remaining`
    );

    // Step 8: Try Manifold validation (informational only)
    try {
        const { Manifold, Mesh } = await getWasm();
        const isManifold = validateWithManifold(Manifold, Mesh, mesh);
        console.log(`[MeshRepair] Manifold validation: ${isManifold ? '✓ PASSED' : '✗ not manifold (some gaps remain)'}`);
    } catch (e) {
        console.log('[MeshRepair] Manifold validation skipped (WASM unavailable)');
    }

    // Step 9: Convert back to BufferGeometry
    const result = indexedMeshToBufferGeometry(mesh);
    console.log(`[MeshRepair] ✓ Repair complete`);
    return result;
}
