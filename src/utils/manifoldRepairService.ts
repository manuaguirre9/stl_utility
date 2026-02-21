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
 * Build indexed mesh from a BufferGeometry, respecting its index if it exists.
 * Welds coincident vertices into a unified topological structure.
 */
function buildIndexedMesh(geometry: THREE.BufferGeometry): IndexedMesh {
    const posAttr = geometry.attributes.position;
    const index = geometry.index;
    const PREC = 1e3;
    const vertexMap = new Map<string, number>();
    const vertices: number[] = [];
    const faces: number[] = [];
    let vertCount = 0;

    const getVertIdx = (vIdx: number) => {
        const x = posAttr.getX(vIdx);
        const y = posAttr.getY(vIdx);
        const z = posAttr.getZ(vIdx);

        // Coalesce NaNs/Infs to zero to prevent downstream geometric failures
        const safeX = Number.isFinite(x) ? x : 0;
        const safeY = Number.isFinite(y) ? y : 0;
        const safeZ = Number.isFinite(z) ? z : 0;

        const key = `${Math.round(safeX * PREC)},${Math.round(safeY * PREC)},${Math.round(safeZ * PREC)}`;
        let idx = vertexMap.get(key);
        if (idx === undefined) {
            idx = vertCount++;
            vertexMap.set(key, idx);
            vertices.push(safeX, safeY, safeZ);
        }
        return idx;
    };

    if (index) {
        for (let i = 0; i < index.count; i++) {
            faces.push(getVertIdx(index.getX(i)));
        }
    } else {
        for (let i = 0; i < posAttr.count; i++) {
            faces.push(getVertIdx(i));
        }
    }

    return { vertices, faces };
}

/**
 * Resolves T-Junctions purely topologically by splitting large triangles that have 
 * multiple smaller boundary vertices lying collinearly along their open edge.
 * This natively stitches non-conformal mesh boundaries like knurl edges against flat caps.
 */
function resolveTJunctions(mesh: IndexedMesh): number {
    let splits = 0;
    let changed = true;
    let safety = 0;

    while (changed && safety < 100) {
        changed = false;
        safety++;

        const edgeFaceCount = new Map<string, number>();
        for (let fi = 0; fi < mesh.faces.length / 3; fi++) {
            const v0 = mesh.faces[fi * 3], v1 = mesh.faces[fi * 3 + 1], v2 = mesh.faces[fi * 3 + 2];
            for (const [a, b] of [[v0, v1], [v1, v2], [v2, v0]]) {
                const key = edgeKey(a, b);
                edgeFaceCount.set(key, (edgeFaceCount.get(key) || 0) + 1);
            }
        }

        const boundaryEdges: { u: number, v: number, faceIdx: number, opp: number }[] = [];
        for (let fi = 0; fi < mesh.faces.length / 3; fi++) {
            const v0 = mesh.faces[fi * 3], v1 = mesh.faces[fi * 3 + 1], v2 = mesh.faces[fi * 3 + 2];
            const edges = [[v0, v1, v2], [v1, v2, v0], [v2, v0, v1]];
            for (const [u, v, opp] of edges) {
                if (edgeFaceCount.get(edgeKey(u, v)) === 1) {
                    boundaryEdges.push({ u, v, faceIdx: fi, opp });
                }
            }
        }

        const bVerts = new Set<number>();
        for (const e of boundaryEdges) {
            bVerts.add(e.u);
            bVerts.add(e.v);
        }
        const bVertsArr = Array.from(bVerts);

        const grid = new Map<string, number[]>();
        const gridSize = 2.0;
        for (const p of bVertsArr) {
            const gx = Math.floor(mesh.vertices[p * 3] / gridSize);
            const gy = Math.floor(mesh.vertices[p * 3 + 1] / gridSize);
            const gz = Math.floor(mesh.vertices[p * 3 + 2] / gridSize);
            const key = `${gx},${gy},${gz}`;
            if (!grid.has(key)) grid.set(key, []);
            grid.get(key)!.push(p);
        }

        const modifiedFaces = new Set<number>();

        for (const e of boundaryEdges) {
            if (modifiedFaces.has(e.faceIdx)) continue;

            const ux = mesh.vertices[e.u * 3], uy = mesh.vertices[e.u * 3 + 1], uz = mesh.vertices[e.u * 3 + 2];
            const vx = mesh.vertices[e.v * 3], vy = mesh.vertices[e.v * 3 + 1], vz = mesh.vertices[e.v * 3 + 2];

            const minX = Math.min(ux, vx) - 0.01, maxX = Math.max(ux, vx) + 0.01;
            const minY = Math.min(uy, vy) - 0.01, maxY = Math.max(uy, vy) + 0.01;
            const minZ = Math.min(uz, vz) - 0.01, maxZ = Math.max(uz, vz) + 0.01;

            const uvx = vx - ux, uvy = vy - uy, uvz = vz - uz;
            const edgeLenSq = uvx * uvx + uvy * uvy + uvz * uvz;
            if (edgeLenSq < 1e-8) continue;

            // Find ALL collinear points for this edge
            const collinearVerts: { p: number, t: number }[] = [];

            // Spatial query
            const gsx = Math.floor(minX / gridSize), gex = Math.floor(maxX / gridSize);
            const gsy = Math.floor(minY / gridSize), gey = Math.floor(maxY / gridSize);
            const gsz = Math.floor(minZ / gridSize), gez = Math.floor(maxZ / gridSize);

            for (let gx = gsx; gx <= gex; gx++) {
                for (let gy = gsy; gy <= gey; gy++) {
                    for (let gz = gsz; gz <= gez; gz++) {
                        const candidates = grid.get(`${gx},${gy},${gz}`);
                        if (!candidates) continue;

                        for (const p of candidates) {
                            if (p === e.u || p === e.v || p === e.opp) continue;

                            const px = mesh.vertices[p * 3], py = mesh.vertices[p * 3 + 1], pz = mesh.vertices[p * 3 + 2];
                            if (px < minX || px > maxX || py < minY || py > maxY || pz < minZ || pz > maxZ) continue;

                            const upx = px - ux, upy = py - uy, upz = pz - uz;
                            const t = (upx * uvx + upy * uvy + upz * uvz) / edgeLenSq;

                            if (t > 1e-4 && t < 1 - 1e-4) {
                                const projX = ux + t * uvx, projY = uy + t * uvy, projZ = uz + t * uvz;
                                const distSq = (px - projX) ** 2 + (py - projY) ** 2 + (pz - projZ) ** 2;
                                if (distSq < 1e-8) {
                                    collinearVerts.push({ p, t });
                                }
                            }
                        }
                    }
                }
            }

            if (collinearVerts.length > 0) {
                // Sort by distance from u to v
                collinearVerts.sort((a, b) => a.t - b.t);

                // Replace the original face with the first slice of the fan
                mesh.faces[e.faceIdx * 3 + 0] = e.u;
                mesh.faces[e.faceIdx * 3 + 1] = collinearVerts[0].p;
                mesh.faces[e.faceIdx * 3 + 2] = e.opp;
                modifiedFaces.add(e.faceIdx);
                splits++;

                // Add intermediate slices
                for (let i = 0; i < collinearVerts.length - 1; i++) {
                    mesh.faces.push(collinearVerts[i].p, collinearVerts[i + 1].p, e.opp);
                    splits++;
                }

                // Add the final slice
                mesh.faces.push(collinearVerts[collinearVerts.length - 1].p, e.v, e.opp);

                changed = true;
            }
        }
    }
    return splits;
}

/**
 * Iteratively removes "flaps": faces that have at least one boundary edge (count=1)
 * AND at least one non-manifold edge (count > 2).
 * These are typically zero-volume walls or spikes extending off the main mesh.
 */
function peelLooseFlaps(mesh: IndexedMesh): number {
    let totalRemoved = 0;
    let changed = true;
    let passCount = 0;
    const MAX_PASSES = 20; // Reduced from 100 to prevent hang on complex failures

    while (changed && passCount < MAX_PASSES) {
        changed = false;
        passCount++;
        const edgeCount = new Map<string, number>();
        for (let i = 0; i < mesh.faces.length; i += 3) {
            for (const [a, b] of [[mesh.faces[i], mesh.faces[i + 1]], [mesh.faces[i + 1], mesh.faces[i + 2]], [mesh.faces[i + 2], mesh.faces[i]]]) {
                const k = edgeKey(a, b);
                edgeCount.set(k, (edgeCount.get(k) || 0) + 1);
            }
        }

        const clean: number[] = [];
        let removedThisPass = 0;
        for (let i = 0; i < mesh.faces.length; i += 3) {
            const keys = [
                edgeKey(mesh.faces[i], mesh.faces[i + 1]),
                edgeKey(mesh.faces[i + 1], mesh.faces[i + 2]),
                edgeKey(mesh.faces[i + 2], mesh.faces[i])
            ];

            const counts = keys.map(k => edgeCount.get(k)!);
            const hasBoundary = counts.some(c => c === 1);
            const hasNonManifold = counts.some(c => c > 2);
            const isIsolated = counts.every(c => c === 1);

            if ((hasBoundary && hasNonManifold) || isIsolated) {
                removedThisPass++;
            } else {
                clean.push(mesh.faces[i], mesh.faces[i + 1], mesh.faces[i + 2]);
            }
        }

        if (removedThisPass > 0) {
            mesh.faces = clean;
            totalRemoved += removedThisPass;
            changed = true;
        }
    }
    return totalRemoved;
}

/**
 * Trims un-closeable open boundaries (orange edges that hit dead ends).
 * If a boundary path cannot form a closed loop, the triangles touching it are eroded
 * iteratively until the crack merges into a valid manifold hole or disappears.
 */
function removeUnloopableBoundaries(mesh: IndexedMesh): number {
    let totalRemoved = 0;
    let changed = true;
    let passCount = 0;
    const MAX_PASSES = 20; // Hard fail-safe to prevent infinite topological oscillation

    while (changed && passCount < MAX_PASSES) {
        changed = false;
        passCount++;
        const { boundaryCount, adj } = findBoundaryEdges(mesh);
        if (boundaryCount === 0) break;

        const visitedEdges = new Set<string>();
        const badEdges = new Set<string>();

        for (const startNode of adj.keys()) {
            const neighbors = adj.get(startNode);
            if (!neighbors) continue;

            for (const nextNode of neighbors) {
                const startEdgeKey = `${startNode}_${nextNode}`;
                if (visitedEdges.has(startEdgeKey)) continue;

                const path: string[] = [];
                let curr = startNode;
                let next = nextNode;
                let stuck = false;

                while (true) {
                    const edgeKey = `${curr}_${next}`;
                    path.push(edgeKey);
                    visitedEdges.add(edgeKey);

                    if (next === startNode) break; // closed loop
                    if (path.length > 50000) { stuck = true; break; }

                    const nexts = adj.get(next);
                    if (!nexts) { stuck = true; break; } // dead end

                    let foundNext = -1;
                    for (const candidate of nexts) {
                        if (!visitedEdges.has(`${next}_${candidate}`)) {
                            foundNext = candidate;
                            break;
                        }
                    }
                    if (foundNext === -1) { stuck = true; break; }

                    curr = next;
                    next = foundNext;
                }

                if (stuck || path.length < 3) {
                    for (const ek of path) badEdges.add(ek);
                }
            }
        }

        if (badEdges.size > 0) {
            const clean: number[] = [];
            let removedPass = 0;
            for (let i = 0; i < mesh.faces.length; i += 3) {
                const a = mesh.faces[i], b = mesh.faces[i + 1], c = mesh.faces[i + 2];
                if (
                    badEdges.has(`${a}_${b}`) || badEdges.has(`${b}_${a}`) ||
                    badEdges.has(`${b}_${c}`) || badEdges.has(`${c}_${b}`) ||
                    badEdges.has(`${c}_${a}`) || badEdges.has(`${a}_${c}`)
                ) {
                    removedPass++;
                    changed = true;
                } else {
                    clean.push(a, b, c);
                }
            }
            mesh.faces = clean;
            totalRemoved += removedPass;
        }
    }
    return totalRemoved;
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
function traceBoundaryLoops(adj: Map<number, number[]>, maxEdges: number = 50000): number[][] {
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
    rawBoundaries: Float32Array;
} {
    const mesh = buildIndexedMesh(geometry);
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

    // Collect all raw edges that are boundaries
    const rawSegments: number[] = [];
    for (const [u, list] of adj.entries()) {
        for (const v of list) {
            rawSegments.push(
                mesh.vertices[u * 3 + 0], mesh.vertices[u * 3 + 1], mesh.vertices[u * 3 + 2],
                mesh.vertices[v * 3 + 0], mesh.vertices[v * 3 + 1], mesh.vertices[v * 3 + 2]
            );
        }
    }

    return { avgEdgeLength, holes, rawBoundaries: new Float32Array(rawSegments) };
}

export async function repairWithManifold(
    geometry: THREE.BufferGeometry,
    options: RepairOptions = {}
): Promise<THREE.BufferGeometry> {
    const {
        edgeMultiplier = 3,
        holeDiameterMultiplier = 10,
    } = options;

    const inputTriCount = geometry.index ? geometry.index.count / 3 : geometry.attributes.position.count / 3;

    console.log(`[MeshRepair] === Starting conservative repair: ${inputTriCount} triangles ===`);
    console.log(`[MeshRepair] Options: edgeMultiplier=${edgeMultiplier}, holeDiameterMultiplier=${holeDiameterMultiplier}`);

    // Step 1: Build indexed mesh (weld coincident vertices)
    const mesh = buildIndexedMesh(geometry);
    const origFaceCount = mesh.faces.length / 3;

    // Step 2: Remove only truly degenerate faces (same vertex indices)
    const degRemoved = removeDegenerateFaces(mesh);
    if (degRemoved > 0) {
        console.log(`[MeshRepair] Removed ${degRemoved} degenerate faces`);
    }

    // Step 2.5: Topologically resolve T-junctions
    const tSplits = resolveTJunctions(mesh);
    if (tSplits > 0) {
        console.log(`[MeshRepair] Resolved ${tSplits} T-junction instances (split boundary faces).`);
    }

    // Step 2.6: Peel non-manifold flaps and unloopable boundaries
    const flapsRemoved = peelLooseFlaps(mesh);
    if (flapsRemoved > 0) console.log(`[MeshRepair] Peeled ${flapsRemoved} non-manifold loose flaps.`);

    const unloopableRemoved = removeUnloopableBoundaries(mesh);
    if (unloopableRemoved > 0) console.log(`[MeshRepair] Eroded ${unloopableRemoved} triangles along un-closeable open cracks.`);

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
    const annularGapThreshold = avgEdge * 3;
    const { remainingLoops, totalFacesAdded: bridgedFaces } = detectAndBridgeAnnularGaps(mesh, loops, annularGapThreshold);
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
