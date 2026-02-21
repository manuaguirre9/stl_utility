import * as THREE from 'three';
import { getContiguousIslands, subdivideTriangles, subdivideSelectedFacesToSize, fitCylinderToSelection } from './meshUtils';
import { repairWithManifold } from './manifoldRepairService';

export type KnurlPattern = 'diamond' | 'straight' | 'diagonal' | 'square';

export interface KnurlingParams {
    type: 'knurling';
    pitch: number;
    depth: number;
    angle: number;
    pattern: KnurlPattern;
    holeFillThreshold: number;
    holeFillEnabled: boolean;
}

export interface HoneycombParams {
    type: 'honeycomb';
    cellSize: number;
    wallThickness: number;
    depth: number;
    angle: number;
    direction: 'inward' | 'outward';
    holeFillThreshold: number;
    holeFillEnabled: boolean;
}

export interface FuzzySkinParams {
    type: 'fuzzy';
    thickness: number;
    pointDistance: number;
    holeFillThreshold: number;
    holeFillEnabled: boolean;
}

// --------------------------------------------------------------------------------
// SHARED UTILITIES
// --------------------------------------------------------------------------------

function weldGeometry(geo: THREE.BufferGeometry, prec: number = 1000): THREE.BufferGeometry {
    const pos = geo.attributes.position, vMap = new Map<string, number>(), nPos: number[] = [], indices: number[] = [];
    for (let i = 0; i < pos.count; i++) {
        const x = pos.getX(i), y = pos.getY(i), z = pos.getZ(i), h = `${Math.round(x * prec)},${Math.round(y * prec)},${Math.round(z * prec)}`;
        if (!vMap.has(h)) { vMap.set(h, nPos.length / 3); nPos.push(x, y, z); }
        indices.push(vMap.get(h)!);
    }
    const fIdx: number[] = [];
    for (let f = 0; f < indices.length / 3; f++) {
        const i1 = indices[f * 3], i2 = indices[f * 3 + 1], i3 = indices[f * 3 + 2];
        if (i1 !== i2 && i2 !== i3 && i3 !== i1) fIdx.push(i1, i2, i3);
    }
    const res = new THREE.BufferGeometry();
    res.setAttribute('position', new THREE.BufferAttribute(new Float32Array(nPos), 3));
    res.setIndex(fIdx);
    return res.toNonIndexed();
}

function getBarycentric(p: { u: number, v: number }, a: { u: number, v: number }, b: { u: number, v: number }, c: { u: number, v: number }) {
    const v0u = b.u - a.u, v0v = b.v - a.v, v1u = c.u - a.u, v1v = c.v - a.v, v2u = p.u - a.u, v2v = p.v - a.v;
    const den = v0u * v1v - v1u * v0v;
    if (Math.abs(den) < 1e-12) return { wA: 1, wB: 0, wC: 0 };
    const v = (v2u * v1v - v1u * v2v) / den;
    const w = (v0u * v2v - v2u * v0v) / den;
    return { wA: 1 - v - w, wB: v, wC: w };
}

function distSq(x1: number, y1: number, z1: number, x2: number, y2: number, z2: number): number {
    return (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2;
}

// --------------------------------------------------------------------------------
// BOUNDARY UTILS
// --------------------------------------------------------------------------------



/**
 * Returns an array of boundary line segments (start, end) for geometric checking.
 * Allows fixing gaps by snapping/stitching to ANY point on a boundary edge.
 */
function getBoundarySegments(geometry: THREE.BufferGeometry, faceIndices: number[]): { a: THREE.Vector3, b: THREE.Vector3, lenSq: number, ab: THREE.Vector3 }[] {
    const posAttr = geometry.attributes.position;
    const PREC = 1e3;
    const vKey = (idx: number) => `${Math.round(posAttr.getX(idx) * PREC)},${Math.round(posAttr.getY(idx) * PREC)},${Math.round(posAttr.getZ(idx) * PREC)}`;

    const edgeMap = new Map<string, { a: THREE.Vector3, b: THREE.Vector3, count: number }>();

    const index = geometry.index;

    faceIndices.forEach(fIdx => {
        const i0 = index ? index.getX(fIdx * 3 + 0) : fIdx * 3 + 0;
        const i1 = index ? index.getX(fIdx * 3 + 1) : fIdx * 3 + 1;
        const i2 = index ? index.getX(fIdx * 3 + 2) : fIdx * 3 + 2;

        const v0 = new THREE.Vector3(posAttr.getX(i0), posAttr.getY(i0), posAttr.getZ(i0));
        const v1 = new THREE.Vector3(posAttr.getX(i1), posAttr.getY(i1), posAttr.getZ(i1));
        const v2 = new THREE.Vector3(posAttr.getX(i2), posAttr.getY(i2), posAttr.getZ(i2));

        const k0 = vKey(i0), k1 = vKey(i1), k2 = vKey(i2);

        const processEdge = (keyA: string, keyB: string, vecA: THREE.Vector3, vecB: THREE.Vector3) => {
            const key = keyA < keyB ? `${keyA}:${keyB}` : `${keyB}:${keyA}`;
            if (!edgeMap.has(key)) {
                edgeMap.set(key, { a: vecA, b: vecB, count: 0 });
            }
            edgeMap.get(key)!.count++;
        };

        processEdge(k0, k1, v0, v1);
        processEdge(k1, k2, v1, v2);
        processEdge(k2, k0, v2, v0);
    });

    const segments: { a: THREE.Vector3, b: THREE.Vector3, lenSq: number, ab: THREE.Vector3 }[] = [];
    for (const data of edgeMap.values()) {
        if (data.count === 1) {
            const ab = new THREE.Vector3().subVectors(data.b, data.a);
            segments.push({ a: data.a, b: data.b, lenSq: ab.lengthSq(), ab });
        }
    }
    return segments;
}

// --------------------------------------------------------------------------------
// SPATIAL HASH GRID FOR PERFORMANCE OPTIMIZATION
// --------------------------------------------------------------------------------
// Problem: Checking if a point is on a boundary (isOnAnySegment) scales as O(N) where N is the number of boundary segments.
// When generating textures, we test millions of points against thousands of boundary segments. This O(P * N) complexity
// would freeze the browser for minutes.
// Solution: We partition the 3D space into a grid of cells (Spatial Hashing).
// Each cell stores only the boundary segments that pass through it.
// Finding if a point is on a boundary becomes an O(1) operation: we just look up the cell the point is in and check
// only the few segments (if any) inside that specific cell.

interface BoundaryGrid {
    cellSize: number;
    cells: Map<string, { a: THREE.Vector3, b: THREE.Vector3, lenSq: number, ab: THREE.Vector3 }[]>;
}

/**
 * Creates a Spatial Hash Grid from a list of boundary segments.
 * @param segments The list of boundary segments extracted from the original mesh selection.
 * @param cellSize The size of each 3D cubic cell. 2.0mm is a good default for 3D printing tasks.
 */
function createBoundaryGrid(segments: { a: THREE.Vector3, b: THREE.Vector3, lenSq: number, ab: THREE.Vector3 }[], cellSize: number = 2.0): BoundaryGrid {
    const grid: BoundaryGrid = { cellSize, cells: new Map() };
    for (const seg of segments) {
        // Find the bounding box of the segment, expanded by a tiny buffer
        const minX = Math.min(seg.a.x, seg.b.x) - 0.1;
        const maxX = Math.max(seg.a.x, seg.b.x) + 0.1;
        const minY = Math.min(seg.a.y, seg.b.y) - 0.1;
        const maxY = Math.max(seg.a.y, seg.b.y) + 0.1;
        const minZ = Math.min(seg.a.z, seg.b.z) - 0.1;
        const maxZ = Math.max(seg.a.z, seg.b.z) + 0.1;

        // Convert the continuous 3D bounding box into integer grid coordinates
        const startX = Math.floor(minX / cellSize);
        const endX = Math.floor(maxX / cellSize);
        const startY = Math.floor(minY / cellSize);
        const endY = Math.floor(maxY / cellSize);
        const startZ = Math.floor(minZ / cellSize);
        const endZ = Math.floor(maxZ / cellSize);

        // Register the segment in every cell it touches
        for (let x = startX; x <= endX; x++) {
            for (let y = startY; y <= endY; y++) {
                for (let z = startZ; z <= endZ; z++) {
                    const key = `${x},${y},${z}`;
                    if (!grid.cells.has(key)) grid.cells.set(key, []);
                    grid.cells.get(key)!.push(seg);
                }
            }
        }
    }
    return grid;
}

const _ap = new THREE.Vector3(); // Reused vector to avoid garbage collection pressure in tight loops

/**
 * Checks if a given 3D point `p` lies on ANY of the segments stored in the Spatial Hash Grid.
 * @param p The 3D point to test (usually a midpoint of an outer polygon edge).
 * @param grid The pre-computed Spatial Hash Grid.
 * @param eps Tolerance (in mm) for floating point inaccuracies. Default 1e-2. 1e-4 used for stricter exact matching.
 */
function isOnAnySegmentGrid(p: THREE.Vector3, grid: BoundaryGrid, eps: number = 1e-2): boolean {
    const { cellSize, cells } = grid;

    // 1. Locate the exact grid cell the point falls into
    const x = Math.floor(p.x / cellSize);
    const y = Math.floor(p.y / cellSize);
    const z = Math.floor(p.z / cellSize);

    // 2. Fetch only the segments registered in that cell
    const segments = cells.get(`${x},${y},${z}`);

    // If no segments exist in this region, the point is definitely not on the boundary (O(1) early exit)
    if (!segments) return false;

    // 3. Perform the exact point-to-line-segment distance calculation for the few candidates
    const epsSq = eps * eps;
    for (const seg of segments) {
        if (seg.lenSq < 1e-12) {
            // Handle edge-case: degenerate segments (points)
            if (p.distanceToSquared(seg.a) < epsSq) return true;
            continue;
        }

        // Vector from segment start A to point P
        _ap.subVectors(p, seg.a);

        // Project AP onto AB to find normalized position t along the segment
        const t = _ap.dot(seg.ab) / seg.lenSq;

        // Is the projection strictly between pointA (t=0) and pointB (t=1)?
        if (t >= -1e-4 && t <= 1 + 1e-4) {
            // Calculate squared distance from point to the line segment
            // d^2 = |AP|^2 - (|AP| cos(theta))^2
            const dSq = Math.max(0, _ap.lengthSq() - (_ap.dot(seg.ab) ** 2) / seg.lenSq);
            if (dSq < epsSq) return true;
        }
    }
    return false;
}

// --------------------------------------------------------------------------------
// GENERIC TEXTURE ENGINE
// --------------------------------------------------------------------------------

interface ProjectionData {
    avgR: number;
    circPhys: number;
    depthSign: number;
    isPlanar: boolean;
    project: (p: THREE.Vector3) => { u: number; v: number };
    unproject: (u: number, v: number, h: number) => THREE.Vector3;
    toRot: (u: number, v: number) => { ur: number; vr: number };
    fromRot: (ur: number, vr: number) => { u: number; v: number };
    getNormal: (u: number, v: number) => THREE.Vector3;
}

function createProjectionData(geometry: THREE.BufferGeometry, faceIndices: number[], angleDeg: number): ProjectionData | null {
    const fit = fitCylinderToSelection(geometry, faceIndices);
    if (!fit) return null;

    const angleRad = (angleDeg * Math.PI) / 180;
    const cosA = Math.cos(angleRad), sinA = Math.sin(angleRad);
    const toRot = (u: number, v: number) => ({ ur: u * cosA - v * sinA, vr: u * sinA + v * cosA });
    const fromRot = (ur: number, vr: number) => ({ u: ur * cosA + vr * sinA, v: -ur * sinA + vr * cosA });

    const avgR = fit.avgR!;
    const circPhys = 2 * Math.PI * avgR;
    const depthSign = fit.isFlipped ? -1 : 1;
    const isPlanar = fit.isPlanar;

    const upDir = isPlanar ? fit.avgN! : fit.up!;
    const origin = isPlanar ? fit.origin! : fit.axisOrigin!;
    const rightDir = new THREE.Vector3();
    if (Math.abs(upDir.y) < 0.9) rightDir.set(0, 1, 0).cross(upDir).normalize(); else rightDir.set(1, 0, 0).cross(upDir).normalize();
    const fwdDir = new THREE.Vector3().crossVectors(upDir, rightDir).normalize();

    const project = (p: THREE.Vector3) => {
        const v = p.clone().sub(origin);
        if (isPlanar) return { u: v.dot(rightDir), v: v.dot(fwdDir) };
        let uRad = Math.atan2(v.dot(fwdDir), v.dot(rightDir));
        if (uRad < 0) uRad += 2 * Math.PI;
        return { u: uRad * avgR, v: v.dot(upDir) };
    };

    const unproject = (u: number, v: number, h: number) => {
        let r = 1;
        if (!isPlanar) {
            r = (fit.m! * v + fit.b!) + h * depthSign;
        }
        if (isPlanar) {
            return rightDir.clone().multiplyScalar(u).add(fwdDir.clone().multiplyScalar(v)).add(upDir.clone().multiplyScalar(h)).add(origin);
        }
        const ang = u / avgR;
        const toS = rightDir.clone().multiplyScalar(Math.cos(ang)).add(fwdDir.clone().multiplyScalar(Math.sin(ang)));
        return toS.multiplyScalar(r).add(upDir.clone().multiplyScalar(v)).add(origin);
    };

    const getNormal = (u: number, _v: number) => {
        if (isPlanar) return upDir.clone();
        const ang = u / avgR;
        const radial = rightDir.clone().multiplyScalar(Math.cos(ang)).add(fwdDir.clone().multiplyScalar(Math.sin(ang)));
        const m = fit.m || 0;
        const norm = radial.clone().addScaledVector(upDir, -m).normalize();
        return norm.multiplyScalar(depthSign);
    };

    return { avgR, circPhys, depthSign, isPlanar, project, unproject, toRot, fromRot, getNormal };
}

// --------------------------------------------------------------------------------
// KNURLING IMPLEMENTATION
// --------------------------------------------------------------------------------

export async function applyKnurling(
    geometry: THREE.BufferGeometry,
    faceIndices: number[],
    params: KnurlingParams
): Promise<THREE.BufferGeometry> {
    if (!faceIndices || faceIndices.length === 0) return geometry;

    const index = geometry.index;
    const posAttr = geometry.attributes.position;
    const islands = getContiguousIslands(geometry, faceIndices);

    // --- TOPOLOGICAL INDEXING SETUP ---
    const vertexMap = new Map<string, number>();
    const vertices: number[] = [];
    const newIndices: number[] = [];

    const getVertIndex = (x: number, y: number, z: number): number => {
        const key = `${Math.round(x * 1e5)},${Math.round(y * 1e5)},${Math.round(z * 1e5)}`;
        let idx = vertexMap.get(key);
        if (idx === undefined) {
            idx = vertices.length / 3;
            vertexMap.set(key, idx);
            vertices.push(x, y, z);
        }
        return idx;
    };

    // Pre-compute boundary segments for wall generation
    const boundarySegments = getBoundarySegments(geometry, faceIndices);
    const boundaryGrid = createBoundaryGrid(boundarySegments, 2.0);
    console.log(`[Knurling] Found ${boundarySegments.length} boundary segments.`);
    let wallCount = 0;

    islands.forEach(islandIndices => {
        const islandPos = new Float32Array(islandIndices.length * 9);
        islandIndices.forEach((fIdx, i) => {
            const i0 = index ? index.getX(fIdx * 3 + 0) : fIdx * 3 + 0;
            const i1 = index ? index.getX(fIdx * 3 + 1) : fIdx * 3 + 1;
            const i2 = index ? index.getX(fIdx * 3 + 2) : fIdx * 3 + 2;
            islandPos[i * 9 + 0] = posAttr.getX(i0); islandPos[i * 9 + 1] = posAttr.getY(i0); islandPos[i * 9 + 2] = posAttr.getZ(i0);
            islandPos[i * 9 + 3] = posAttr.getX(i1); islandPos[i * 9 + 4] = posAttr.getY(i1); islandPos[i * 9 + 5] = posAttr.getZ(i1);
            islandPos[i * 9 + 6] = posAttr.getX(i2); islandPos[i * 9 + 7] = posAttr.getY(i2); islandPos[i * 9 + 8] = posAttr.getZ(i2);
        });

        // STEP 1: Step-based Subdivision (to prevent exponential bloat on large faces)
        // We use 2 steps for tight resolution, 1 step for large macro textures
        const steps = params.pitch > 2.5 ? 1 : 2;
        const subPos = subdivideTriangles(islandPos, steps);
        const subCount = subPos.length / 9;

        const proj = createProjectionData(geometry, islandIndices, params.angle);
        if (!proj) return;

        // Calculate pitch variables based on whether the geometry is planar or cylindrical
        const nD_nom = Math.round(proj.circPhys / params.pitch);
        const nD = Math.max(2, nD_nom % 2 === 0 ? nD_nom : nD_nom + 1);
        const pitchP = proj.isPlanar ? params.pitch : (proj.circPhys / nD);

        // Define rotation angle transformation for the UV grid mapping
        const angRad = (params.angle * Math.PI) / 180;
        const factor = Math.abs(Math.cos(angRad)) + Math.abs(Math.sin(angRad));
        const pU = pitchP * factor > 0.01 ? pitchP * factor : pitchP;
        const pV = pitchP * factor > 0.01 ? pitchP * factor : pitchP;

        // Line intersection logic for Sutherland-Hodgman Polygon Clipping
        // This cuts overlapping 2D texture polygons strictly to the bounds of the original 3D triangles
        const intersect = (a: any, b: any, e1: any, e2: any) => {
            const da = (e2.u - e1.u) * (a.v - e1.v) - (e2.v - e1.v) * (a.u - e1.u), db = (e2.u - e1.u) * (b.v - e1.v) - (e2.v - e1.v) * (b.u - e1.u);
            const t = Math.abs(da) / (Math.max(1e-10, Math.abs(da) + Math.abs(db)));
            return { u: a.u + t * (b.u - a.u), v: a.v + t * (b.v - a.v), h: a.h + t * (b.h - a.h) };
        };
        const clip = (poly: any[], e1: any, e2: any) => {
            const out: any[] = [];
            for (let i = 0; i < poly.length; i++) {
                const a = poly[i], b = poly[(i + 1) % poly.length];
                const aIn = (e2.u - e1.u) * (a.v - e1.v) - (e2.v - e1.v) * (a.u - e1.u) >= -1e-8, bIn = (e2.u - e1.u) * (b.v - e1.v) - (e2.v - e1.v) * (b.u - e1.u) >= -1e-8;
                if (aIn) { if (bIn) out.push(b); else out.push(intersect(a, b, e1, e2)); } else if (bIn) { out.push(intersect(a, b, e1, e2)); out.push(b); }
            }
            return out;
        };

        // STEP 2: Iterate over all subdivided base triangles.
        for (let t_idx = 0; t_idx < subCount; t_idx++) {
            const off = t_idx * 9;
            const pA_vec = new THREE.Vector3(subPos[off], subPos[off + 1], subPos[off + 2]);
            const pB_vec = new THREE.Vector3(subPos[off + 3], subPos[off + 4], subPos[off + 5]);
            const pC_vec = new THREE.Vector3(subPos[off + 6], subPos[off + 7], subPos[off + 8]);

            // Map 3D points to 2D UV Space
            let pA_uv = proj.project(pA_vec), pB_uv = proj.project(pB_vec), pC_uv = proj.project(pC_vec);

            // Fix unrolling seam for cylindrical coordinates (ensure continuous UV space across 360 wrap)
            if (!proj.isPlanar) {
                const uVals = [pA_uv.u, pB_uv.u, pC_uv.u];
                const minU = Math.min(...uVals), maxU = Math.max(...uVals);
                if (maxU - minU > proj.circPhys * 0.5) {
                    if (pA_uv.u < proj.circPhys * 0.5) pA_uv.u += proj.circPhys;
                    if (pB_uv.u < proj.circPhys * 0.5) pB_uv.u += proj.circPhys;
                    if (pC_uv.u < proj.circPhys * 0.5) pC_uv.u += proj.circPhys;
                }
            }

            // Rotate mapped points to desired knurl angle (0, 30, 45, etc.) using `ur` and `vr`
            const rotA = proj.toRot(pA_uv.u, pA_uv.v), rotB = proj.toRot(pB_uv.u, pB_uv.v), rotC = proj.toRot(pC_uv.u, pC_uv.v);

            // Reorient triangulation to be Counter-Clockwise in UV space
            const area = (rotB.ur - rotA.ur) * (rotC.vr - rotA.vr) - (rotB.vr - rotA.vr) * (rotC.ur - rotA.ur);
            const tri = area < 0 ? [{ u: rotA.ur, v: rotA.vr }, { u: rotC.ur, v: rotC.vr }, { u: rotB.ur, v: rotB.vr }] : [{ u: rotA.ur, v: rotA.vr }, { u: rotB.ur, v: rotB.vr }, { u: rotC.ur, v: rotC.vr }];
            const triVecs = area < 0 ? [pA_vec, pC_vec, pB_vec] : [pA_vec, pB_vec, pC_vec];

            // STEP 3: Find bounding box of this triangle in the 2D Knurl Grid
            const iS = Math.floor(Math.min(tri[0].u, tri[1].u, tri[2].u) / pU) - 1, iE = Math.ceil(Math.max(tri[0].u, tri[1].u, tri[2].u) / pU) + 1;
            const jS = Math.floor(Math.min(tri[0].v, tri[1].v, tri[2].v) / pV) - 1, jE = Math.ceil(Math.max(tri[0].v, tri[1].v, tri[2].v) / pV) + 1;

            // Iterate through every applicable 2D pattern cell that intersects this triangle
            for (let j = jS; j <= jE; j++) {
                for (let i = iS; i <= iE; i++) {
                    const u0 = i * pU, u1 = (i + 1) * pU, v0 = j * pV, v1 = (j + 1) * pV, um = u0 + pU * 0.5, vm = v0 + pV * 0.5;

                    // Generate 2D parametric polygons based on the chosen pattern shape
                    // `h` determines the height map (0 is base, params.depth is peak)
                    let cell: { u: number; v: number; h: number }[][] = [];
                    if (params.pattern === 'diamond') {
                        const mid = { u: um, v: vm, h: ((Math.abs(i) + Math.abs(j)) % 2 === 0) ? params.depth : 0 };
                        cell = [[{ u: u0, v: v0, h: 0 }, { u: u1, v: v0, h: 0 }, mid], [{ u: u1, v: v0, h: 0 }, { u: u1, v: v1, h: 0 }, mid], [{ u: u1, v: v1, h: 0 }, { u: u0, v: v1, h: 0 }, mid], [{ u: u0, v: v1, h: 0 }, { u: u0, v: v0, h: 0 }, mid]];
                    } else if (params.pattern === 'straight' || params.pattern === 'diagonal') {
                        cell = [[{ u: u0, v: v0, h: 0 }, { u: um, v: v0, h: params.depth }, { u: u0, v: v1, h: 0 }], [{ u: um, v: v0, h: params.depth }, { u: um, v: v1, h: params.depth }, { u: u0, v: v1, h: 0 }], [{ u: um, v: v0, h: params.depth }, { u: u1, v: v0, h: 0 }, { u: um, v: v1, h: params.depth }], [{ u: u1, v: v0, h: 0 }, { u: u1, v: v1, h: 0 }, { u: um, v: v1, h: params.depth }]];
                    } else if (params.pattern === 'square') {
                        const mid = { u: um, v: vm, h: 0 };
                        cell = [[{ u: u0, v: v0, h: 0 }, { u: um, v: v0, h: params.depth }, mid], [{ u: um, v: v0, h: params.depth }, { u: u1, v: v0, h: 0 }, mid], [{ u: u1, v: v0, h: 0 }, { u: u1, v: vm, h: params.depth }, mid], [{ u: u1, v: vm, h: params.depth }, { u: u1, v: v1, h: 0 }, mid], [{ u: u1, v: v1, h: 0 }, { u: um, v: v1, h: params.depth }, mid], [{ u: um, v: v1, h: params.depth }, { u: u0, v: v1, h: 0 }, mid], [{ u: u0, v: v1, h: 0 }, { u: u0, v: vm, h: params.depth }, mid], [{ u: u0, v: vm, h: params.depth }, { u: u0, v: v0, h: 0 }, mid]];
                    }

                    // STEP 4: Process and extrude the generated shape
                    cell.forEach(kTri => {
                        // Sutherland-Hodgman clipping: Clip the 2D knurl polygon to the borders of the base 3D triangle
                        // This prevents texturing from bleeding outside the strict triangle boundaries.
                        let poly = kTri; poly = clip(poly, tri[0], tri[1]); poly = clip(poly, tri[1], tri[2]); poly = clip(poly, tri[2], tri[0]);

                        if (poly.length >= 3) {

                            // STEP 5: Extrude the Inner Face Geometry
                            // Convert the clipped 2D UV polygon back into 3D using Barycentric Coordinates
                            for (let v = 1; v < poly.length - 1; v++) {
                                const pts = [poly[0], poly[v], poly[v + 1]];
                                const tempTri: number[] = [-1, -1, -1];

                                pts.forEach((p, pIdx) => {
                                    const w = getBarycentric(p, tri[0], tri[1], tri[2]);
                                    const pBase = new THREE.Vector3().addScaledVector(triVecs[0], w.wA).addScaledVector(triVecs[1], w.wB).addScaledVector(triVecs[2], w.wC);
                                    const p0_u = proj.fromRot(p.u, p.v);
                                    const norm = proj.getNormal(p0_u.u, p0_u.v);
                                    const vFinal = pBase.addScaledVector(norm, (p as any).h || 0);

                                    const vIdx = getVertIndex(vFinal.x, vFinal.y, vFinal.z);

                                    if (proj.depthSign < 0) {
                                        // Reverse order for internal surfaces
                                        if (pIdx === 0) tempTri[0] = vIdx;
                                        else if (pIdx === 1) tempTri[1] = vIdx;
                                        else tempTri[2] = vIdx;
                                    } else {
                                        tempTri[pIdx] = vIdx;
                                    }
                                });

                                if (proj.depthSign < 0) {
                                    newIndices.push(tempTri[2], tempTri[1], tempTri[0]);
                                } else {
                                    newIndices.push(tempTri[0], tempTri[1], tempTri[2]);
                                }
                            }

                            // --- Outer Perimeter Wall Extrusion via Midpoints ---
                            for (let v = 0; v < poly.length; v++) {
                                const pA = poly[v], pB = poly[(v + 1) % poly.length];

                                const wA = getBarycentric(pA, tri[0], tri[1], tri[2]);
                                const vBaseA = new THREE.Vector3().addScaledVector(triVecs[0], wA.wA).addScaledVector(triVecs[1], wA.wB).addScaledVector(triVecs[2], wA.wC);

                                const wB = getBarycentric(pB, tri[0], tri[1], tri[2]);
                                const vBaseB = new THREE.Vector3().addScaledVector(triVecs[0], wB.wA).addScaledVector(triVecs[1], wB.wB).addScaledVector(triVecs[2], wB.wC);

                                const vMid = new THREE.Vector3().lerpVectors(vBaseA, vBaseB, 0.5);

                                if (isOnAnySegmentGrid(vMid, boundaryGrid, 1e-4)) {
                                    const pA_u = proj.fromRot(pA.u, pA.v);
                                    const pB_u = proj.fromRot(pB.u, pB.v);

                                    const nA = proj.getNormal(pA_u.u, pA_u.v);
                                    const nB = proj.getNormal(pB_u.u, pB_u.v);

                                    const valA = (pA as any).h !== undefined ? (pA as any).h : 0;
                                    const valB = (pB as any).h !== undefined ? (pB as any).h : 0;

                                    if (Math.abs(valA) < 1e-4 && Math.abs(valB) < 1e-4) continue;

                                    const vExtA = vBaseA.clone().addScaledVector(nA, valA);
                                    const vExtB = vBaseB.clone().addScaledVector(nB, valB);

                                    const idxBaseA = getVertIndex(vBaseA.x, vBaseA.y, vBaseA.z);
                                    const idxBaseB = getVertIndex(vBaseB.x, vBaseB.y, vBaseB.z);
                                    const idxExtA = getVertIndex(vExtA.x, vExtA.y, vExtA.z);
                                    const idxExtB = getVertIndex(vExtB.x, vExtB.y, vExtB.z);

                                    // For manifoldness, the wall base must be (BaseB, BaseA) 
                                    // because the selection face uses (BaseA, BaseB).
                                    newIndices.push(idxBaseB, idxBaseA, idxExtA);
                                    newIndices.push(idxBaseB, idxExtA, idxExtB);
                                    wallCount++;
                                }
                            }
                        }
                    });
                }
            }
        }
    });

    console.log(`[Knurling] Generated ${wallCount} wall quads.`);

    // Map all unselected (original base) faces into the same vertexMap
    const sel = new Set(faceIndices);
    const faceCount = index ? index.count / 3 : posAttr.count / 3;

    for (let f = 0; f < faceCount; f++) {
        if (!sel.has(f)) {
            const i0 = index ? index.getX(f * 3 + 0) : f * 3 + 0;
            const i1 = index ? index.getX(f * 3 + 1) : f * 3 + 1;
            const i2 = index ? index.getX(f * 3 + 2) : f * 3 + 2;

            const a = getVertIndex(posAttr.getX(i0), posAttr.getY(i0), posAttr.getZ(i0));
            const b = getVertIndex(posAttr.getX(i1), posAttr.getY(i1), posAttr.getZ(i1));
            const c = getVertIndex(posAttr.getX(i2), posAttr.getY(i2), posAttr.getZ(i2));
            newIndices.push(a, b, c);
        }
    }

    const finalGeo = new THREE.BufferGeometry();
    finalGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(vertices), 3));
    finalGeo.setIndex(newIndices);

    // Repair manifoldness if absolutely necessary, but topology *should* be perfect now
    const repaired = params.holeFillEnabled ? await repairWithManifold(finalGeo, { holeDiameterMultiplier: params.holeFillThreshold }) : finalGeo;
    repaired.computeVertexNormals();
    return repaired;
}

// --------------------------------------------------------------------------------
// HONEYCOMB IMPLEMENTATION
// --------------------------------------------------------------------------------

export async function applyHoneycomb(
    geometry: THREE.BufferGeometry,
    faceIndices: number[],
    params: HoneycombParams
): Promise<THREE.BufferGeometry> {
    if (!faceIndices || faceIndices.length === 0) return geometry;

    const index = geometry.index;
    const posAttr = geometry.attributes.position;
    const islands = getContiguousIslands(geometry, faceIndices);

    // --- TOPOLOGICAL INDEXING SETUP ---
    const vertexMap = new Map<string, number>();
    const vertices: number[] = [];
    const newIndices: number[] = [];

    const getVertIndex = (x: number, y: number, z: number): number => {
        const key = `${Math.round(x * 1e3)},${Math.round(y * 1e3)},${Math.round(z * 1e3)}`;
        let idx = vertexMap.get(key);
        if (idx === undefined) {
            idx = vertices.length / 3;
            vertexMap.set(key, idx);
            vertices.push(x, y, z);
        }
        return idx;
    };

    // Boundary segments for wall stitching
    const boundarySegments = getBoundarySegments(geometry, faceIndices);
    const boundaryGrid = createBoundaryGrid(boundarySegments, 2.0);
    console.log(`[Honeycomb] Found ${boundarySegments.length} boundary segments.`);
    let wallCount = 0;

    islands.forEach(islandIndices => {
        const islandPos = new Float32Array(islandIndices.length * 9);
        islandIndices.forEach((fIdx, i) => {
            const i0 = index ? index.getX(fIdx * 3 + 0) : fIdx * 3 + 0;
            const i1 = index ? index.getX(fIdx * 3 + 1) : fIdx * 3 + 1;
            const i2 = index ? index.getX(fIdx * 3 + 2) : fIdx * 3 + 2;
            islandPos[i * 9 + 0] = posAttr.getX(i0); islandPos[i * 9 + 1] = posAttr.getY(i0); islandPos[i * 9 + 2] = posAttr.getZ(i0);
            islandPos[i * 9 + 3] = posAttr.getX(i1); islandPos[i * 9 + 4] = posAttr.getY(i1); islandPos[i * 9 + 5] = posAttr.getZ(i1);
            islandPos[i * 9 + 6] = posAttr.getX(i2); islandPos[i * 9 + 7] = posAttr.getY(i2); islandPos[i * 9 + 8] = posAttr.getZ(i2);
        });

        // Step-based subdivision to prevent exponential face generation on large initial triangles
        const steps = params.cellSize > 3.0 ? 1 : 2;
        const subPos = subdivideTriangles(islandPos, steps);
        const subCount = subPos.length / 9;

        const proj = createProjectionData(geometry, islandIndices, params.angle);
        if (!proj) return;

        let { cellSize: W, wallThickness: t, depth, direction, angle: userAngle } = params;

        if (!proj.isPlanar && proj.circPhys > 0) {
            const ang = ((userAngle % 180) + 180) % 180;
            const symmetries = [0, 30, 60, 90, 120, 150, 180];
            const bestSym = symmetries.reduce((p, c) => Math.abs(c - ang) < Math.abs(p - ang) ? c : p);
            const periodFactor = (bestSym % 60 === 0) ? 1.0 : Math.sqrt(3);
            const nUnits = Math.round(proj.circPhys / (W * periodFactor));
            W = proj.circPhys / (Math.max(1, nUnits) * periodFactor);
        }

        const appliedDepth = direction === 'inward' ? -depth : depth;
        const s = (W / 2) / (Math.sqrt(3) / 2); // Side length
        const H = 1.5 * s; // Vertical distance between rows

        const intersect = (a: any, b: any, e1: any, e2: any) => {
            const da = (e2.u - e1.u) * (a.v - e1.v) - (e2.v - e1.v) * (a.u - e1.u), db = (e2.u - e1.u) * (b.v - e1.v) - (e2.v - e1.v) * (b.u - e1.u);
            const t = Math.abs(da) / (Math.max(1e-10, Math.abs(da) + Math.abs(db)));
            return { u: a.u + t * (b.u - a.u), v: a.v + t * (b.v - a.v), h: a.h + t * (b.h - a.h) };
        };
        const clip = (poly: any[], e1: any, e2: any) => {
            const out: any[] = [];
            for (let i = 0; i < poly.length; i++) {
                const a = poly[i], b = poly[(i + 1) % poly.length];
                const aIn = (e2.u - e1.u) * (a.v - e1.v) - (e2.v - e1.v) * (a.u - e1.u) >= -1e-8, bIn = (e2.u - e1.u) * (b.v - e1.v) - (e2.v - e1.v) * (b.u - e1.u) >= -1e-8;
                if (aIn) { if (bIn) out.push(b); else out.push(intersect(a, b, e1, e2)); } else if (bIn) { out.push(intersect(a, b, e1, e2)); out.push(b); }
            }
            return out;
        };

        for (let t_idx = 0; t_idx < subCount; t_idx++) {
            const off = t_idx * 9;
            const pA_vec = new THREE.Vector3(subPos[off], subPos[off + 1], subPos[off + 2]);
            const pB_vec = new THREE.Vector3(subPos[off + 3], subPos[off + 4], subPos[off + 5]);
            const pC_vec = new THREE.Vector3(subPos[off + 6], subPos[off + 7], subPos[off + 8]);
            let pA_uv = proj.project(pA_vec), pB_uv = proj.project(pB_vec), pC_uv = proj.project(pC_vec);
            if (!proj.isPlanar) {
                const uVals = [pA_uv.u, pB_uv.u, pC_uv.u];
                const minU = Math.min(...uVals), maxU = Math.max(...uVals);
                if (maxU - minU > proj.circPhys * 0.5) {
                    if (pA_uv.u < proj.circPhys * 0.5) pA_uv.u += proj.circPhys;
                    if (pB_uv.u < proj.circPhys * 0.5) pB_uv.u += proj.circPhys;
                    if (pC_uv.u < proj.circPhys * 0.5) pC_uv.u += proj.circPhys;
                }
            }
            const rotA = proj.toRot(pA_uv.u, pA_uv.v), rotB = proj.toRot(pB_uv.u, pB_uv.v), rotC = proj.toRot(pC_uv.u, pC_uv.v);
            const area = (rotB.ur - rotA.ur) * (rotC.vr - rotA.vr) - (rotB.vr - rotA.vr) * (rotC.ur - rotA.ur);
            const tri = area < 0 ? [{ u: rotA.ur, v: rotA.vr }, { u: rotC.ur, v: rotC.vr }, { u: rotB.ur, v: rotB.vr }] : [{ u: rotA.ur, v: rotA.vr }, { u: rotB.ur, v: rotB.vr }, { u: rotC.ur, v: rotC.vr }];
            const triVecs = area < 0 ? [pA_vec, pC_vec, pB_vec] : [pA_vec, pB_vec, pC_vec];

            const uMin = Math.min(tri[0].u, tri[1].u, tri[2].u), uMax = Math.max(tri[0].u, tri[1].u, tri[2].u);
            const vMin = Math.min(tri[0].v, tri[1].v, tri[2].v), vMax = Math.max(tri[0].v, tri[1].v, tri[2].v);

            const rIn = Math.max(0, (W / 2) - t);
            const sIn = rIn / (Math.sqrt(3) / 2);

            for (let j = Math.floor(vMin / H) - 1; j <= Math.ceil(vMax / H) + 1; j++) {
                for (let i = Math.floor(uMin / W) - 1; i <= Math.ceil(uMax / W) + 1; i++) {
                    const uc = (i + (Math.abs(j) % 2 === 1 ? 0.5 : 0)) * W;
                    const vc = j * H;

                    const vOut = [
                        { u: uc, v: vc + s, h: 0 },
                        { u: uc + W / 2, v: vc + s / 2, h: 0 },
                        { u: uc + W / 2, v: vc - s / 2, h: 0 },
                        { u: uc, v: vc - s, h: 0 },
                        { u: uc - W / 2, v: vc - s / 2, h: 0 },
                        { u: uc - W / 2, v: vc + s / 2, h: 0 }
                    ];

                    const vInn = [
                        { u: uc, v: vc + sIn, h: appliedDepth },
                        { u: uc + rIn, v: vc + sIn / 2, h: appliedDepth },
                        { u: uc + rIn, v: vc - sIn / 2, h: appliedDepth },
                        { u: uc, v: vc - sIn, h: appliedDepth },
                        { u: uc - rIn, v: vc - sIn / 2, h: appliedDepth },
                        { u: uc - rIn, v: vc + sIn / 2, h: appliedDepth }
                    ];

                    const hexTris: { u: number, v: number, h: number }[][] = [];
                    const center = { u: uc, v: vc, h: appliedDepth };

                    for (let k = 0; k < 6; k++) {
                        hexTris.push([vInn[k], vInn[(k + 1) % 6], center]);
                        hexTris.push([vOut[k], vOut[(k + 1) % 6], vInn[k]]);
                        hexTris.push([vOut[(k + 1) % 6], vInn[(k + 1) % 6], vInn[k]]);
                    }

                    hexTris.forEach(kTri => {
                        let poly = kTri; poly = clip(poly, tri[0], tri[1]); poly = clip(poly, tri[1], tri[2]); poly = clip(poly, tri[2], tri[0]);
                        if (poly.length >= 3) {
                            // STEP 5: Extrude the Inner Face Geometry
                            for (let v = 1; v < poly.length - 1; v++) {
                                const pts = [poly[0], poly[v], poly[v + 1]];
                                const tempTri: number[] = [-1, -1, -1];

                                pts.forEach((p, pIdx) => {
                                    const w = getBarycentric(p, tri[0], tri[1], tri[2]);
                                    const pBase = new THREE.Vector3().addScaledVector(triVecs[0], w.wA).addScaledVector(triVecs[1], w.wB).addScaledVector(triVecs[2], w.wC);
                                    const p0_u = proj.fromRot(p.u, p.v);
                                    const norm = proj.getNormal(p0_u.u, p0_u.v);
                                    const vFinal = pBase.addScaledVector(norm, (p as any).h || 0);

                                    const vIdx = getVertIndex(vFinal.x, vFinal.y, vFinal.z);

                                    if (proj.depthSign < 0) {
                                        if (pIdx === 0) tempTri[0] = vIdx;
                                        else if (pIdx === 1) tempTri[1] = vIdx;
                                        else tempTri[2] = vIdx;
                                    } else {
                                        tempTri[pIdx] = vIdx;
                                    }
                                });

                                if (proj.depthSign < 0) {
                                    newIndices.push(tempTri[2], tempTri[1], tempTri[0]);
                                } else {
                                    newIndices.push(tempTri[0], tempTri[1], tempTri[2]);
                                }
                            }

                            // STEP 6: Outer Perimeter Wall Extrusion via Exact Geometric Midpoints
                            // Similar to Knurling, we iterate the inner polygon structure near the boundaries.
                            // If the midpoint of an edge rests precisely on the user-selected mesh boundary,
                            // we duplicate those two vertices and extrude them down/up to seal the volume.
                            for (let v = 0; v < poly.length; v++) {
                                const pA = poly[v], pB = poly[(v + 1) % poly.length];

                                const wA = getBarycentric(pA, tri[0], tri[1], tri[2]);
                                const vBaseA = new THREE.Vector3().addScaledVector(triVecs[0], wA.wA).addScaledVector(triVecs[1], wA.wB).addScaledVector(triVecs[2], wA.wC);

                                const wB = getBarycentric(pB, tri[0], tri[1], tri[2]);
                                const vBaseB = new THREE.Vector3().addScaledVector(triVecs[0], wB.wA).addScaledVector(triVecs[1], wB.wB).addScaledVector(triVecs[2], wB.wC);

                                const vMid = new THREE.Vector3().lerpVectors(vBaseA, vBaseB, 0.5);

                                if (isOnAnySegmentGrid(vMid, boundaryGrid, 1e-4)) {
                                    const pA_u = proj.fromRot(pA.u, pA.v);
                                    const pB_u = proj.fromRot(pB.u, pB.v);

                                    const nA = proj.getNormal(pA_u.u, pA_u.v);
                                    const nB = proj.getNormal(pB_u.u, pB_u.v);

                                    // Gather exact extruding parameters (normal and target height)
                                    // Height is `valA` and `valB`, which depends on whether the hexagon edge
                                    // belongs to the inner cavity or outer perimeter
                                    const valA = (pA as any).h !== undefined ? (pA as any).h : 0;
                                    const valB = (pB as any).h !== undefined ? (pB as any).h : 0;

                                    if (Math.abs(valA) < 1e-4 && Math.abs(valB) < 1e-4) continue;

                                    const vExtA = vBaseA.clone().addScaledVector(nA, valA);
                                    const vExtB = vBaseB.clone().addScaledVector(nB, valB);

                                    const idxBaseA = getVertIndex(vBaseA.x, vBaseA.y, vBaseA.z);
                                    const idxBaseB = getVertIndex(vBaseB.x, vBaseB.y, vBaseB.z);
                                    const idxExtA = getVertIndex(vExtA.x, vExtA.y, vExtA.z);
                                    const idxExtB = getVertIndex(vExtB.x, vExtB.y, vExtB.z);

                                    // For manifoldness, the wall base must be (BaseB, BaseA) 
                                    // because the selection face uses (BaseA, BaseB).
                                    newIndices.push(idxBaseB, idxBaseA, idxExtA);
                                    newIndices.push(idxBaseB, idxExtA, idxExtB);
                                    wallCount++;
                                }
                            }
                        }
                    });
                }
            }
        }
    });

    console.log(`[Honeycomb] Generated ${wallCount} wall quads.`);

    // Map all unselected (original base) faces into the same vertexMap
    const sel = new Set(faceIndices);
    const faceCount = index ? index.count / 3 : posAttr.count / 3;

    for (let f = 0; f < faceCount; f++) {
        if (!sel.has(f)) {
            const i0 = index ? index.getX(f * 3 + 0) : f * 3 + 0;
            const i1 = index ? index.getX(f * 3 + 1) : f * 3 + 1;
            const i2 = index ? index.getX(f * 3 + 2) : f * 3 + 2;

            const a = getVertIndex(posAttr.getX(i0), posAttr.getY(i0), posAttr.getZ(i0));
            const b = getVertIndex(posAttr.getX(i1), posAttr.getY(i1), posAttr.getZ(i1));
            const c = getVertIndex(posAttr.getX(i2), posAttr.getY(i2), posAttr.getZ(i2));
            newIndices.push(a, b, c);
        }
    }

    const finalGeo = new THREE.BufferGeometry();
    finalGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(vertices), 3));
    finalGeo.setIndex(newIndices);

    const repaired = params.holeFillEnabled ? await repairWithManifold(finalGeo, { holeDiameterMultiplier: params.holeFillThreshold }) : finalGeo;
    repaired.computeVertexNormals();
    return repaired;
}

// --------------------------------------------------------------------------------
// FUZZY SKIN IMPLEMENTATION
// --------------------------------------------------------------------------------

export async function applyFuzzySkin(
    geometry: THREE.BufferGeometry,
    faceIndices: number[],
    params: FuzzySkinParams
): Promise<THREE.BufferGeometry> {
    const islands = getContiguousIslands(geometry, faceIndices);
    const posAttr = geometry.attributes.position;

    // --- TOPOLOGICAL INDEXING SETUP ---
    const vertexMap = new Map<string, number>();
    const vertices: number[] = [];
    const newIndices: number[] = [];

    const getVertIndex = (x: number, y: number, z: number): number => {
        const key = `${Math.round(x * 1e3)},${Math.round(y * 1e3)},${Math.round(z * 1e3)}`;
        let idx = vertexMap.get(key);
        if (idx === undefined) {
            idx = vertices.length / 3;
            vertexMap.set(key, idx);
            vertices.push(x, y, z);
        }
        return idx;
    };

    // Fuzzy Skin does not need complex perimeter wall generation like knurling/honeycomb.
    // Instead, to keep boundaries watertight with the original, we just CLAMP the applied noise
    // to 0 for any points that fall smoothly within the fuzzy skin edge seam.
    const boundarySegments = getBoundarySegments(geometry, faceIndices);
    const boundaryGrid = createBoundaryGrid(boundarySegments, 2.0);
    const clampDist = params.holeFillThreshold;

    const index = geometry.index;
    islands.forEach(islandIndices => {
        const islandPos = new Float32Array(islandIndices.length * 9);
        islandIndices.forEach((fIdx, i) => {
            const i0 = index ? index.getX(fIdx * 3 + 0) : fIdx * 3 + 0;
            const i1 = index ? index.getX(fIdx * 3 + 1) : fIdx * 3 + 1;
            const i2 = index ? index.getX(fIdx * 3 + 2) : fIdx * 3 + 2;
            islandPos[i * 9 + 0] = posAttr.getX(i0); islandPos[i * 9 + 1] = posAttr.getY(i0); islandPos[i * 9 + 2] = posAttr.getZ(i0);
            islandPos[i * 9 + 3] = posAttr.getX(i1); islandPos[i * 9 + 4] = posAttr.getY(i1); islandPos[i * 9 + 5] = posAttr.getZ(i1);
            islandPos[i * 9 + 6] = posAttr.getX(i2); islandPos[i * 9 + 7] = posAttr.getY(i2); islandPos[i * 9 + 8] = posAttr.getZ(i2);
        });

        // 1. Homogeneous Density Calculation
        let currentSubPos: Float32Array<ArrayBufferLike> = islandPos;
        // Cap the point distance to avoid infinite/crashing subdivisions
        const effectivePtDist = Math.max(0.1, params.pointDistance);
        const targetSq = effectivePtDist * effectivePtDist;
        const MAX_STEPS = 6;

        for (let s = 0; s < MAX_STEPS; s++) {
            let needsMore = false;
            const triCount = currentSubPos.length / 9;

            // Safety cap: don't subdivide if we already have > 1,000,000 triangles to prevent memory crash
            if (triCount > 1000000) break;

            for (let i = 0; i < triCount; i++) {
                const o = i * 9;
                const d12 = distSq(currentSubPos[o], currentSubPos[o + 1], currentSubPos[o + 2], currentSubPos[o + 3], currentSubPos[o + 4], currentSubPos[o + 5]);
                const d23 = distSq(currentSubPos[o + 3], currentSubPos[o + 4], currentSubPos[o + 5], currentSubPos[o + 6], currentSubPos[o + 7], currentSubPos[o + 8]);
                const d31 = distSq(currentSubPos[o + 6], currentSubPos[o + 7], currentSubPos[o + 8], currentSubPos[o], currentSubPos[o + 1], currentSubPos[o + 2]);
                if (Math.max(d12, d23, d31) > targetSq) {
                    needsMore = true;
                    break;
                }
            }
            if (!needsMore) break;
            currentSubPos = subdivideTriangles(currentSubPos, 1) as Float32Array;
        }

        // 2. Mesh Context for jitter
        let subGeo = new THREE.BufferGeometry();
        subGeo.setAttribute('position', new THREE.BufferAttribute(currentSubPos, 3));
        subGeo = weldGeometry(subGeo, 10000);
        subGeo.computeVertexNormals();

        const positions = subGeo.attributes.position.array as any as Float32Array;
        const normals = subGeo.attributes.normal.array as any as Float32Array;

        // 3. Application of Noise (Random Displacement) with Boundary Clamp
        for (let i = 0; i < positions.length / 3; i++) {
            let noise = (Math.random() - 0.5) * 2 * params.thickness;
            const px = positions[i * 3], py = positions[i * 3 + 1], pz = positions[i * 3 + 2];
            const p = new THREE.Vector3(px, py, pz);

            // Check border collision using the high-performance Spatial Hash grid.
            // If the point is within the `clampDist` distance of a boundary, we apply 0 noise
            // so the generated geometry remains perfectly flush and manifold with the rest of the model.
            if (isOnAnySegmentGrid(p, boundaryGrid, clampDist)) {
                noise = 0;
            }

            positions[i * 3] += normals[i * 3] * noise;
            positions[i * 3 + 1] += normals[i * 3 + 1] * noise;
            positions[i * 3 + 2] += normals[i * 3 + 2] * noise;
        }

        // 4. Conversion to final indexed geometry
        const finalIslandPos = subGeo.attributes.position.array as any as Float32Array;
        const finalIslandIdx = subGeo.getIndex()?.array || Array.from({ length: finalIslandPos.length / 3 }, (_, k) => k);

        for (let i = 0; i < finalIslandIdx.length; i += 3) {
            const i0 = finalIslandIdx[i] * 3, i1 = finalIslandIdx[i + 1] * 3, i2 = finalIslandIdx[i + 2] * 3;

            const px0 = finalIslandPos[i0], py0 = finalIslandPos[i0 + 1], pz0 = finalIslandPos[i0 + 2];
            const px1 = finalIslandPos[i1], py1 = finalIslandPos[i1 + 1], pz1 = finalIslandPos[i1 + 2];
            const px2 = finalIslandPos[i2], py2 = finalIslandPos[i2 + 1], pz2 = finalIslandPos[i2 + 2];

            newIndices.push(
                getVertIndex(px0, py0, pz0),
                getVertIndex(px1, py1, pz1),
                getVertIndex(px2, py2, pz2)
            );
        }
    });

    // Map all unselected (original base) faces into the same vertexMap
    const sel = new Set(faceIndices);
    const faceCount = index ? index.count / 3 : posAttr.count / 3;

    for (let f = 0; f < faceCount; f++) {
        if (!sel.has(f)) {
            const i0 = index ? index.getX(f * 3 + 0) : f * 3 + 0;
            const i1 = index ? index.getX(f * 3 + 1) : f * 3 + 1;
            const i2 = index ? index.getX(f * 3 + 2) : f * 3 + 2;

            const a = getVertIndex(posAttr.getX(i0), posAttr.getY(i0), posAttr.getZ(i0));
            const b = getVertIndex(posAttr.getX(i1), posAttr.getY(i1), posAttr.getZ(i1));
            const c = getVertIndex(posAttr.getX(i2), posAttr.getY(i2), posAttr.getZ(i2));
            newIndices.push(a, b, c);
        }
    }

    const finalGeo = new THREE.BufferGeometry();
    finalGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(vertices), 3));
    finalGeo.setIndex(newIndices);

    const repaired = params.holeFillEnabled ? await repairWithManifold(finalGeo, { holeDiameterMultiplier: params.holeFillThreshold }) : finalGeo;
    repaired.computeVertexNormals();
    return repaired;
}

export function applyDecimate(geo: THREE.BufferGeometry, _faceIndices: number[], _params: any) { return geo; }
