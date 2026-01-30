import * as THREE from 'three';
import { getContiguousIslands, subdivideTriangles, fitCylinderToSelection } from './meshUtils';

export type KnurlPattern = 'diamond' | 'straight' | 'diagonal' | 'square';

export interface KnurlingParams {
    type: 'knurling';
    pitch: number;
    depth: number;
    angle: number;
    pattern: KnurlPattern;
}

export interface HoneycombParams {
    type: 'honeycomb';
    cellSize: number;
    wallThickness: number;
    depth: number;
    angle: number;
    direction: 'inward' | 'outward';
}

// --------------------------------------------------------------------------------
// SHARED UTILITIES
// --------------------------------------------------------------------------------

function getVKey(x: number, y: number, z: number): string {
    return `${Math.round(x * 100000)},${Math.round(y * 100000)},${Math.round(z * 100000)}`;
}

function weldGeometry(geo: THREE.BufferGeometry, prec: number = 100000): THREE.BufferGeometry {
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

function repairMesh(geo: THREE.BufferGeometry): THREE.BufferGeometry {
    let welded = weldGeometry(geo, 100000);
    const pos = welded.attributes.position, vertexCount = pos.count, edgeCounts = new Map<string, number>();
    for (let i = 0; i < vertexCount; i += 3) {
        for (let j = 0; j < 3; j++) {
            const v1 = i + j, v2 = i + (j + 1) % 3, h1 = getVKey(pos.getX(v1), pos.getY(v1), pos.getZ(v1)), h2 = getVKey(pos.getX(v2), pos.getY(v2), pos.getZ(v2));
            const edgeKey = h1 < h2 ? `${h1}|${h2}` : `${h2}|${h1}`;
            edgeCounts.set(edgeKey, (edgeCounts.get(edgeKey) || 0) + 1);
        }
    }
    const bIdx = new Set<number>();
    for (let i = 0; i < vertexCount; i += 3) {
        for (let j = 0; j < 3; j++) {
            const v1 = i + j, v2 = i + (j + 1) % 3, h1 = getVKey(pos.getX(v1), pos.getY(v1), pos.getZ(v1)), h2 = getVKey(pos.getX(v2), pos.getY(v2), pos.getZ(v2));
            const edgeKey = h1 < h2 ? `${h1}|${h2}` : `${h2}|${h1}`;
            if (edgeCounts.get(edgeKey) === 1) { bIdx.add(v1); bIdx.add(v2); }
        }
    }
    if (bIdx.size > 0) {
        const snapT = 0.2, positions = welded.attributes.position.array as Float32Array, grid = new Map<string, number>();
        bIdx.forEach(idx => {
            const x = positions[idx * 3], y = positions[idx * 3 + 1], z = positions[idx * 3 + 2];
            const gx = Math.floor(x / snapT), gy = Math.floor(y / snapT), gz = Math.floor(z / snapT);
            let snapped = false;
            for (let dx = -1; dx <= 1 && !snapped; dx++) {
                for (let dy = -1; dy <= 1 && !snapped; dy++) {
                    for (let dz = -1; dz <= 1 && !snapped; dz++) {
                        const key = `${gx + dx},${gy + dy},${gz + dz}`;
                        if (grid.has(key)) {
                            const tIdx = grid.get(key)!; const tx = positions[tIdx * 3], ty = positions[tIdx * 3 + 1], tz = positions[tIdx * 3 + 2];
                            if ((x - tx) ** 2 + (y - ty) ** 2 + (z - tz) ** 2 < snapT ** 2) {
                                positions[idx * 3] = tx; positions[idx * 3 + 1] = ty; positions[idx * 3 + 2] = tz;
                                snapped = true;
                            }
                        }
                    }
                }
            }
            if (!snapped) grid.set(`${gx},${gy},${gz}`, idx);
        });
        welded = weldGeometry(welded, 100000);
    }
    return welded;
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

    return { avgR, circPhys, depthSign, isPlanar, project, unproject, toRot, fromRot };
}

// --------------------------------------------------------------------------------
// KNURLING IMPLEMENTATION
// --------------------------------------------------------------------------------

export function applyKnurling(
    geometry: THREE.BufferGeometry,
    faceIndices: number[],
    params: KnurlingParams
): THREE.BufferGeometry {
    const proj = createProjectionData(geometry, faceIndices, params.angle);
    if (!proj) return geometry;

    const posAttr = geometry.attributes.position;
    const islands = getContiguousIslands(geometry, faceIndices);
    const newPos: number[] = [];

    islands.forEach(islandIndices => {
        const islandPos = new Float32Array(islandIndices.length * 9);
        islandIndices.forEach((fIdx, i) => {
            const off = fIdx * 3;
            islandPos[i * 9 + 0] = posAttr.getX(off); islandPos[i * 9 + 1] = posAttr.getY(off); islandPos[i * 9 + 2] = posAttr.getZ(off);
            islandPos[i * 9 + 3] = posAttr.getX(off + 1); islandPos[i * 9 + 4] = posAttr.getY(off + 1); islandPos[i * 9 + 5] = posAttr.getZ(off + 1);
            islandPos[i * 9 + 6] = posAttr.getX(off + 2); islandPos[i * 9 + 7] = posAttr.getY(off + 2); islandPos[i * 9 + 8] = posAttr.getZ(off + 2);
        });
        const subPos = subdivideTriangles(islandPos, 2);
        const subCount = subPos.length / 9;

        const nD_nom = Math.round(proj.circPhys / params.pitch);
        const nD = Math.max(2, nD_nom % 2 === 0 ? nD_nom : nD_nom + 1);
        const pitchP = proj.isPlanar ? params.pitch : (proj.circPhys / nD);
        const angRad = (params.angle * Math.PI) / 180;
        const factor = Math.cos(angRad) + Math.sin(angRad);
        const pU = pitchP * factor, pV = pitchP * factor;

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

        for (let t = 0; t < subCount; t++) {
            const off = t * 9;
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

            const iS = Math.floor(Math.min(tri[0].u, tri[1].u, tri[2].u) / pU) - 1, iE = Math.ceil(Math.max(tri[0].u, tri[1].u, tri[2].u) / pU) + 1;
            const jS = Math.floor(Math.min(tri[0].v, tri[1].v, tri[2].v) / pV) - 1, jE = Math.ceil(Math.max(tri[0].v, tri[1].v, tri[2].v) / pV) + 1;

            for (let j = jS; j <= jE; j++) {
                for (let i = iS; i <= iE; i++) {
                    const u0 = i * pU, u1 = (i + 1) * pU, v0 = j * pV, v1 = (j + 1) * pV, um = u0 + pU * 0.5, vm = v0 + pV * 0.5;
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

                    cell.forEach(kTri => {
                        let poly = kTri; poly = clip(poly, tri[0], tri[1]); poly = clip(poly, tri[1], tri[2]); poly = clip(poly, tri[2], tri[0]);
                        if (poly.length >= 3) {
                            for (let v = 1; v < poly.length - 1; v++) {
                                const p0_u = proj.fromRot(poly[0].u, poly[0].v), p1_u = proj.fromRot(poly[v].u, poly[v].v), p2_u = proj.fromRot(poly[v + 1].u, poly[v + 1].v);
                                const vA = proj.unproject(p0_u.u, p0_u.v, poly[0].h || 0), vB = proj.unproject(p1_u.u, p1_u.v, poly[v].h || 0), vC = proj.unproject(p2_u.u, p2_u.v, poly[v + 1].h || 0);
                                newPos.push(vA.x, vA.y, vA.z, vB.x, vB.y, vB.z, vC.x, vC.y, vC.z);
                            }
                        }
                    });
                }
            }
        }
    });

    const finalR: number[] = []; const sel = new Set(faceIndices);
    for (let f = 0; f < posAttr.count / 3; f++) { if (!sel.has(f)) { for (let v = 0; v < 3; v++) { const idx = f * 3 + v; finalR.push(posAttr.getX(idx), posAttr.getY(idx), posAttr.getZ(idx)); } } }
    const combined = new Float32Array(finalR.length + newPos.length); combined.set(finalR); combined.set(newPos, finalR.length);
    const finalGeo = new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(combined, 3));
    const repaired = repairMesh(finalGeo); repaired.computeVertexNormals(); return repaired;
}

// --------------------------------------------------------------------------------
// HONEYCOMB IMPLEMENTATION
// --------------------------------------------------------------------------------

export function applyHoneycomb(
    geometry: THREE.BufferGeometry,
    faceIndices: number[],
    params: HoneycombParams
): THREE.BufferGeometry {
    const proj = createProjectionData(geometry, faceIndices, params.angle);
    if (!proj) return geometry;

    let { cellSize: W, wallThickness: t, depth, direction, angle: userAngle } = params;

    // 1. Technical Wrap Adjustment (Cylinders/Cones)
    if (!proj.isPlanar && proj.circPhys > 0) {
        // Find nearest technical symmetry angle
        const ang = ((userAngle % 180) + 180) % 180;
        const symmetries = [0, 30, 60, 90, 120, 150, 180];
        const bestSym = symmetries.reduce((p, c) => Math.abs(c - ang) < Math.abs(p - ang) ? c : p);

        // If we are close to a symmetry angle or if we want a clean wrap anyway
        // For honeycomb, the period along the circumference depends on the orientation:
        // - Flat edges vertical (0, 60, 120): Period = cellSize (W)
        // - Points vertical (30, 90, 150): Period = W * sqrt(3)
        const periodFactor = (bestSym % 60 === 0) ? 1.0 : Math.sqrt(3);

        const nUnits = Math.round(proj.circPhys / (W * periodFactor));
        W = proj.circPhys / (Math.max(1, nUnits) * periodFactor);
    }

    const appliedDepth = direction === 'inward' ? -depth : depth;
    // Hexagon height from center to vertices
    const s = (W / 2) / (Math.sqrt(3) / 2); // Side length
    const H = 1.5 * s; // Vertical distance between rows

    const posAttr = geometry.attributes.position;
    const islands = getContiguousIslands(geometry, faceIndices);
    const newPos: number[] = [];

    islands.forEach(islandIndices => {
        const islandPos = new Float32Array(islandIndices.length * 9);
        islandIndices.forEach((fIdx, i) => {
            const off = fIdx * 3;
            islandPos[i * 9 + 0] = posAttr.getX(off); islandPos[i * 9 + 1] = posAttr.getY(off); islandPos[i * 9 + 2] = posAttr.getZ(off);
            islandPos[i * 9 + 3] = posAttr.getX(off + 1); islandPos[i * 9 + 4] = posAttr.getY(off + 1); islandPos[i * 9 + 5] = posAttr.getZ(off + 1);
            islandPos[i * 9 + 6] = posAttr.getX(off + 2); islandPos[i * 9 + 7] = posAttr.getY(off + 2); islandPos[i * 9 + 8] = posAttr.getZ(off + 2);
        });
        const subPos = subdivideTriangles(islandPos, 2);
        const subCount = subPos.length / 9;

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

            const uMin = Math.min(tri[0].u, tri[1].u, tri[2].u), uMax = Math.max(tri[0].u, tri[1].u, tri[2].u);
            const vMin = Math.min(tri[0].v, tri[1].v, tri[2].v), vMax = Math.max(tri[0].v, tri[1].v, tri[2].v);

            // Honeycomb inner radius (for the walls)
            const rIn = Math.max(0, (W / 2) - t);
            const sIn = rIn / (Math.sqrt(3) / 2);

            for (let j = Math.floor(vMin / H) - 1; j <= Math.ceil(vMax / H) + 1; j++) {
                for (let i = Math.floor(uMin / W) - 1; i <= Math.ceil(uMax / W) + 1; i++) {
                    const uc = (i + (Math.abs(j) % 2 === 1 ? 0.5 : 0)) * W;
                    const vc = j * H;

                    // Hexagon outer vertices
                    const vOut = [
                        { u: uc, v: vc + s, h: 0 },
                        { u: uc + W / 2, v: vc + s / 2, h: 0 },
                        { u: uc + W / 2, v: vc - s / 2, h: 0 },
                        { u: uc, v: vc - s, h: 0 },
                        { u: uc - W / 2, v: vc - s / 2, h: 0 },
                        { u: uc - W / 2, v: vc + s / 2, h: 0 }
                    ];

                    // Hexagon inner vertices (displaced)
                    const vInn = [
                        { u: uc, v: vc + sIn, h: appliedDepth },
                        { u: uc + rIn, v: vc + sIn / 2, h: appliedDepth },
                        { u: uc + rIn, v: vc - sIn / 2, h: appliedDepth },
                        { u: uc, v: vc - sIn, h: appliedDepth },
                        { u: uc - rIn, v: vc - sIn / 2, h: appliedDepth },
                        { u: uc - rIn, v: vc + sIn / 2, h: appliedDepth }
                    ];

                    // Triangles for the hex pocket
                    const hexTris: { u: number, v: number, h: number }[][] = [];
                    const center = { u: uc, v: vc, h: appliedDepth };

                    for (let k = 0; k < 6; k++) {
                        // Floor
                        hexTris.push([vInn[k], vInn[(k + 1) % 6], center]);
                        // Walls
                        hexTris.push([vOut[k], vOut[(k + 1) % 6], vInn[k]]);
                        hexTris.push([vOut[(k + 1) % 6], vInn[(k + 1) % 6], vInn[k]]);
                    }

                    hexTris.forEach(kTri => {
                        let poly = kTri; poly = clip(poly, tri[0], tri[1]); poly = clip(poly, tri[1], tri[2]); poly = clip(poly, tri[2], tri[0]);
                        if (poly.length >= 3) {
                            for (let v = 1; v < poly.length - 1; v++) {
                                const p0_u = proj.fromRot(poly[0].u, poly[0].v), p1_u = proj.fromRot(poly[v].u, poly[v].v), p2_u = proj.fromRot(poly[v + 1].u, poly[v + 1].v);
                                const vA = proj.unproject(p0_u.u, p0_u.v, poly[0].h || 0), vB = proj.unproject(p1_u.u, p1_u.v, poly[v].h || 0), vC = proj.unproject(p2_u.u, p2_u.v, poly[v + 1].h || 0);
                                newPos.push(vA.x, vA.y, vA.z, vB.x, vB.y, vB.z, vC.x, vC.y, vC.z);
                            }
                        }
                    });
                }
            }
        }
    });

    const finalR: number[] = []; const sel = new Set(faceIndices);
    for (let f = 0; f < posAttr.count / 3; f++) { if (!sel.has(f)) { for (let v = 0; v < 3; v++) { const idx = f * 3 + v; finalR.push(posAttr.getX(idx), posAttr.getY(idx), posAttr.getZ(idx)); } } }
    const combined = new Float32Array(finalR.length + newPos.length); combined.set(finalR); combined.set(newPos, finalR.length);
    const finalGeo = new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(combined, 3));
    const repaired = repairMesh(finalGeo); repaired.computeVertexNormals(); return repaired;
}

export function applyDecimate(geo: THREE.BufferGeometry, _faceIndices: number[], _params: any) { return geo; }
