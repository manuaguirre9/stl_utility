import * as THREE from 'three';
import { getContiguousIslands, subdivideTriangles, fitCylinderToSelection } from './meshUtils';

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

function getBarycentric(p: { u: number, v: number }, a: { u: number, v: number }, b: { u: number, v: number }, c: { u: number, v: number }) {
    const v0u = b.u - a.u, v0v = b.v - a.v, v1u = c.u - a.u, v1v = c.v - a.v, v2u = p.u - a.u, v2v = p.v - a.v;
    const den = v0u * v1v - v1u * v0v;
    if (Math.abs(den) < 1e-12) return { wA: 1, wB: 0, wC: 0 };
    const v = (v2u * v1v - v1u * v2v) / den;
    const w = (v0u * v2v - v2u * v0v) / den;
    return { wA: 1 - v - w, wB: v, wC: w };
}

function repairMesh(geo: THREE.BufferGeometry, threshold: number, enabled: boolean): THREE.BufferGeometry {
    // Pass 1: Initial weld to bridge microscopic gaps
    let welded = weldGeometry(geo, 10000);
    let pos = welded.attributes.position;

    const findBoundaries = (positions: THREE.BufferAttribute | THREE.InterleavedBufferAttribute) => {
        const edgeMap = new Map<string, { a: number[], b: number[], count: number }>();
        for (let i = 0; i < positions.count; i += 3) {
            for (let j = 0; j < 3; j++) {
                const idx1 = i + j;
                const idx2 = i + (j + 1) % 3;
                const v1 = [positions.getX(idx1), positions.getY(idx1), positions.getZ(idx1)];
                const v2 = [positions.getX(idx2), positions.getY(idx2), positions.getZ(idx2)];
                const h1 = getVKey(v1[0], v1[1], v1[2]);
                const h2 = getVKey(v2[0], v2[1], v2[2]);
                const key = h1 < h2 ? `${h1}|${h2}` : `${h2}|${h1}`;
                const entry = edgeMap.get(key) || { a: v1, b: v2, count: 0 };
                entry.count++;
                edgeMap.set(key, entry);
            }
        }
        const boundaryEdges: { a: number[], b: number[] }[] = [];
        edgeMap.forEach(e => { if (e.count === 1) boundaryEdges.push({ a: e.a, b: e.b }); });
        return boundaryEdges;
    };

    // Pass 2: Snapping boundary vertices to close small T-junctions
    let bEdges = findBoundaries(pos);
    if (bEdges.length > 0) {
        const snapT = 0.15;
        const positions = welded.attributes.position.array as Float32Array;
        const grid = new Map<string, number[]>();

        // Index all boundary vertices
        const bVerts = new Set<string>();
        bEdges.forEach(e => { bVerts.add(getVKey(e.a[0], e.a[1], e.a[2])); bVerts.add(getVKey(e.b[0], e.b[1], e.b[2])); });

        for (let i = 0; i < pos.count; i++) {
            const vx = positions[i * 3], vy = positions[i * 3 + 1], vz = positions[i * 3 + 2];
            const vh = getVKey(vx, vy, vz);
            if (!bVerts.has(vh)) continue;

            const gx = Math.floor(vx / snapT), gy = Math.floor(vy / snapT), gz = Math.floor(vz / snapT);
            const key = `${gx},${gy},${gz}`;
            const cell = grid.get(key) || [];
            let snapped = false;
            for (const otherIdx of cell) {
                const ox = positions[otherIdx * 3], oy = positions[otherIdx * 3 + 1], oz = positions[otherIdx * 3 + 2];
                if ((vx - ox) ** 2 + (vy - oy) ** 2 + (vz - oz) ** 2 < snapT * snapT) {
                    positions[i * 3] = ox; positions[i * 3 + 1] = oy; positions[i * 3 + 2] = oz;
                    snapped = true; break;
                }
            }
            if (!snapped) cell.push(i);
            grid.set(key, cell);
        }
        welded = weldGeometry(welded, 10000);
        pos = welded.attributes.position;
    }

    // Pass 3: Hole Filling - detect and triangulate open loops
    if (enabled) {
        bEdges = findBoundaries(pos);
        if (bEdges.length > 3) {
            const adj = new Map<string, string[]>();
            bEdges.forEach(e => {
                const h1 = getVKey(e.a[0], e.a[1], e.a[2]), h2 = getVKey(e.b[0], e.b[1], e.b[2]);
                const a1 = adj.get(h1) || []; a1.push(h2); adj.set(h1, a1);
                const a2 = adj.get(h2) || []; a2.push(h1); adj.set(h2, a2);
            });

            const loops: string[][] = [];
            const visited = new Set<string>();
            const hToV = new Map<string, number[]>();
            bEdges.forEach(e => {
                hToV.set(getVKey(e.a[0], e.a[1], e.a[2]), e.a);
                hToV.set(getVKey(e.b[0], e.b[1], e.b[2]), e.b);
            });

            adj.forEach((neighbors, startNode) => {
                if (visited.has(startNode) || neighbors.length !== 2) return;
                const loop = [startNode];
                visited.add(startNode);
                let curr = neighbors[0];
                while (curr && curr !== startNode && !visited.has(curr)) {
                    visited.add(curr);
                    loop.push(curr);
                    const nexts = adj.get(curr) || [];
                    curr = nexts.find(n => !visited.has(n)) || (nexts.length === 2 && nexts.includes(startNode) ? startNode : "");
                }
                if (loop.length >= 3) loops.push(loop);
            });

            if (loops.length > 0) {
                const extraTris: number[] = [];
                loops.forEach(loop => {
                    const points = loop.map(h => hToV.get(h)!);
                    const center = [0, 0, 0];
                    points.forEach(p => { center[0] += p[0]; center[1] += p[1]; center[2] += p[2]; });
                    center[0] /= points.length; center[1] /= points.length; center[2] /= points.length;

                    // Only fill small holes (e.g., radius < 2mm)
                    let maxDistSq = 0;
                    for (const p of points) {
                        const dSq = (p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2 + (p[2] - center[2]) ** 2;
                        if (dSq > maxDistSq) maxDistSq = dSq;
                    }

                    if (maxDistSq < threshold * threshold) { // Configurable threshold
                        for (let i = 0; i < points.length; i++) {
                            const p1 = points[i], p2 = points[(i + 1) % points.length];
                            extraTris.push(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], center[0], center[1], center[2]);
                        }
                    }
                });
                if (extraTris.length > 0) {
                    const oldPos = welded.attributes.position.array as Float32Array;
                    const combined = new Float32Array(oldPos.length + extraTris.length);
                    combined.set(oldPos); combined.set(extraTris, oldPos.length);
                    welded = new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(combined, 3));
                }
            }
        }
    }

    return weldGeometry(welded, 10000);
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
        // True Normal for a cone (r = m*v + b) is (1, -m) in (radial, up) space
        const m = fit.m || 0;
        const norm = radial.clone().addScaledVector(upDir, -m).normalize();
        return norm.multiplyScalar(depthSign);
    };

    return { avgR, circPhys, depthSign, isPlanar, project, unproject, toRot, fromRot, getNormal };
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
            const triVecs = area < 0 ? [pA_vec, pC_vec, pB_vec] : [pA_vec, pB_vec, pC_vec];

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
                                const pts = [poly[0], poly[v], poly[v + 1]];
                                const tempTri: THREE.Vector3[] = [new THREE.Vector3(), new THREE.Vector3(), new THREE.Vector3()];
                                pts.forEach((p, pIdx) => {
                                    const w = getBarycentric(p, tri[0], tri[1], tri[2]);
                                    const pBase = new THREE.Vector3().addScaledVector(triVecs[0], w.wA).addScaledVector(triVecs[1], w.wB).addScaledVector(triVecs[2], w.wC);
                                    const p0_u = proj.fromRot(p.u, p.v);
                                    const norm = proj.getNormal(p0_u.u, p0_u.v);
                                    const vFinal = pBase.addScaledVector(norm, p.h);
                                    if (proj.depthSign < 0) {
                                        // Reverse order for internal surfaces
                                        if (pIdx === 0) tempTri[0] = vFinal;
                                        else if (pIdx === 1) tempTri[1] = vFinal;
                                        else tempTri[2] = vFinal;
                                    } else {
                                        newPos.push(vFinal.x, vFinal.y, vFinal.z);
                                    }
                                });
                                if (proj.depthSign < 0) {
                                    newPos.push(tempTri[2].x, tempTri[2].y, tempTri[2].z, tempTri[1].x, tempTri[1].y, tempTri[1].z, tempTri[0].x, tempTri[0].y, tempTri[0].z);
                                }
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
    const repaired = repairMesh(finalGeo, params.holeFillThreshold, params.holeFillEnabled); repaired.computeVertexNormals(); return repaired;
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

        // Period along the circumference depends on the orientation
        const periodFactor = (bestSym % 60 === 0) ? 1.0 : Math.sqrt(3);

        const nUnits = Math.round(proj.circPhys / (W * periodFactor));
        W = proj.circPhys / (Math.max(1, nUnits) * periodFactor);
    }

    const appliedDepth = direction === 'inward' ? -depth : depth;
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
                            for (let v = 1; v < poly.length - 1; v++) {
                                const pts = [poly[0], poly[v], poly[v + 1]];
                                const tempTri: THREE.Vector3[] = [new THREE.Vector3(), new THREE.Vector3(), new THREE.Vector3()];
                                pts.forEach((p, pIdx) => {
                                    const w = getBarycentric(p, tri[0], tri[1], tri[2]);
                                    const pBase = new THREE.Vector3().addScaledVector(triVecs[0], w.wA).addScaledVector(triVecs[1], w.wB).addScaledVector(triVecs[2], w.wC);
                                    const p0_u = proj.fromRot(p.u, p.v);
                                    const norm = proj.getNormal(p0_u.u, p0_u.v);
                                    const vFinal = pBase.addScaledVector(norm, (p as any).h || 0);
                                    if (proj.depthSign < 0) {
                                        if (pIdx === 0) tempTri[0] = vFinal;
                                        else if (pIdx === 1) tempTri[1] = vFinal;
                                        else tempTri[2] = vFinal;
                                    } else {
                                        newPos.push(vFinal.x, vFinal.y, vFinal.z);
                                    }
                                });
                                if (proj.depthSign < 0) {
                                    newPos.push(tempTri[2].x, tempTri[2].y, tempTri[2].z, tempTri[1].x, tempTri[1].y, tempTri[1].z, tempTri[0].x, tempTri[0].y, tempTri[0].z);
                                }
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
    const repaired = repairMesh(finalGeo, params.holeFillThreshold, params.holeFillEnabled); repaired.computeVertexNormals(); return repaired;
}

// --------------------------------------------------------------------------------
// FUZZY SKIN IMPLEMENTATION
// --------------------------------------------------------------------------------

export function applyFuzzySkin(
    geometry: THREE.BufferGeometry,
    faceIndices: number[],
    params: FuzzySkinParams
): THREE.BufferGeometry {
    const islands = getContiguousIslands(geometry, faceIndices);
    const posAttr = geometry.attributes.position;
    const newPos: number[] = [];

    const distSq = (x1: number, y1: number, z1: number, x2: number, y2: number, z2: number) => {
        return (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2;
    };

    islands.forEach(islandIndices => {
        const islandPos = new Float32Array(islandIndices.length * 9);
        islandIndices.forEach((fIdx, i) => {
            const off = fIdx * 3;
            islandPos[i * 9 + 0] = posAttr.getX(off); islandPos[i * 9 + 1] = posAttr.getY(off); islandPos[i * 9 + 2] = posAttr.getZ(off);
            islandPos[i * 9 + 3] = posAttr.getX(off + 1); islandPos[i * 9 + 4] = posAttr.getY(off + 1); islandPos[i * 9 + 5] = posAttr.getZ(off + 1);
            islandPos[i * 9 + 6] = posAttr.getX(off + 2); islandPos[i * 9 + 7] = posAttr.getY(off + 2); islandPos[i * 9 + 8] = posAttr.getZ(off + 2);
        });

        // 1. Homogeneous Density Calculation
        let currentSubPos: Float32Array<ArrayBufferLike> = islandPos;
        const targetSq = params.pointDistance * params.pointDistance;
        const MAX_STEPS = 6;

        for (let s = 0; s < MAX_STEPS; s++) {
            let needsMore = false;
            const triCount = currentSubPos.length / 9;
            for (let i = 0; i < triCount; i++) {
                const o = i * 9;
                const d12 = distSq(currentSubPos[o], currentSubPos[o + 1], currentSubPos[o + 2], currentSubPos[o + 3], currentSubPos[o + 4], currentSubPos[o + 5]);
                const d23 = distSq(currentSubPos[o + 3], currentSubPos[o + 4], currentSubPos[o + 5], currentSubPos[o + 6], currentSubPos[o + 7], currentSubPos[o + 8]);
                const d31 = distSq(currentSubPos[o + 6], currentSubPos[o + 7], currentSubPos[o + 8], currentSubPos[o + 0], currentSubPos[o + 1], currentSubPos[o + 2]);
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

        // 3. Application of Noise (Random Displacement)
        for (let i = 0; i < positions.length / 3; i++) {
            const noise = (Math.random() - 0.5) * 2 * params.thickness;
            positions[i * 3] += normals[i * 3] * noise;
            positions[i * 3 + 1] += normals[i * 3 + 1] * noise;
            positions[i * 3 + 2] += normals[i * 3 + 2] * noise;
        }

        // 4. Conversion to final positions
        const finalIslandPos = subGeo.toNonIndexed().attributes.position.array as any as Float32Array;
        for (let i = 0; i < finalIslandPos.length; i++) {
            newPos.push(finalIslandPos[i]);
        }
    });

    const finalR: number[] = []; const sel = new Set(faceIndices);
    for (let f = 0; f < posAttr.count / 3; f++) { if (!sel.has(f)) { for (let v = 0; v < 3; v++) { const idx = f * 3 + v; finalR.push(posAttr.getX(idx), posAttr.getY(idx), posAttr.getZ(idx)); } } }
    const combined = new Float32Array(finalR.length + newPos.length); combined.set(finalR); combined.set(newPos, finalR.length);
    const finalGeo = new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(combined, 3));
    const repaired = repairMesh(finalGeo, params.holeFillThreshold, params.holeFillEnabled); repaired.computeVertexNormals(); return repaired;
}

export function applyDecimate(geo: THREE.BufferGeometry, _faceIndices: number[], _params: any) { return geo; }
