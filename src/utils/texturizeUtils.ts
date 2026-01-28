import * as THREE from 'three';

export type KnurlPattern = 'diamond' | 'straight' | 'diagonal' | 'square';

interface KnurlingParams {
    pitch: number;
    depth: number;
    angle: number;
    pattern: KnurlPattern;
}

export function applyKnurling(
    geometry: THREE.BufferGeometry,
    faceIndices: number[],
    params: KnurlingParams
): THREE.BufferGeometry {
    const { pitch: targetPitch, depth } = params;
    const posAttr = geometry.attributes.position;

    // 1. DATA COLLECTION
    const points: THREE.Vector3[] = [];
    const normals: THREE.Vector3[] = [];
    const centroids: THREE.Vector3[] = [];
    faceIndices.forEach(fIdx => {
        const off = fIdx * 3;
        const vA = new THREE.Vector3().fromBufferAttribute(posAttr, off + 0);
        const vB = new THREE.Vector3().fromBufferAttribute(posAttr, off + 1);
        const vC = new THREE.Vector3().fromBufferAttribute(posAttr, off + 2);
        const center = new THREE.Vector3().add(vA).add(vB).add(vC).divideScalar(3);
        const n = new THREE.Vector3().crossVectors(new THREE.Vector3().subVectors(vB, vA), new THREE.Vector3().subVectors(vC, vA)).normalize();
        points.push(vA, vB, vC); normals.push(n); centroids.push(center);
    });

    if (centroids.length === 0) return geometry;

    // ROBUST AXIS DETECTION (Centered Covariance PCA)
    // The axis of a revolution surface is the direction with minimum variance of normal projections
    const avgN = new THREE.Vector3(); normals.forEach(n => avgN.add(n)); avgN.divideScalar(normals.length);
    const avgN_norm = avgN.clone().normalize();
    let nVar = 0; normals.forEach(n => nVar += (1 - n.dot(avgN_norm))); nVar /= normals.length;
    const isPlanar = nVar < 0.01;

    let axis = new THREE.Vector3(1, 1, 1).normalize();
    if (!isPlanar) {
        let mxx = 0, mxy = 0, mxz = 0, myy = 0, myz = 0, mzz = 0;
        normals.forEach(n => {
            const dx = n.x - avgN.x, dy = n.y - avgN.y, dz = n.z - avgN.z;
            mxx += dx * dx; mxy += dx * dy; mxz += dx * dz; myy += dy * dy; myz += dy * dz; mzz += dz * dz;
        });
        const tr = mxx + myy + mzz;
        for (let i = 0; i < 20; i++) {
            const nx = (tr - mxx) * axis.x - mxy * axis.y - mxz * axis.z;
            const ny = -mxy * axis.x + (tr - myy) * axis.y - myz * axis.z;
            const nz = -mxz * axis.x - myz * axis.y + (tr - mzz) * axis.z;
            axis.set(nx, ny, nz).normalize();
        }
    }

    interface UVProjection {
        project: (p: THREE.Vector3) => { u: number, v: number, h: number };
        unproject: (u: number, v: number, h: number) => THREE.Vector3;
        bounds: { uMin: number, uMax: number, vMin: number, vMax: number };
        isFullCircle?: boolean;
    }

    let mapper: UVProjection;
    let rStart_global = 0, rEnd_global = 0, depthSign_global = 1;

    if (isPlanar) {
        const upP = avgN.clone().normalize(), rightP = new THREE.Vector3();
        if (Math.abs(upP.y) < 0.9) rightP.set(0, 1, 0).cross(upP).normalize(); else rightP.set(1, 0, 0).cross(upP).normalize();
        const fwdP = new THREE.Vector3().crossVectors(rightP, upP).normalize();
        const center = new THREE.Vector3(); centroids.forEach(c => center.add(c)); center.divideScalar(centroids.length);
        mapper = {
            project: (p) => { const v = p.clone().sub(center); return { u: v.dot(rightP), v: v.dot(fwdP), h: v.dot(upP) }; },
            unproject: (u, v, h) => rightP.clone().multiplyScalar(u).add(fwdP.clone().multiplyScalar(v)).add(upP.clone().multiplyScalar(h)).add(center),
            bounds: { uMin: 0, uMax: 0, vMin: 0, vMax: 0 }
        };
        let uMin = Infinity, uMax = -Infinity, vMin = Infinity, vMax = -Infinity;
        points.forEach(p => { const pr = mapper.project(p); uMin = Math.min(uMin, pr.u); uMax = Math.max(uMax, pr.u); vMin = Math.min(vMin, pr.v); vMax = Math.max(vMax, pr.v); });
        mapper.bounds = { uMin, uMax, vMin, vMax };
    } else {
        const up = axis, right = new THREE.Vector3();
        if (Math.abs(up.y) < 0.9) right.set(0, 1, 0).cross(up).normalize(); else right.set(1, 0, 0).cross(up).normalize();
        const fwd = new THREE.Vector3().crossVectors(up, right).normalize();

        // 2D Origin Solver (Intersection of projected normal lines)
        let sUU = 0, sUV = 0, sVV = 0, bU = 0, bV = 0;
        centroids.forEach((q, idx) => {
            const n = normals[idx].clone().projectOnPlane(up).normalize();
            if (n.lengthSq() < 0.001) return;
            const qu = q.dot(right), qv = q.dot(fwd), nu = n.dot(right), nv = n.dot(fwd);
            const a = 1 - nu * nu, b = -nu * nv, c = 1 - nv * nv;
            sUU += a; sUV += b; sVV += c; bU += a * qu + b * qv; bV += b * qu + c * qv;
        });
        const det = sUU * sVV - sUV * sUV;
        const axisOrigin = (Math.abs(det) > 1e-8)
            ? right.clone().multiplyScalar((sVV * bU - sUV * bV) / det).add(fwd.clone().multiplyScalar((sUU * bV - sUV * bU) / det))
            : centroids[0].clone();

        const pCyl = points.map(p => {
            const v = p.clone().sub(axisOrigin);
            return { h: v.dot(up), a: Math.atan2(v.dot(fwd), v.dot(right)), r: v.clone().projectOnPlane(up).length() };
        });

        // Slanted profile analysis (r = m*h + b)
        let sumH = 0, sumR = 0, sumHH = 0, sumHR = 0, count = pCyl.length;
        pCyl.forEach(p => { sumH += p.h; sumR += p.r; sumHH += p.h * p.h; sumHR += p.h * p.r; });
        const den = (count * sumHH - sumH * sumH);
        const m = Math.abs(den) > 1e-9 ? (count * sumHR - sumH * sumR) / den : 0;
        const b = (sumR - m * sumH) / count;

        let hMin = Infinity, hMax = -Infinity, aMin = Infinity, aMax = -Infinity;
        pCyl.forEach(p => { hMin = Math.min(hMin, p.h); hMax = Math.max(hMax, p.h); aMin = Math.min(aMin, p.a); aMax = Math.max(aMax, p.a); });
        const isFull = (aMax - aMin) > 5.5; if (isFull) { aMin = -Math.PI; aMax = Math.PI; }
        const hR = hMax - hMin;
        rStart_global = m * hMin + b; rEnd_global = m * hMax + b;

        // Unified Normal Analysis
        let nrSum = 0, nvSum = 0;
        centroids.forEach((c, i) => {
            const toS = c.clone().sub(axisOrigin).projectOnPlane(up).normalize();
            nrSum += normals[i].dot(toS); nvSum += normals[i].dot(up);
        });
        const nLen = Math.sqrt(nrSum * nrSum + nvSum * nvSum);
        const unitNR = nrSum / (nLen || 1), unitNV = nvSum / (nLen || 1);
        depthSign_global = unitNR < 0 ? -1 : 1;

        mapper = {
            project: (p) => { const v = p.clone().sub(axisOrigin); return { u: Math.atan2(v.dot(fwd), v.dot(right)), v: v.dot(up), h: 0 }; },
            unproject: (u, v, h) => {
                const r = (m * v + b) + h * unitNR, hO = h * unitNV;
                return right.clone().multiplyScalar(Math.cos(u) * r).add(fwd.clone().multiplyScalar(Math.sin(u) * r)).add(up.clone().multiplyScalar(v + hO)).add(axisOrigin);
            },
            bounds: { uMin: aMin, uMax: aMax, vMin: hMin, vMax: hMax }, isFullCircle: isFull
        };
    }

    // 2. TILING
    const uRange = mapper.bounds.uMax - mapper.bounds.uMin, vRange = mapper.bounds.vMax - mapper.bounds.vMin;
    let physULen = isPlanar ? uRange : uRange * ((rStart_global + rEnd_global) / 2);
    let physVLen = isPlanar ? vRange : Math.sqrt(vRange * vRange + (rEnd_global - rStart_global) ** 2);
    let nD = Math.max(2, Math.min(1500, targetPitch > 1e-6 ? Math.round(physULen / targetPitch) : 2));
    if (!isPlanar && mapper.isFullCircle && nD % 2 !== 0) nD++;
    const mD = Math.max(1, Math.min(1500, targetPitch > 1e-6 ? Math.round(physVLen / targetPitch) : 1)), pU = uRange / nD, pV = vRange / mD;

    // 3. GEOMETRY ENGINE
    const newPos: number[] = [];
    if (nD > 2000 || mD > 2000 || isNaN(nD) || isNaN(mD)) return geometry; // Extra safety bail
    const origTriangles: { u: number, v: number }[][] = [];
    for (let t = 0; t < points.length; t += 3) {
        let pA = mapper.project(points[t]), pB = mapper.project(points[t + 1]), pC = mapper.project(points[t + 2]);
        if (!isPlanar && mapper.isFullCircle && Math.max(pA.u, pB.u, pC.u) - Math.min(pA.u, pB.u, pC.u) > Math.PI) {
            if (pA.u < 0) pA.u += 2 * Math.PI; if (pB.u < 0) pB.u += 2 * Math.PI; if (pC.u < 0) pC.u += 2 * Math.PI;
        }
        const area = (pB.u - pA.u) * (pC.v - pA.v) - (pB.v - pA.v) * (pC.u - pA.u);
        origTriangles.push(area < 0 ? [pA, pC, pB] : [pA, pB, pC]);
    }

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

    origTriangles.forEach(tri => {
        const iS = Math.floor((Math.min(tri[0].u, tri[1].u, tri[2].u) - mapper.bounds.uMin) / pU), iE = Math.ceil((Math.max(tri[0].u, tri[1].u, tri[2].u) - mapper.bounds.uMin) / pU);
        const jS = Math.floor((Math.min(tri[0].v, tri[1].v, tri[2].v) - mapper.bounds.vMin) / pV), jE = Math.ceil((Math.max(tri[0].v, tri[1].v, tri[2].v) - mapper.bounds.vMin) / pV);
        for (let j = jS; j <= jE; j++) {
            for (let i = iS; i <= iE; i++) {
                const u0 = mapper.bounds.uMin + i * pU, u1 = mapper.bounds.uMin + (i + 1) * pU, v0 = mapper.bounds.vMin + j * pV, v1 = mapper.bounds.vMin + (j + 1) * pV;
                const mid = { u: u0 + pU * 0.5, v: v0 + pV * 0.5, h: ((Math.abs(i) + Math.abs(j)) % 2 === 0) ? depth : 0 };
                const cell = [[{ u: u0, v: v0, h: 0 }, { u: u1, v: v0, h: 0 }, mid], [{ u: u1, v: v0, h: 0 }, { u: u1, v: v1, h: 0 }, mid], [{ u: u1, v: v1, h: 0 }, { u: u0, v: v1, h: 0 }, mid], [{ u: u0, v: v1, h: 0 }, { u: u0, v: v0, h: 0 }, mid]];
                cell.forEach(kTri => {
                    let poly = kTri; poly = clip(poly, tri[0], tri[1]); poly = clip(poly, tri[1], tri[2]); poly = clip(poly, tri[2], tri[0]);
                    if (poly.length >= 3) {
                        for (let v = 1; v < poly.length - 1; v++) {
                            const p0 = mapper.unproject(poly[0].u, poly[0].v, poly[0].h), p1 = mapper.unproject(poly[v].u, poly[v].v, poly[v].h), p2 = mapper.unproject(poly[v + 1].u, poly[v + 1].v, poly[v + 1].h);
                            if (depthSign_global < 0) newPos.push(p0.x, p0.y, p0.z, p2.x, p2.y, p2.z, p1.x, p1.y, p1.z);
                            else newPos.push(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
                        }
                    }
                });
            }
        }
    });

    const finalR: number[] = []; const sel = new Set(faceIndices);
    for (let f = 0; f < posAttr.count / 3; f++) { if (!sel.has(f)) { for (let v = 0; v < 3; v++) { const idx = f * 3 + v; finalR.push(posAttr.getX(idx), posAttr.getY(idx), posAttr.getZ(idx)); } } }
    const combined = new Float32Array(finalR.length + newPos.length); combined.set(finalR); combined.set(newPos, finalR.length);
    const finalGeo = new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(combined, 3));
    const repaired = repairMesh(finalGeo); repaired.computeVertexNormals(); return repaired;
}

function repairMesh(geo: THREE.BufferGeometry): THREE.BufferGeometry {
    let welded = weldGeometry(geo, 100000);
    const pos = welded.attributes.position, vertexCount = pos.count, edgeCounts = new Map<string, number>();
    for (let i = 0; i < vertexCount; i += 3) {
        for (let j = 0; j < 3; j++) {
            const v1 = i + j, v2 = i + (j + 1) % 3, h1 = getVKey(pos.getX(v1), pos.getY(v1), pos.getZ(v1)), h2 = getVKey(pos.getX(v2), pos.getY(v2), pos.getZ(v2));
            const edgeKey = h1 < h2 ? `${h1}|${h2}` : `${h2}|${h1}`; edgeCounts.set(edgeKey, (edgeCounts.get(edgeKey) || 0) + 1);
        }
    }
    const bIdx = new Set<number>();
    for (let i = 0; i < vertexCount; i += 3) {
        for (let j = 0; j < 3; j++) {
            const v1 = i + j, v2 = i + (j + 1) % 3, h1 = getVKey(pos.getX(v1), pos.getY(v1), pos.getZ(v1)), h2 = getVKey(pos.getX(v2), pos.getY(v2), pos.getZ(v2));
            const edgeKey = h1 < h2 ? `${h1}|${h2}` : `${h2}|${h1}`; if (edgeCounts.get(edgeKey) === 1) { bIdx.add(v1); bIdx.add(v2); }
        }
    }
    if (bIdx.size > 0) {
        const snapT = 0.05, positions = welded.attributes.position.array as Float32Array, grid = new Map<string, number>();
        bIdx.forEach(idx => {
            const x = positions[idx * 3], y = positions[idx * 3 + 1], z = positions[idx * 3 + 2], gx = Math.floor(x / snapT), gy = Math.floor(y / snapT), gz = Math.floor(z / snapT);
            let snapped = false;
            for (let dx = -1; dx <= 1 && !snapped; dx++) {
                for (let dy = -1; dy <= 1 && !snapped; dy++) {
                    for (let dz = -1; dz <= 1 && !snapped; dz++) {
                        const key = `${gx + dx},${gy + dy},${gz + dz}`;
                        if (grid.has(key)) {
                            const tIdx = grid.get(key)!; const tx = positions[tIdx * 3], ty = positions[tIdx * 3 + 1], tz = positions[tIdx * 3 + 2];
                            if ((x - tx) ** 2 + (y - ty) ** 2 + (z - tz) ** 2 < snapT ** 2) { positions[idx * 3] = tx; positions[idx * 3 + 1] = ty; positions[idx * 3 + 2] = tz; snapped = true; }
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

function getVKey(x: number, y: number, z: number): string { return `${Math.round(x * 100000)},${Math.round(y * 100000)},${Math.round(z * 100000)}`; }

function weldGeometry(geo: THREE.BufferGeometry, prec: number = 100000): THREE.BufferGeometry {
    const pos = geo.attributes.position, vMap = new Map<string, number>(), nPos: number[] = [], indices: number[] = [];
    for (let i = 0; i < pos.count; i++) {
        const x = pos.getX(i), y = pos.getY(i), z = pos.getZ(i), h = `${Math.round(x * prec)},${Math.round(y * prec)},${Math.round(z * prec)}`;
        if (!vMap.has(h)) { vMap.set(h, nPos.length / 3); nPos.push(x, y, z); }
        indices.push(vMap.get(h)!);
    }
    const fIdx: number[] = [];
    for (let f = 0; f < indices.length / 3; f++) {
        const i1 = indices[f * 3], i2 = indices[f * 3 + 1], i3 = indices[f * 3 + 2]; if (i1 !== i2 && i2 !== i3 && i3 !== i1) fIdx.push(i1, i2, i3);
    }
    const res = new THREE.BufferGeometry(); res.setAttribute('position', new THREE.BufferAttribute(new Float32Array(nPos), 3)); res.setIndex(fIdx); return res.toNonIndexed();
}

export function applyHoneycomb(geo: THREE.BufferGeometry, _faceIndices: number[], _params: any) { return geo; }
export function applyDecimate(geo: THREE.BufferGeometry, _faceIndices: number[], _params: any) { return geo; }
