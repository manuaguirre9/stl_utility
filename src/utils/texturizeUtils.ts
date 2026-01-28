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
    const { pitch: targetPitch, depth, angle: userAngle } = params;
    const posAttr = geometry.attributes.position;

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

    const avgN = new THREE.Vector3(); normals.forEach(n => avgN.add(n)); avgN.divideScalar(normals.length);
    const avgN_norm = avgN.clone().normalize();
    let nVar = 0; normals.forEach(n => nVar += (1 - n.dot(avgN_norm))); nVar /= normals.length;
    const isPlanar = nVar < 0.01;

    let axis = new THREE.Vector3(1, 1, 1).normalize();
    let axisOrigin = new THREE.Vector3();
    let rStart_global = 0, rEnd_global = 0, depthSign_global = 1, m = 0, b = 0;

    if (!isPlanar) {
        let mxx = 0, mxy = 0, mxz = 0, myy = 0, myz = 0, mzz = 0;
        normals.forEach(n => {
            const dx = n.x - avgN.x, dy = n.y - avgN.y, dz = n.z - avgN.z;
            mxx += dx * dx; mxy += dx * dy; mxz += dx * dz; myy += dy * dy; myz += dy * dz; mzz += dz * dz;
        });
        const trace = mxx + myy + mzz;
        for (let i = 0; i < 20; i++) {
            const nx = (trace - mxx) * axis.x - mxy * axis.y - mxz * axis.z;
            const ny = -mxy * axis.x + (trace - myy) * axis.y - myz * axis.z;
            const nz = -mxz * axis.x - myz * axis.y + (trace - mzz) * axis.z;
            axis.set(nx, ny, nz).normalize();
        }
        const up = axis, right = new THREE.Vector3();
        if (Math.abs(up.y) < 0.9) right.set(0, 1, 0).cross(up).normalize(); else right.set(1, 0, 0).cross(up).normalize();
        const fwd = new THREE.Vector3().crossVectors(up, right).normalize();

        let sUU = 0, sUV = 0, sVV = 0, bU = 0, bV = 0;
        centroids.forEach((q, idx) => {
            const n = normals[idx].clone().projectOnPlane(up).normalize();
            if (n.lengthSq() < 0.001) return;
            const qu = q.dot(right), qv = q.dot(fwd), nu = n.dot(right), nv = n.dot(fwd);
            const aa = 1 - nu * nu, ab = -nu * nv, ac = 1 - nv * nv;
            sUU += aa; sUV += ab; sVV += ac; bU += aa * qu + ab * qv; bV += ab * qu + ac * qv;
        });
        const det = sUU * sVV - sUV * sUV;
        axisOrigin = (Math.abs(det) > 1e-8) ? right.clone().multiplyScalar((sVV * bU - sUV * bV) / det).add(fwd.clone().multiplyScalar((sUU * bV - sUV * bU) / det)) : centroids[0].clone();

        const pCyl = points.map(p => { const v = p.clone().sub(axisOrigin); return { h: v.dot(up), r: v.clone().projectOnPlane(up).length() }; });
        let sumH = 0, sumR = 0, sumHH = 0, sumHR = 0, cnt = pCyl.length;
        pCyl.forEach(p => { sumH += p.h; sumR += p.r; sumHH += p.h * p.h; sumHR += p.h * p.r; });
        const den = (cnt * sumHH - sumH * sumH);
        m = Math.abs(den) > 1e-9 ? (cnt * sumHR - sumH * sumR) / den : 0; b = (sumR - m * sumH) / cnt;

        let nrSum = 0, nvSum = 0;
        centroids.forEach((c, i) => { const toS = c.clone().sub(axisOrigin).projectOnPlane(up).normalize(); nrSum += normals[i].dot(toS); nvSum += normals[i].dot(up); });
        const nLen = Math.sqrt(nrSum * nrSum + nvSum * nvSum); depthSign_global = (nrSum / (nLen || 1)) < 0 ? -1 : 1;
        const hPoints = pCyl.map(p => p.h); const hMin = Math.min(...hPoints), hMax = Math.max(...hPoints);
        rStart_global = m * hMin + b; rEnd_global = m * hMax + b;
    }

    const upDir = isPlanar ? avgN_norm : axis;
    const rightDir = new THREE.Vector3();
    if (Math.abs(upDir.y) < 0.9) rightDir.set(0, 1, 0).cross(upDir).normalize(); else rightDir.set(1, 0, 0).cross(upDir).normalize();
    const fwdDir = new THREE.Vector3().crossVectors(upDir, rightDir).normalize();
    const origin = isPlanar ? centroids.reduce((acc, c) => acc.add(c), new THREE.Vector3()).divideScalar(centroids.length) : axisOrigin;

    // SEAMLESS UV LOGIC
    const avgR = isPlanar ? 1 : (rStart_global + rEnd_global) / 2;
    const circPhys = 2 * Math.PI * avgR;

    const project = (p: THREE.Vector3) => {
        const v = p.clone().sub(origin);
        if (isPlanar) return { u: v.dot(rightDir), v: v.dot(fwdDir) };
        let uRad = Math.atan2(v.dot(fwdDir), v.dot(rightDir));
        if (uRad < 0) uRad += 2 * Math.PI;
        return { u: uRad * avgR, v: v.dot(upDir) };
    };

    const unproject = (u: number, v: number, h: number) => {
        const r = isPlanar ? 1 : (m * v + b) + h * depthSign_global;
        if (isPlanar) return rightDir.clone().multiplyScalar(u).add(fwdDir.clone().multiplyScalar(v)).add(upDir.clone().multiplyScalar(h)).add(origin);
        const ang = u / avgR;
        const toS = rightDir.clone().multiplyScalar(Math.cos(ang)).add(fwdDir.clone().multiplyScalar(Math.sin(ang)));
        return toS.multiplyScalar(r).add(upDir.clone().multiplyScalar(v)).add(origin);
    };

    // ROTATION & TILING
    const angRad = (userAngle * Math.PI) / 180;
    const cosA = Math.cos(angRad), sinA = Math.sin(angRad);
    const toRot = (u: number, v: number) => ({ ur: u * cosA - v * sinA, vr: u * sinA + v * cosA });
    const fromRot = (ur: number, vr: number) => ({ u: ur * cosA + vr * sinA, v: -ur * sinA + vr * cosA });

    // Step adjustment for seamless wrap
    const nD_nom = Math.round(circPhys / targetPitch);
    const nD = Math.max(2, nD_nom % 2 === 0 ? nD_nom : nD_nom + 1);
    const pitchP = circPhys / nD;
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

    const newPos: number[] = [];
    for (let t = 0; t < points.length; t += 3) {
        let pA_uv = project(points[t]), pB_uv = project(points[t + 1]), pC_uv = project(points[t + 2]);

        // SEAM WRAP CORRECTION
        if (!isPlanar) {
            const uVals = [pA_uv.u, pB_uv.u, pC_uv.u];
            const minU = Math.min(...uVals), maxU = Math.max(...uVals);
            if (maxU - minU > circPhys * 0.5) {
                if (pA_uv.u < circPhys * 0.5) pA_uv.u += circPhys;
                if (pB_uv.u < circPhys * 0.5) pB_uv.u += circPhys;
                if (pC_uv.u < circPhys * 0.5) pC_uv.u += circPhys;
            }
        }

        const pA = toRot(pA_uv.u, pA_uv.v), pB = toRot(pB_uv.u, pB_uv.v), pC = toRot(pC_uv.u, pC_uv.v);
        const area = (pB.ur - pA.ur) * (pC.vr - pA.vr) - (pB.vr - pA.vr) * (pC.ur - pA.ur);
        const tri = area < 0 ? [{ u: pA.ur, v: pA.vr }, { u: pC.ur, v: pC.vr }, { u: pB.ur, v: pB.vr }] : [{ u: pA.ur, v: pA.vr }, { u: pB.ur, v: pB.vr }, { u: pC.ur, v: pC.vr }];

        const iS = Math.floor((Math.min(tri[0].u, tri[1].u, tri[2].u)) / pU) - 1, iE = Math.ceil((Math.max(tri[0].u, tri[1].u, tri[2].u)) / pU) + 1;
        const jS = Math.floor((Math.min(tri[0].v, tri[1].v, tri[2].v)) / pV) - 1, jE = Math.ceil((Math.max(tri[0].v, tri[1].v, tri[2].v)) / pV) + 1;

        for (let j = jS; j <= jE; j++) {
            for (let i = iS; i <= iE; i++) {
                const u0 = i * pU, u1 = (i + 1) * pU, v0 = j * pV, v1 = (j + 1) * pV;
                const mid = { u: u0 + pU * 0.5, v: v0 + pV * 0.5, h: ((Math.abs(i) + Math.abs(j)) % 2 === 0) ? depth : 0 };
                const cell = [[{ u: u0, v: v0, h: 0 }, { u: u1, v: v0, h: 0 }, mid], [{ u: u1, v: v0, h: 0 }, { u: u1, v: v1, h: 0 }, mid], [{ u: u1, v: v1, h: 0 }, { u: u0, v: v1, h: 0 }, mid], [{ u: u0, v: v1, h: 0 }, { u: u0, v: v0, h: 0 }, mid]];
                cell.forEach(kTri => {
                    let poly = kTri; poly = clip(poly, tri[0], tri[1]); poly = clip(poly, tri[1], tri[2]); poly = clip(poly, tri[2], tri[0]);
                    if (poly.length >= 3) {
                        for (let v = 1; v < poly.length - 1; v++) {
                            const p0_u = fromRot(poly[0].u, poly[0].v), p1_u = fromRot(poly[v].u, poly[v].v), p2_u = fromRot(poly[v + 1].u, poly[v + 1].v);
                            const vA = unproject(p0_u.u, p0_u.v, poly[0].h || 0), vB = unproject(p1_u.u, p1_u.v, poly[v].h || 0), vC = unproject(p2_u.u, p2_u.v, poly[v + 1].h || 0);
                            newPos.push(vA.x, vA.y, vA.z, vB.x, vB.y, vB.z, vC.x, vC.y, vC.z);
                        }
                    }
                });
            }
        }
    }

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
