/** Reusable incremental Delaunay output. Only the first tetrahedronCount entries are valid. */
export class DelaunayResult {
  tetrahedra = new Int32Array(0);
  circumcenters = new Float64Array(0);
  tetrahedronCount = 0;

  reserve(count: number): void {
    if (this.tetrahedra.length >= count * 4) return;
    let size = Math.max(16, this.tetrahedra.length / 4);
    while (size < count) size *= 2;
    const tets = new Int32Array(size * 4);
    const centers = new Float64Array(size * 3);
    tets.set(this.tetrahedra);
    centers.set(this.circumcenters);
    this.tetrahedra = tets;
    this.circumcenters = centers;
  }

  removeSwapBack(index: number): void {
    const last = --this.tetrahedronCount;
    const dst = index * 4, src = last * 4;
    for (let i = 0; i < 4; i++) this.tetrahedra[dst + i] = this.tetrahedra[src + i];
    const cdst = index * 3, csrc = last * 3;
    for (let i = 0; i < 3; i++) this.circumcenters[cdst + i] = this.circumcenters[csrc + i];
  }
}

/** Reusable workspace. Use a separate instance for concurrent builds. */
export class DelaunayScratch {
  points = new Float64Array(0);
  faceA = new Int32Array(0);
  faceB = new Int32Array(0);
  faceC = new Int32Array(0);
  faceActive = new Uint8Array(0);
  faceCount = 0;

  hashA = new Int32Array(0);
  hashB = new Int32Array(0);
  hashC = new Int32Array(0);
  hashFace = new Int32Array(0);
  hashStamp = new Uint32Array(0);
  generation = 0;

  reservePoints(siteCount: number): void {
    const needed = (siteCount + 4) * 3;
    if (this.points.length < needed) this.points = new Float64Array(needed);
  }

  resetFaces(maxFaces: number): void {
    this.faceCount = 0;
    if (this.faceA.length < maxFaces) {
      this.faceA = new Int32Array(maxFaces);
      this.faceB = new Int32Array(maxFaces);
      this.faceC = new Int32Array(maxFaces);
      this.faceActive = new Uint8Array(maxFaces);
    }
    let hashSize = 16;
    while (hashSize < maxFaces * 2) hashSize *= 2;
    if (this.hashA.length < hashSize) {
      this.hashA = new Int32Array(hashSize);
      this.hashB = new Int32Array(hashSize);
      this.hashC = new Int32Array(hashSize);
      this.hashFace = new Int32Array(hashSize);
      this.hashStamp = new Uint32Array(hashSize);
      this.generation = 0;
    }
    this.generation = (this.generation + 1) >>> 0;
    if (this.generation === 0) {
      this.hashStamp.fill(0);
      this.generation = 1;
    }
  }

  toggleFace(a: number, b: number, c: number): void {
    if (a > b) { const swap = a; a = b; b = swap; }
    if (b > c) { const swap = b; b = c; c = swap; }
    if (a > b) { const swap = a; a = b; b = swap; }
    const mask = this.hashA.length - 1;
    let slot = (Math.imul(a, 0x9e3779b1) ^ Math.imul(b, 0x85ebca6b) ^
                Math.imul(c, 0xc2b2ae35)) & mask;
    while (true) {
      if (this.hashStamp[slot] !== this.generation) {
        const face = this.faceCount++;
        this.hashStamp[slot] = this.generation;
        this.hashA[slot] = a; this.hashB[slot] = b; this.hashC[slot] = c;
        this.hashFace[slot] = face;
        this.faceA[face] = a; this.faceB[face] = b; this.faceC[face] = c;
        this.faceActive[face] = 1;
        return;
      }
      if (this.hashA[slot] === a && this.hashB[slot] === b && this.hashC[slot] === c) {
        this.faceActive[this.hashFace[slot]] = 0;
        return;
      }
      slot = (slot + 1) & mask;
    }
  }
}

// Solve three perpendicular-bisector equations in double precision.
function writeCircumcenter(points: Float64Array, a: number, b: number, c: number, d: number,
                           out: Float64Array, index: number): boolean {
  const p = a * 3, q = b * 3, r = c * 3, s = d * 3;
  const ux = points[q] - points[p], uy = points[q+1] - points[p+1], uz = points[q+2] - points[p+2];
  const vx = points[r] - points[p], vy = points[r+1] - points[p+1], vz = points[r+2] - points[p+2];
  const wx = points[s] - points[p], wy = points[s+1] - points[p+1], wz = points[s+2] - points[p+2];
  const vwx = vy*wz-vz*wy, vwy = vz*wx-vx*wz, vwz = vx*wy-vy*wx;
  const wux = wy*uz-wz*uy, wuy = wz*ux-wx*uz, wuz = wx*uy-wy*ux;
  const uvx = uy*vz-uz*vy, uvy = uz*vx-ux*vz, uvz = ux*vy-uy*vx;
  const det2 = 2 * (ux*vwx + uy*vwy + uz*vwz);
  const scale = Math.max(Math.hypot(ux,uy,uz), Math.hypot(vx,vy,vz), Math.hypot(wx,wy,wz));
  if (Math.abs(det2) <= 1e-12 * scale * scale * scale) return false;
  const u2 = ux*ux+uy*uy+uz*uz, v2 = vx*vx+vy*vy+vz*vz, w2 = wx*wx+wy*wy+wz*wz;
  const dest = index * 3;
  out[dest] = points[p] + (u2*vwx+v2*wux+w2*uvx)/det2;
  out[dest+1] = points[p+1] + (u2*vwy+v2*wuy+w2*uvy)/det2;
  out[dest+2] = points[p+2] + (u2*vwz+v2*wuz+w2*uvz)/det2;
  return Number.isFinite(out[dest]) && Number.isFinite(out[dest+1]) && Number.isFinite(out[dest+2]);
}

/**
 * Incremental Bowyer-Watson 3D Delaunay tetrahedralization. Sites must be
 * finite, distinct, and in general position. Reuse workspace/output to avoid
 * steady-state heap allocations; capacity growth may allocate.
 */
export function buildDelaunay3D(sites: ArrayLike<number>, scratch: DelaunayScratch,
                                output: DelaunayResult): DelaunayResult {
  if (sites.length % 3 || sites.length < 12) throw new RangeError("At least four 3D sites are required");
  const n = sites.length / 3;
  scratch.reservePoints(n);
  const points = scratch.points;
  let minX = sites[0], minY = sites[1], minZ = sites[2];
  let maxX = minX, maxY = minY, maxZ = minZ;
  for (let i = 0; i < sites.length; i += 3) {
    const x = sites[i], y = sites[i+1], z = sites[i+2];
    if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z))
      throw new RangeError("Sites must be finite");
    points[i] = x; points[i+1] = y; points[i+2] = z;
    minX = Math.min(minX,x); minY = Math.min(minY,y); minZ = Math.min(minZ,z);
    maxX = Math.max(maxX,x); maxY = Math.max(maxY,y); maxZ = Math.max(maxZ,z);
  }
  const midX = (minX+maxX)/2, midY = (minY+maxY)/2, midZ = (minZ+maxZ)/2;
  const radius = Math.max(maxX-minX, maxY-minY, maxZ-minZ) * 16;
  if (!(radius > 0) || !Number.isFinite(radius)) throw new RangeError("Site bounds must be finite and nonzero");
  let sp = n*3;
  points[sp++] = midX+radius; points[sp++] = midY+radius; points[sp++] = midZ+radius;
  points[sp++] = midX-radius; points[sp++] = midY-radius; points[sp++] = midZ+radius;
  points[sp++] = midX-radius; points[sp++] = midY+radius; points[sp++] = midZ-radius;
  points[sp++] = midX+radius; points[sp++] = midY-radius; points[sp++] = midZ-radius;

  output.tetrahedronCount = 0;
  output.reserve(1);
  output.tetrahedra[0] = n; output.tetrahedra[1] = n+1;
  output.tetrahedra[2] = n+2; output.tetrahedra[3] = n+3;
  if (!writeCircumcenter(points,n,n+1,n+2,n+3,output.circumcenters,0))
    throw new Error("Could not initialize the enclosing tetrahedron");
  output.tetrahedronCount = 1;

  for (let point = 0; point < n; point++) {
    scratch.resetFaces(output.tetrahedronCount * 4);
    let badCount = 0;
    for (let t = 0; t < output.tetrahedronCount;) {
      const ti = t*4, ci = t*3;
      const a = output.tetrahedra[ti], b = output.tetrahedra[ti+1],
            c = output.tetrahedra[ti+2], d = output.tetrahedra[ti+3];
      const x = output.circumcenters[ci], y = output.circumcenters[ci+1], z = output.circumcenters[ci+2];
      const ai = a*3, pi = point*3;
      const ax = x-points[ai], ay = y-points[ai+1], az = z-points[ai+2];
      const px = x-points[pi], py = y-points[pi+1], pz = z-points[pi+2];
      const radiusSq = ax*ax+ay*ay+az*az;
      if (px*px+py*py+pz*pz <= radiusSq*(1+1e-6)) {
        scratch.toggleFace(a,b,c); scratch.toggleFace(a,b,d);
        scratch.toggleFace(a,c,d); scratch.toggleFace(b,c,d);
        output.removeSwapBack(t);
        badCount++;
      } else t++;
    }
    if (!badCount) throw new Error("No Delaunay cavity found; check site precision and degeneracy");
    for (let f = 0; f < scratch.faceCount; f++) {
      if (!scratch.faceActive[f]) continue;
      output.reserve(output.tetrahedronCount+1);
      const ti = output.tetrahedronCount * 4;
      const a = scratch.faceA[f], b = scratch.faceB[f], c = scratch.faceC[f];
      output.tetrahedra[ti] = a; output.tetrahedra[ti+1] = b;
      output.tetrahedra[ti+2] = c; output.tetrahedra[ti+3] = point;
      if (!writeCircumcenter(points,a,b,c,point,output.circumcenters,output.tetrahedronCount))
        throw new Error("Degenerate tetrahedron; sites must be in general position");
      output.tetrahedronCount++;
    }
  }
  for (let t = 0; t < output.tetrahedronCount;) {
    const i = t*4;
    if (output.tetrahedra[i] >= n || output.tetrahedra[i+1] >= n ||
        output.tetrahedra[i+2] >= n || output.tetrahedra[i+3] >= n)
      output.removeSwapBack(t);
    else t++;
  }
  return output;
}
