/** Flat xyz positions, flat abcd tetrahedra, and one xyz circumcenter per tetrahedron. */
export interface VoronoiInput {
  sites: ArrayLike<number>;
  tetrahedra: ArrayLike<number>;
  circumcenters: ArrayLike<number>;
  filled: ArrayLike<number>; // 0 = empty, nonzero = filled
  emitCell?: ArrayLike<number>; // optional ownership mask for chunked output
  duplicateEpsilon?: number;
}

/** Reused output. Only positions/normals[0..vertexCount*3) and indices[0..indexCount) are valid. */
export class MeshBuffer {
  positions = new Float32Array(0);
  normals = new Float32Array(0);
  indices = new Uint32Array(0);
  vertexCount = 0;
  indexCount = 0;

  reserve(vertices: number, indices: number): void {
    const floats = vertices * 3;
    if (this.positions.length < floats) {
      const size = grow(this.positions.length, floats);
      const nextPositions = new Float32Array(size);
      const nextNormals = new Float32Array(size);
      nextPositions.set(this.positions);
      nextNormals.set(this.normals);
      this.positions = nextPositions;
      this.normals = nextNormals;
    }
    if (this.indices.length < indices) {
      const next = new Uint32Array(grow(this.indices.length, indices));
      next.set(this.indices);
      this.indices = next;
    }
  }
}

/** Reusable topology and polygon workspace. Use a separate instance for concurrent builds. */
export class Scratch {
  edgeA = new Int32Array(0);
  edgeB = new Int32Array(0);
  head = new Int32Array(0);
  next = new Int32Array(0);
  tetIndex = new Int32Array(0);
  faceCount = 0;
  nodeCount = 0;

  hashA = new Int32Array(0);
  hashB = new Int32Array(0);
  hashFace = new Int32Array(0);
  hashStamp = new Uint32Array(0);
  generation = 0;

  points = new Float64Array(0);
  angles = new Float64Array(0);

  reset(maxEdges: number): void {
    this.faceCount = this.nodeCount = 0;
    if (this.edgeA.length < maxEdges) {
      this.edgeA = new Int32Array(maxEdges);
      this.edgeB = new Int32Array(maxEdges);
      this.head = new Int32Array(maxEdges);
      this.next = new Int32Array(maxEdges);
      this.tetIndex = new Int32Array(maxEdges);
    }
    // Keep load under 50%. Generation stamps avoid clearing the hash table per run.
    const hashSize = nextPowerOfTwo(Math.max(16, maxEdges * 2));
    if (this.hashA.length < hashSize) {
      this.hashA = new Int32Array(hashSize);
      this.hashB = new Int32Array(hashSize);
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

  reservePoints(count: number): void {
    if (this.angles.length >= count) return;
    const size = grow(this.angles.length, count);
    const nextAngles = new Float64Array(size);
    const nextPoints = new Float64Array(size * 3);
    nextAngles.set(this.angles);
    nextPoints.set(this.points);
    this.angles = nextAngles;
    this.points = nextPoints;
  }
}

function grow(oldSize: number, needed: number): number {
  let size = Math.max(16, oldSize);
  while (size < needed) size *= 2;
  return size;
}

function nextPowerOfTwo(value: number): number {
  let size = 1;
  while (size < value) size *= 2;
  return size;
}

function addEdge(s: Scratch, a: number, b: number, tet: number): void {
  if (a === b) return;
  if (a > b) { const swap = a; a = b; b = swap; }
  const mask = s.hashA.length - 1;
  let slot = (Math.imul(a, 0x9e3779b1) ^ Math.imul(b, 0x85ebca6b)) & mask;
  let face: number;
  while (true) {
    if (s.hashStamp[slot] !== s.generation) {
      face = s.faceCount++;
      s.hashStamp[slot] = s.generation;
      s.hashA[slot] = a;
      s.hashB[slot] = b;
      s.hashFace[slot] = face;
      s.edgeA[face] = a;
      s.edgeB[face] = b;
      s.head[face] = -1;
      break;
    }
    if (s.hashA[slot] === a && s.hashB[slot] === b) {
      face = s.hashFace[slot];
      break;
    }
    slot = (slot + 1) & mask;
  }
  const node = s.nodeCount++;
  s.tetIndex[node] = tet;
  s.next[node] = s.head[face];
  s.head[face] = node;
}

function countThird(s: Scratch, face: number, tets: ArrayLike<number>, third: number): number {
  let count = 0;
  for (let node = s.head[face]; node >= 0; node = s.next[node]) {
    const base = s.tetIndex[node] * 4;
    if (tets[base] === third || tets[base + 1] === third ||
        tets[base + 2] === third || tets[base + 3] === third) count++;
  }
  return count;
}

// A boundary triangle containing this edge means its dual face extends to infinity.
function isClosed(s: Scratch, face: number, tets: ArrayLike<number>): boolean {
  const a = s.edgeA[face], b = s.edgeB[face];
  for (let node = s.head[face]; node >= 0; node = s.next[node]) {
    const base = s.tetIndex[node] * 4;
    let c = -1, d = -1, count = 0;
    for (let i = 0; i < 4; i++) {
      const vertex = tets[base + i];
      if (vertex !== a && vertex !== b) {
        if (count++ === 0) c = vertex;
        else d = vertex;
      }
    }
    if (count !== 2 || c === d || countThird(s, face, tets, c) !== 2 ||
        countThird(s, face, tets, d) !== 2) return false;
  }
  return true;
}

/**
 * Triangulate bounded Voronoi faces between filled and empty sites. The input
 * tetrahedralization and circumcenters are supplied by the caller; no Delaunay
 * solver is included. Reuse scratch and mesh for allocation-free steady-state builds.
 */
export function buildVoronoiHull(input: VoronoiInput, scratch: Scratch, mesh: MeshBuffer): MeshBuffer {
  const { sites, tetrahedra: tets, circumcenters: centers, filled, emitCell } = input;
  if (sites.length % 3 || tets.length % 4 || centers.length !== tets.length / 4 * 3 ||
      filled.length !== sites.length / 3 || (emitCell && emitCell.length !== filled.length))
    throw new RangeError("Input array lengths do not match");
  const epsilon = input.duplicateEpsilon ?? 1e-5;
  if (epsilon < 0) throw new RangeError("duplicateEpsilon must be nonnegative");
  const epsilonSq = epsilon * epsilon;
  const siteCount = sites.length / 3;
  const tetCount = tets.length / 4;
  mesh.vertexCount = mesh.indexCount = 0;
  scratch.reset(tetCount * 6);

  for (let t = 0; t < tetCount; t++) {
    const i = t * 4;
    const a = tets[i], b = tets[i + 1], c = tets[i + 2], d = tets[i + 3];
    if (!Number.isInteger(a) || !Number.isInteger(b) || !Number.isInteger(c) || !Number.isInteger(d) ||
        Math.min(a, b, c, d) < 0 || Math.max(a, b, c, d) >= siteCount)
      throw new RangeError("Tetrahedron site index is out of range");
    addEdge(scratch, a, b, t); addEdge(scratch, a, c, t); addEdge(scratch, a, d, t);
    addEdge(scratch, b, c, t); addEdge(scratch, b, d, t); addEdge(scratch, c, d, t);
  }

  for (let face = 0; face < scratch.faceCount; face++) {
    const a = scratch.edgeA[face], b = scratch.edgeB[face];
    if (!!filled[a] === !!filled[b]) continue;
    const source = filled[a] ? a : b;
    if (emitCell && !emitCell[source]) continue;
    if (!isClosed(scratch, face, tets)) continue;

    const src = source * 3, dst = (source === a ? b : a) * 3;
    let nx = sites[dst] - sites[src], ny = sites[dst + 1] - sites[src + 1],
        nz = sites[dst + 2] - sites[src + 2];
    const normalLength = Math.hypot(nx, ny, nz);
    if (!(normalLength > 1e-10)) continue;
    nx /= normalLength; ny /= normalLength; nz /= normalLength;

    let count = 0, cx = 0, cy = 0, cz = 0;
    for (let node = scratch.head[face]; node >= 0; node = scratch.next[node]) {
      const ci = scratch.tetIndex[node] * 3;
      const x = centers[ci], y = centers[ci + 1], z = centers[ci + 2];
      if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) continue;
      let duplicate = false;
      for (let i = 0; i < count; i++) {
        const p = i * 3, dx = scratch.points[p] - x, dy = scratch.points[p + 1] - y,
              dz = scratch.points[p + 2] - z;
        if (dx * dx + dy * dy + dz * dz <= epsilonSq) { duplicate = true; break; }
      }
      if (duplicate) continue;
      scratch.reservePoints(count + 1);
      const p = count++ * 3;
      scratch.points[p] = x; scratch.points[p + 1] = y; scratch.points[p + 2] = z;
      cx += x; cy += y; cz += z;
    }
    if (count < 3) continue;
    cx /= count; cy /= count; cz /= count;

    // tangent = normalize(cross(normal, X or Y)); bitangent = cross(normal, tangent)
    const ax = Math.abs(nx) < 0.9 ? 1 : 0, ay = ax ? 0 : 1;
    let tx = -nz * ay, ty = nz * ax, tz = nx * ay - ny * ax;
    const tangentLength = Math.hypot(tx, ty, tz);
    tx /= tangentLength; ty /= tangentLength; tz /= tangentLength;
    const bx = ny * tz - nz * ty, by = nz * tx - nx * tz, bz = nx * ty - ny * tx;
    for (let i = 0; i < count; i++) {
      const p = i * 3, x = scratch.points[p] - cx, y = scratch.points[p + 1] - cy,
            z = scratch.points[p + 2] - cz;
      scratch.angles[i] = Math.atan2(x * bx + y * by + z * bz, x * tx + y * ty + z * tz);
    }
    for (let i = 1; i < count; i++) {
      const p = i * 3, x = scratch.points[p], y = scratch.points[p + 1],
            z = scratch.points[p + 2], angle = scratch.angles[i];
      let j = i - 1;
      while (j >= 0 && scratch.angles[j] > angle) {
        const dstp = (j + 1) * 3, srcp = j * 3;
        scratch.points[dstp] = scratch.points[srcp];
        scratch.points[dstp + 1] = scratch.points[srcp + 1];
        scratch.points[dstp + 2] = scratch.points[srcp + 2];
        scratch.angles[j + 1] = scratch.angles[j];
        j--;
      }
      const dstp = (j + 1) * 3;
      scratch.points[dstp] = x; scratch.points[dstp + 1] = y; scratch.points[dstp + 2] = z;
      scratch.angles[j + 1] = angle;
    }

    let signedArea = 0;
    const x0 = scratch.points[0], y0 = scratch.points[1], z0 = scratch.points[2];
    for (let i = 1; i < count - 1; i++) {
      const p = i * 3, q = (i + 1) * 3;
      const ux = scratch.points[p] - x0, uy = scratch.points[p + 1] - y0, uz = scratch.points[p + 2] - z0;
      const vx = scratch.points[q] - x0, vy = scratch.points[q + 1] - y0, vz = scratch.points[q + 2] - z0;
      signedArea += (uy * vz - uz * vy) * nx + (uz * vx - ux * vz) * ny + (ux * vy - uy * vx) * nz;
    }
    if (Math.abs(signedArea) <= epsilonSq) continue;
    const first = mesh.vertexCount;
    mesh.reserve(first + count, mesh.indexCount + (count - 2) * 3);
    for (let i = 0; i < count; i++) {
      const from = (signedArea < 0 ? count - 1 - i : i) * 3;
      const to = mesh.vertexCount++ * 3;
      mesh.positions[to] = scratch.points[from];
      mesh.positions[to + 1] = scratch.points[from + 1];
      mesh.positions[to + 2] = scratch.points[from + 2];
      mesh.normals[to] = nx; mesh.normals[to + 1] = ny; mesh.normals[to + 2] = nz;
    }
    for (let i = 1; i < count - 1; i++) {
      mesh.indices[mesh.indexCount++] = first;
      mesh.indices[mesh.indexCount++] = first + i;
      mesh.indices[mesh.indexCount++] = first + i + 1;
    }
  }
  return mesh;
}
