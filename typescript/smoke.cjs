const assert = require("node:assert/strict");
const { buildVoronoiHull, buildVoronoiHullFromSites, VoronoiWorkspace, Scratch, MeshBuffer } = require("./dist/voronoiHull.js");

// Eight tetrahedra fill an octahedron around the center site. Its Voronoi cell is a cube.
const sites = new Float32Array([0,0,0, 1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1]);
const tets = [];
const centers = [];
for (const x of [1, 2]) for (const y of [3, 4]) for (const z of [5, 6]) {
  tets.push(0, x, y, z);
  centers.push(sites[x*3] / 2, sites[y*3+1] / 2, sites[z*3+2] / 2);
}
const input = { sites, tetrahedra: new Uint32Array(tets), circumcenters: new Float32Array(centers),
  filled: new Uint8Array([1,0,0,0,0,0,0]) };
const scratch = new Scratch(), mesh = new MeshBuffer();
buildVoronoiHull(input, scratch, mesh);
assert.equal(mesh.vertexCount, 24);
assert.equal(mesh.indexCount, 36);
const firstPositions = mesh.positions, firstIndices = mesh.indices;
for (let i = 0; i < mesh.indexCount; i += 3) {
  const a = mesh.indices[i] * 3, b = mesh.indices[i+1] * 3, c = mesh.indices[i+2] * 3;
  const ux = mesh.positions[b]-mesh.positions[a], uy = mesh.positions[b+1]-mesh.positions[a+1], uz = mesh.positions[b+2]-mesh.positions[a+2];
  const vx = mesh.positions[c]-mesh.positions[a], vy = mesh.positions[c+1]-mesh.positions[a+1], vz = mesh.positions[c+2]-mesh.positions[a+2];
  const dot = (uy*vz-uz*vy)*mesh.normals[a] + (uz*vx-ux*vz)*mesh.normals[a+1] + (ux*vy-uy*vx)*mesh.normals[a+2];
  assert.ok(dot > 0, "triangle winding must match outward normal");
}
buildVoronoiHull(input, scratch, mesh);
assert.equal(mesh.positions, firstPositions);
assert.equal(mesh.indices, firstIndices);
input.filled[1] = 1;
buildVoronoiHull(input, scratch, mesh);
assert.equal(mesh.vertexCount, 20);
assert.equal(mesh.indexCount, 30);
input.emitCell = new Uint8Array(7);
buildVoronoiHull(input, scratch, mesh);
assert.equal(mesh.indexCount, 0);

// The direct site-cloud path includes Delaunay tetrahedralization.
const genericSites = new Float64Array([0,0,0, 1,0.03,0.02, -1.17,0.01,-0.03,
  0.02,1.19,0.04, -0.02,-0.94,0.01, 0.01,-0.02,1.08, -0.02,0.03,-1.26]);
const workspace = new VoronoiWorkspace();
const genericFilled = new Uint8Array([1,0,0,0,0,0,0]);
buildVoronoiHullFromSites(genericSites, genericFilled, workspace);
assert.equal(workspace.delaunay.tetrahedronCount, 8);
assert.equal(workspace.mesh.vertexCount, 24);
assert.equal(workspace.mesh.indexCount, 36);
const tetsBuffer = workspace.delaunay.tetrahedra;
const positionsBuffer = workspace.mesh.positions;
buildVoronoiHullFromSites(genericSites, genericFilled, workspace);
assert.equal(workspace.delaunay.tetrahedra, tetsBuffer);
assert.equal(workspace.mesh.positions, positionsBuffer);
genericFilled[1] = 1;
buildVoronoiHullFromSites(genericSites, genericFilled, workspace);
assert.equal(workspace.mesh.indexCount, 30);
console.log("TypeScript smoke checks passed");
