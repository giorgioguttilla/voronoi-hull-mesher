# Voronoi hull mesher

Small, dependency free C# and TypeScript implementations of a 3D Voronoi **boundary mesh** generator. They take sites, a Delaunay tetrahedralization, and one circumcenter per tetrahedron. They build a polygon around each Delaunay edge whose sites have different `filled` values, then triangulate that polygon. Faces between two filled sites are omitted. This is useful for a voxel or chunk style world where adjacent Voronoi cells should not generate hidden triangles.

The C# code is adapted from the `CellGeneration.VoronoiHull` Burst job in Starlog. The project specific material, liquid, grass, planet, and half cell logic has been removed. These versions use ordinary C# and TypeScript and do **not** include a Delaunay solver.

## Inputs and output

| Input | Meaning |
| --- | --- |
| `sites` | One XYZ position per Voronoi site. |
| `tetrahedra` | Four site indices per Delaunay tetrahedron. |
| `circumcenters` | One XYZ circumcenter per tetrahedron, in the same order. |
| `filled` | Whether each site's cell belongs to the solid region. |
| `emitCell` | Optional ownership mask. Emit a boundary face only if its filled site belongs to this batch or chunk. |

Both implementations return positions, flat normals, triangle indices, and valid element counts. Positions are in the same coordinate space as input sites and circumcenters. The output is reset and reused on each call; copy it before retaining it across builds. Scratch storage is also reused. Capacity growth allocates, but repeated calls at or below established capacity do not create per face arrays or objects. Use one scratch and output pair per concurrent build.

Only **bounded** Voronoi faces are emitted. A Delaunay edge on the convex hull has an unbounded dual face and needs a separate clipping step. For a closed mesh, provide extra sites outside the region of interest so all relevant faces are bounded. Degenerate tetrahedralizations and duplicate tetrahedra are outside the supported input contract. Nearly coincident circumcenters are merged using `duplicateEpsilon` (default `1e-5` in input units).

## C#

`csharp/VoronoiHull.cs` targets .NET Standard 2.1 and uses `System.Numerics.Vector3`. Add it to a project or compile `csharp/VoronoiHullStandalone.csproj`.

```csharp
using System.Numerics;
using VoronoiHullStandalone;

Vector3[] sites = /* your site coordinates */;
Tetrahedron[] tets = /* your Delaunay tetrahedra */;
Vector3[] circumcenters = /* one per tet */;
bool[] filled = /* one per site */;
bool[] emitCell = /* optional: sites owned by this chunk */;

var scratch = new Scratch();
var mesh = new MeshBuffer();
VoronoiHull.Build(sites, tets, circumcenters, filled, scratch, mesh, emitCell);

// Read mesh.Positions[0..mesh.VertexCount], mesh.Normals[0..mesh.VertexCount],
// and mesh.Indices[0..mesh.IndexCount].
```

## TypeScript / Three.js

`typescript/voronoiHull.ts` has no runtime dependencies. Arrays are flat and typed, ready for `BufferGeometry`:

```ts
import * as THREE from "three";
import { buildVoronoiHull, Scratch, MeshBuffer } from "./voronoiHull";

const scratch = new Scratch();
const mesh = new MeshBuffer();
buildVoronoiHull({
  sites: new Float32Array([/* xyz per site */]),
  tetrahedra: new Uint32Array([/* abcd per tet */]),
  circumcenters: new Float32Array([/* xyz per tet */]),
  filled: new Uint8Array([/* 0 or 1 per site */]),
  // emitCell: new Uint8Array([/* ownership mask */]),
}, scratch, mesh);

const geometry = new THREE.BufferGeometry();
// slice() gives Three.js independent storage that survives the next build.
geometry.setAttribute("position", new THREE.BufferAttribute(mesh.positions.slice(0, mesh.vertexCount * 3), 3));
geometry.setAttribute("normal", new THREE.BufferAttribute(mesh.normals.slice(0, mesh.vertexCount * 3), 3));
geometry.setIndex(new THREE.BufferAttribute(mesh.indices.slice(0, mesh.indexCount), 1));
```

If you update a single Three.js geometry in place, you can instead keep the reusable arrays and mark attributes `needsUpdate`, taking care to refresh their draw ranges when counts change.

## Algorithm

1. Index each unique Delaunay edge and gather its incident tetrahedra.
2. Skip edges whose two cells have the same `filled` state, or whose filled cell is outside `emitCell`.
3. Skip unbounded faces by checking whether every triangle around the edge has two incident tetrahedra.
4. Gather unique circumcenters, sort them around the edge axis, and triangulate the polygon with outward winding.

The result has one copy of each exposed polygon. It does not weld vertices across faces, which preserves flat normals and keeps writes simple. Expected work is linear in tetrahedron count plus sorting of each face's small incident ring; the conservative bounded face check scans each ring quadratically.

## Verify

```sh
dotnet run --project csharp/Smoke/Smoke.csproj
tsc -p typescript/tsconfig.json && node typescript/smoke.cjs
```
