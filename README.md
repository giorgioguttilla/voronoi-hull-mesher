# Voronoi hull mesher

Small, dependency free C# and TypeScript implementations of a 3D Voronoi **boundary mesh** generator. The main API takes a point cloud and a `filled` mask, computes a Delaunay tetrahedralization, and emits only the Voronoi faces between filled and empty cells. Faces between two filled cells are omitted. This is useful for a chunked world where adjacent cells should not generate hidden triangles.

The Delaunay step uses incremental Bowyer–Watson tetrahedralization. The hull step gathers the circumcenters of tetrahedra around each Delaunay edge, sorts them into a polygon, and triangulates exposed faces. Both implementations are ordinary C# and TypeScript, with no runtime dependencies beyond their standard libraries. For the circumcenter mathematics, see [Rodolphe Vaillant's derivation](https://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter).

## Inputs and output

| Input | Meaning |
| --- | --- |
| `sites` | One XYZ position per Voronoi site. |
| `filled` | Whether each site's cell belongs to the solid region. |
| `emitCell` | Optional ownership mask. Emit a boundary face only if its filled site belongs to this batch or chunk. |

The main API computes tetrahedra and circumcenters internally. A lower level API accepts precomputed tetrahedra and circumcenters if you want to reuse a triangulation. Both implementations return positions, flat normals, triangle indices, and valid element counts. Positions are in the same coordinate space as the input sites. The output is reset and reused on each call; copy it before retaining it across builds. Capacity growth allocates, but repeated calls at or below established capacity do not create per face arrays or objects. Use one workspace per concurrent build.

Only **bounded** Voronoi faces are emitted. A Delaunay edge on the convex hull has an unbounded dual face and needs a separate clipping step. For a closed mesh, provide extra sites outside the region of interest so all relevant faces are bounded. The incremental Delaunay builder is intended for finite, distinct sites in general 3D position and modest batches; nearly coplanar or cospherical point sets may need a more robust exact-predicate solver. Nearly coincident circumcenters are merged using `duplicateEpsilon` (default `1e-5` in input units).

## C#

The C# project targets .NET Standard 2.1 and uses `System.Numerics.Vector3`. Add `csharp/Delaunay3D.cs` and `csharp/VoronoiHull.cs` to a project or compile `csharp/VoronoiHullStandalone.csproj`.

```csharp
using System.Numerics;
using VoronoiHullStandalone;

Vector3[] sites = {
    Vector3.Zero, new Vector3(1, .03f, .02f),
    new Vector3(-1.17f, .01f, -.03f), new Vector3(.02f, 1.19f, .04f),
    new Vector3(-.02f, -.94f, .01f), new Vector3(.01f, -.02f, 1.08f),
    new Vector3(-.02f, .03f, -1.26f)
};
bool[] filled = { true, false, false, false, false, false, false };

var workspace = new VoronoiWorkspace(); // keep this between builds
MeshBuffer mesh = VoronoiHull.BuildFromSites(sites, filled, workspace);

// Read mesh.Positions[0..mesh.VertexCount], mesh.Normals[0..mesh.VertexCount],
// and mesh.Indices[0..mesh.IndexCount].
```

## TypeScript / Three.js

The TypeScript functions use flat typed arrays and have no runtime dependencies. The mesh output is ready for `BufferGeometry`:

```ts
import * as THREE from "three";
import { buildVoronoiHullFromSites, VoronoiWorkspace } from "./voronoiHull";

const sites = new Float32Array([
  0,0,0, 1,0.03,0.02, -1.17,0.01,-0.03,
  0.02,1.19,0.04, -0.02,-0.94,0.01,
  0.01,-0.02,1.08, -0.02,0.03,-1.26
]);
const filled = new Uint8Array([1,0,0,0,0,0,0]);
const workspace = new VoronoiWorkspace(); // keep this between builds
const mesh = buildVoronoiHullFromSites(sites, filled, workspace);

const geometry = new THREE.BufferGeometry();
// slice() gives Three.js independent storage that survives the next build.
geometry.setAttribute("position", new THREE.BufferAttribute(mesh.positions.slice(0, mesh.vertexCount * 3), 3));
geometry.setAttribute("normal", new THREE.BufferAttribute(mesh.normals.slice(0, mesh.vertexCount * 3), 3));
geometry.setIndex(new THREE.BufferAttribute(mesh.indices.slice(0, mesh.indexCount), 1));
```

If you update a single Three.js geometry in place, you can instead keep the reusable arrays and mark attributes `needsUpdate`, taking care to refresh their draw ranges when counts change.

For precomputed Delaunay data, call `Delaunay3D.Build` and `VoronoiHull.Build` in C#, or `buildDelaunay3D` and `buildVoronoiHull` in TypeScript. Pass the valid tetrahedron count because the result arrays retain extra capacity.

## Algorithm

1. Incrementally tetrahedralize the sites using Bowyer–Watson and compute each tetrahedron's circumcenter.
2. Index each unique Delaunay edge and gather its incident tetrahedra.
3. Skip edges whose two cells have the same `filled` state, or whose filled cell is outside `emitCell`.
4. Skip unbounded faces by checking whether every triangle around the edge has two incident tetrahedra.
5. Gather unique circumcenters, sort them around the edge axis, and triangulate the polygon with outward winding.

The result has one copy of each exposed polygon. It does not weld vertices across faces, which preserves flat normals and keeps writes simple. The simple Delaunay builder scans all current tetrahedra for each inserted site, so it is best for modest batches. The hull pass is linear in tetrahedron count plus sorting of each face's small incident ring; its conservative bounded face check scans each ring quadratically.

## Verify

```sh
dotnet run --project csharp/Smoke/Smoke.csproj
tsc -p typescript/tsconfig.json && node typescript/smoke.cjs
```
