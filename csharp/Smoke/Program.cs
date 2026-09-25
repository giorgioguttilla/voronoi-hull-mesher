using System;
using System.Numerics;
using VoronoiHullStandalone;

// Eight tetrahedra fill an octahedron around the center site. Its Voronoi cell is a cube.
var sites = new[] { Vector3.Zero, Vector3.UnitX, -Vector3.UnitX,
    Vector3.UnitY, -Vector3.UnitY, Vector3.UnitZ, -Vector3.UnitZ };
var tets = new Tetrahedron[8];
var centers = new Vector3[8];
int t = 0;
foreach (int x in new[] { 1, 2 }) foreach (int y in new[] { 3, 4 }) foreach (int z in new[] { 5, 6 })
{
    tets[t] = new Tetrahedron(0, x, y, z);
    centers[t++] = (sites[x] + sites[y] + sites[z]) / 2;
}
var filled = new bool[7]; filled[0] = true;
var scratch = new Scratch(); var mesh = new MeshBuffer();
VoronoiHull.Build(sites, tets, centers, filled, scratch, mesh);
Check(mesh.VertexCount == 24 && mesh.IndexCount == 36, "cube counts");
var originalPositions = mesh.Positions; var originalIndices = mesh.Indices;
for (int i = 0; i < mesh.IndexCount; i += 3)
{
    var a = mesh.Positions[mesh.Indices[i]];
    var b = mesh.Positions[mesh.Indices[i + 1]];
    var c = mesh.Positions[mesh.Indices[i + 2]];
    var normal = mesh.Normals[mesh.Indices[i]];
    Check(Vector3.Dot(Vector3.Cross(b - a, c - a), normal) > 0, "outward winding");
}
VoronoiHull.Build(sites, tets, centers, filled, scratch, mesh);
Check(ReferenceEquals(mesh.Positions, originalPositions) && ReferenceEquals(mesh.Indices, originalIndices), "buffer reuse");
filled[1] = true;
VoronoiHull.Build(sites, tets, centers, filled, scratch, mesh);
Check(mesh.VertexCount == 20 && mesh.IndexCount == 30, "occluded face culling");
VoronoiHull.Build(sites, tets, centers, filled, scratch, mesh, new bool[7]);
Check(mesh.IndexCount == 0, "ownership mask");
Console.WriteLine("C# smoke checks passed");

static void Check(bool condition, string description)
{
    if (!condition) throw new Exception(description);
}
