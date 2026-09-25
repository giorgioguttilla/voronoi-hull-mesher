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

// The direct site-cloud path includes Delaunay tetrahedralization.
var genericSites = new[] { Vector3.Zero, new Vector3(1, .03f, .02f),
    new Vector3(-1.17f, .01f, -.03f), new Vector3(.02f, 1.19f, .04f),
    new Vector3(-.02f, -.94f, .01f), new Vector3(.01f, -.02f, 1.08f),
    new Vector3(-.02f, .03f, -1.26f) };
var genericFilled = new bool[7]; genericFilled[0] = true;
var workspace = new VoronoiWorkspace();
VoronoiHull.BuildFromSites(genericSites, genericFilled, workspace);
Check(workspace.Delaunay.Count == 8 && workspace.Mesh.VertexCount == 24 && workspace.Mesh.IndexCount == 36,
    "point cloud to hull");
var originalTets = workspace.Delaunay.Tetrahedra;
var genericPositions = workspace.Mesh.Positions;
VoronoiHull.BuildFromSites(genericSites, genericFilled, workspace);
Check(ReferenceEquals(workspace.Delaunay.Tetrahedra, originalTets) &&
    ReferenceEquals(workspace.Mesh.Positions, genericPositions), "pipeline buffer reuse");
genericFilled[1] = true;
VoronoiHull.BuildFromSites(genericSites, genericFilled, workspace);
Check(workspace.Mesh.IndexCount == 30, "pipeline occluded face culling");
Console.WriteLine("C# smoke checks passed");

static void Check(bool condition, string description)
{
    if (!condition) throw new Exception(description);
}
