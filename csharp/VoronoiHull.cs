using System;
using System.Collections.Generic;
using System.Numerics;

namespace VoronoiHullStandalone
{
    public readonly struct Tetrahedron
    {
        public readonly int A, B, C, D;
        public Tetrahedron(int a, int b, int c, int d) { A = a; B = b; C = c; D = d; }
    }

    /// <summary>Reusable, grow-only mesh storage. Only the first VertexCount and IndexCount entries are valid.</summary>
    public sealed class MeshBuffer
    {
        public Vector3[] Positions = Array.Empty<Vector3>();
        public Vector3[] Normals = Array.Empty<Vector3>();
        public int[] Indices = Array.Empty<int>();
        public int VertexCount { get; internal set; }
        public int IndexCount { get; internal set; }

        internal void Reserve(int vertices, int indices)
        {
            if (Positions.Length < vertices)
            {
                int capacity = Grow(Positions.Length, vertices);
                Array.Resize(ref Positions, capacity);
                Array.Resize(ref Normals, capacity);
            }
            if (Indices.Length < indices) Array.Resize(ref Indices, Grow(Indices.Length, indices));
        }

        private static int Grow(int oldCapacity, int needed)
        {
            int capacity = Math.Max(16, oldCapacity);
            while (capacity < needed) capacity = checked(capacity * 2);
            return capacity;
        }
    }

    /// <summary>Reusable topology and polygon workspace. Use one instance per concurrent build.</summary>
    public sealed class Scratch
    {
        internal readonly Dictionary<ulong, int> EdgeToFace = new Dictionary<ulong, int>();
        internal int[] EdgeA = Array.Empty<int>(), EdgeB = Array.Empty<int>();
        internal int[] Head = Array.Empty<int>(), Next = Array.Empty<int>(), TetIndex = Array.Empty<int>();
        internal Vector3[] Points = Array.Empty<Vector3>();
        internal float[] Angles = Array.Empty<float>();
        internal int FaceCount, NodeCount;

        internal void Reset(int maxEdges)
        {
            EdgeToFace.Clear();
            FaceCount = NodeCount = 0;
            if (EdgeA.Length < maxEdges)
            {
                Array.Resize(ref EdgeA, maxEdges);
                Array.Resize(ref EdgeB, maxEdges);
                Array.Resize(ref Head, maxEdges);
                Array.Resize(ref Next, maxEdges);
                Array.Resize(ref TetIndex, maxEdges);
            }
        }

        internal void ReservePoints(int count)
        {
            if (Points.Length >= count) return;
            int capacity = Math.Max(8, Points.Length);
            while (capacity < count) capacity = checked(capacity * 2);
            Array.Resize(ref Points, capacity);
            Array.Resize(ref Angles, capacity);
        }
    }

    public static class VoronoiHull
    {
        /// <summary>
        /// Builds the boundary between filled and empty sites from an existing 3D Delaunay
        /// tetrahedralization and its circumcenters. A face is emitted only when its filled
        /// endpoint is selected by emitCell (or for every filled endpoint when emitCell is null).
        /// Reuse scratch and mesh to avoid allocations after capacities have warmed up.
        /// </summary>
        public static MeshBuffer Build(
            Vector3[] sites, Tetrahedron[] tetrahedra, Vector3[] circumcenters,
            bool[] filled, Scratch scratch, MeshBuffer mesh, bool[] emitCell = null,
            float duplicateEpsilon = 1e-5f)
        {
            if (sites == null || tetrahedra == null || circumcenters == null || filled == null || scratch == null || mesh == null)
                throw new ArgumentNullException("Build inputs cannot be null");
            if (sites.Length != filled.Length || (emitCell != null && emitCell.Length != sites.Length) || tetrahedra.Length != circumcenters.Length)
                throw new ArgumentException("Input array lengths do not match");
            if (duplicateEpsilon < 0) throw new ArgumentOutOfRangeException(nameof(duplicateEpsilon));

            mesh.VertexCount = mesh.IndexCount = 0;
            scratch.Reset(checked(tetrahedra.Length * 6));
            for (int t = 0; t < tetrahedra.Length; t++)
            {
                Tetrahedron tet = tetrahedra[t];
                CheckSite(tet.A, sites.Length); CheckSite(tet.B, sites.Length);
                CheckSite(tet.C, sites.Length); CheckSite(tet.D, sites.Length);
                AddEdge(scratch, tet.A, tet.B, t); AddEdge(scratch, tet.A, tet.C, t);
                AddEdge(scratch, tet.A, tet.D, t); AddEdge(scratch, tet.B, tet.C, t);
                AddEdge(scratch, tet.B, tet.D, t); AddEdge(scratch, tet.C, tet.D, t);
            }

            float epsilonSq = duplicateEpsilon * duplicateEpsilon;
            for (int face = 0; face < scratch.FaceCount; face++)
            {
                int a = scratch.EdgeA[face], b = scratch.EdgeB[face];
                if (filled[a] == filled[b]) continue;
                int source = filled[a] ? a : b;
                if (emitCell != null && !emitCell[source]) continue;
                if (!IsClosed(scratch, face, tetrahedra)) continue; // unbounded Voronoi face

                Vector3 direction = sites[filled[a] ? b : a] - sites[source];
                float lengthSq = direction.LengthSquared();
                if (lengthSq <= 1e-20f) continue;
                Vector3 normal = direction / MathF.Sqrt(lengthSq);

                int count = 0;
                Vector3 centroid = Vector3.Zero;
                for (int node = scratch.Head[face]; node >= 0; node = scratch.Next[node])
                {
                    Vector3 point = circumcenters[scratch.TetIndex[node]];
                    if (!Finite(point)) continue;
                    bool duplicate = false;
                    for (int i = 0; i < count; i++)
                        if (Vector3.DistanceSquared(scratch.Points[i], point) <= epsilonSq) { duplicate = true; break; }
                    if (duplicate) continue;
                    scratch.ReservePoints(count + 1);
                    scratch.Points[count++] = point;
                    centroid += point;
                }
                if (count < 3) continue;
                centroid /= count;

                Vector3 tangent = Vector3.Normalize(Vector3.Cross(normal,
                    MathF.Abs(normal.X) < 0.9f ? Vector3.UnitX : Vector3.UnitY));
                Vector3 bitangent = Vector3.Cross(normal, tangent);
                for (int i = 0; i < count; i++)
                {
                    Vector3 delta = scratch.Points[i] - centroid;
                    scratch.Angles[i] = MathF.Atan2(Vector3.Dot(delta, bitangent), Vector3.Dot(delta, tangent));
                }
                // Insertion sort is efficient for the small polygons around Delaunay edges.
                for (int i = 1; i < count; i++)
                {
                    Vector3 point = scratch.Points[i]; float angle = scratch.Angles[i];
                    int j = i - 1;
                    while (j >= 0 && scratch.Angles[j] > angle)
                    {
                        scratch.Points[j + 1] = scratch.Points[j];
                        scratch.Angles[j + 1] = scratch.Angles[j];
                        j--;
                    }
                    scratch.Points[j + 1] = point; scratch.Angles[j + 1] = angle;
                }

                // A valid polygon can still contain nearly collinear circumcenters.
                float signedArea = 0;
                for (int i = 1; i < count - 1; i++)
                    signedArea += Vector3.Dot(Vector3.Cross(scratch.Points[i] - scratch.Points[0],
                        scratch.Points[i + 1] - scratch.Points[0]), normal);
                if (MathF.Abs(signedArea) <= epsilonSq) continue;
                if (signedArea < 0) Array.Reverse(scratch.Points, 0, count);

                int first = mesh.VertexCount;
                mesh.Reserve(checked(first + count), checked(mesh.IndexCount + (count - 2) * 3));
                for (int i = 0; i < count; i++)
                {
                    mesh.Positions[mesh.VertexCount] = scratch.Points[i];
                    mesh.Normals[mesh.VertexCount++] = normal;
                }
                for (int i = 1; i < count - 1; i++)
                {
                    mesh.Indices[mesh.IndexCount++] = first;
                    mesh.Indices[mesh.IndexCount++] = first + i;
                    mesh.Indices[mesh.IndexCount++] = first + i + 1;
                }
            }
            return mesh;
        }

        private static void CheckSite(int index, int count)
        {
            if ((uint)index >= (uint)count) throw new ArgumentOutOfRangeException(nameof(index), "Tetrahedron site index is out of range");
        }

        private static bool Finite(Vector3 p) =>
            !float.IsNaN(p.X) && !float.IsNaN(p.Y) && !float.IsNaN(p.Z) &&
            !float.IsInfinity(p.X) && !float.IsInfinity(p.Y) && !float.IsInfinity(p.Z);

        private static void AddEdge(Scratch s, int a, int b, int tet)
        {
            if (a == b) return;
            if (a > b) { int swap = a; a = b; b = swap; }
            ulong key = ((ulong)(uint)a << 32) | (uint)b;
            if (!s.EdgeToFace.TryGetValue(key, out int face))
            {
                face = s.FaceCount++;
                s.EdgeToFace.Add(key, face);
                s.EdgeA[face] = a; s.EdgeB[face] = b; s.Head[face] = -1;
            }
            int node = s.NodeCount++;
            s.TetIndex[node] = tet;
            s.Next[node] = s.Head[face];
            s.Head[face] = node;
        }

        // Every triangle containing this edge must belong to two tetrahedra.
        // Otherwise the dual Voronoi face extends to infinity and needs clipping.
        private static bool IsClosed(Scratch s, int face, Tetrahedron[] tetrahedra)
        {
            int a = s.EdgeA[face], b = s.EdgeB[face];
            for (int node = s.Head[face]; node >= 0; node = s.Next[node])
            {
                Tetrahedron t = tetrahedra[s.TetIndex[node]];
                if (!OtherVertices(t, a, b, out int c, out int d)) return false;
                if (CountThird(s, face, tetrahedra, c) != 2 || CountThird(s, face, tetrahedra, d) != 2)
                    return false;
            }
            return true;
        }

        private static int CountThird(Scratch s, int face, Tetrahedron[] tetrahedra, int third)
        {
            int count = 0;
            for (int node = s.Head[face]; node >= 0; node = s.Next[node])
            {
                Tetrahedron t = tetrahedra[s.TetIndex[node]];
                if (t.A == third || t.B == third || t.C == third || t.D == third) count++;
            }
            return count;
        }

        private static bool OtherVertices(Tetrahedron t, int a, int b, out int c, out int d)
        {
            c = d = -1; int count = 0;
            if (t.A != a && t.A != b) { c = t.A; count++; }
            if (t.B != a && t.B != b) { if (count++ == 0) c = t.B; else d = t.B; }
            if (t.C != a && t.C != b) { if (count++ == 0) c = t.C; else d = t.C; }
            if (t.D != a && t.D != b) { if (count++ == 0) c = t.D; else d = t.D; }
            return count == 2 && c != d;
        }
    }
}
