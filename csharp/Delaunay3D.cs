using System;
using System.Collections.Generic;
using System.Numerics;

namespace VoronoiHullStandalone
{
    /// <summary>Reusable Delaunay output. Only the first Count entries are valid.</summary>
    public sealed class DelaunayResult
    {
        public Tetrahedron[] Tetrahedra = Array.Empty<Tetrahedron>();
        public Vector3[] Circumcenters = Array.Empty<Vector3>();
        public int Count { get; internal set; }

        internal void Reserve(int count)
        {
            if (Tetrahedra.Length >= count) return;
            int capacity = Math.Max(16, Tetrahedra.Length);
            while (capacity < count) capacity = checked(capacity * 2);
            Array.Resize(ref Tetrahedra, capacity);
            Array.Resize(ref Circumcenters, capacity);
        }

        internal void RemoveSwapBack(int index)
        {
            int last = --Count;
            Tetrahedra[index] = Tetrahedra[last];
            Circumcenters[index] = Circumcenters[last];
        }
    }

    /// <summary>Reusable workspace for the incremental tetrahedralizer.</summary>
    public sealed class DelaunayScratch
    {
        internal Vector3[] Points = Array.Empty<Vector3>();
        internal readonly Dictionary<FaceKey, int> FaceToIndex = new Dictionary<FaceKey, int>();
        internal FaceKey[] Faces = Array.Empty<FaceKey>();
        internal bool[] FaceActive = Array.Empty<bool>();
        internal int FaceCount;

        internal void ReservePoints(int count)
        {
            if (Points.Length < count) Array.Resize(ref Points, count);
        }

        internal void ResetFaces(int maxFaces)
        {
            FaceToIndex.Clear();
            FaceCount = 0;
            if (Faces.Length < maxFaces)
            {
                Array.Resize(ref Faces, maxFaces);
                Array.Resize(ref FaceActive, maxFaces);
            }
        }

        internal void ToggleFace(int a, int b, int c)
        {
            FaceKey key = new FaceKey(a, b, c);
            if (FaceToIndex.TryGetValue(key, out int index))
                FaceActive[index] = false;
            else
            {
                index = FaceCount++;
                FaceToIndex.Add(key, index);
                Faces[index] = key;
                FaceActive[index] = true;
            }
        }
    }

    internal readonly struct FaceKey : IEquatable<FaceKey>
    {
        internal readonly int A, B, C;

        internal FaceKey(int a, int b, int c)
        {
            if (a > b) { int swap = a; a = b; b = swap; }
            if (b > c) { int swap = b; b = c; c = swap; }
            if (a > b) { int swap = a; a = b; b = swap; }
            A = a; B = b; C = c;
        }

        public bool Equals(FaceKey other) => A == other.A && B == other.B && C == other.C;
        public override bool Equals(object other) => other is FaceKey key && Equals(key);
        public override int GetHashCode()
        {
            unchecked { return ((A * 397) ^ B) * 397 ^ C; }
        }
    }

    /// <summary>
    /// Incremental Bowyer-Watson Delaunay tetrahedralization. Intended for finite,
    /// distinct 3D sites in general position and modest batch sizes. Reuse scratch
    /// and output to avoid steady-state heap allocations; capacity growth may allocate.
    /// </summary>
    public static class Delaunay3D
    {
        public static DelaunayResult Build(Vector3[] sites, DelaunayScratch scratch, DelaunayResult output)
        {
            if (sites == null || scratch == null || output == null)
                throw new ArgumentNullException("Build inputs cannot be null");
            if (sites.Length < 4) throw new ArgumentException("At least four noncoplanar sites are required", nameof(sites));
            output.Count = 0;
            scratch.ReservePoints(checked(sites.Length + 4));
            Array.Copy(sites, scratch.Points, sites.Length);

            Vector3 min = sites[0], max = sites[0];
            for (int i = 1; i < sites.Length; i++)
            {
                min = Vector3.Min(min, sites[i]);
                max = Vector3.Max(max, sites[i]);
            }
            Vector3 mid = (min + max) * 0.5f;
            Vector3 span = max - min;
            float radius = MathF.Max(span.X, MathF.Max(span.Y, span.Z)) * 16f;
            if (!(radius > 0) || float.IsInfinity(radius))
                throw new ArgumentException("Site bounds must be finite and nonzero", nameof(sites));
            int n = sites.Length;
            scratch.Points[n] = mid + radius * new Vector3(1, 1, 1);
            scratch.Points[n + 1] = mid + radius * new Vector3(-1, -1, 1);
            scratch.Points[n + 2] = mid + radius * new Vector3(-1, 1, -1);
            scratch.Points[n + 3] = mid + radius * new Vector3(1, -1, -1);

            output.Reserve(1);
            output.Tetrahedra[0] = new Tetrahedron(n, n + 1, n + 2, n + 3);
            if (!TryCircumcenter(scratch.Points, output.Tetrahedra[0], out output.Circumcenters[0]))
                throw new InvalidOperationException("Could not initialize the enclosing tetrahedron");
            output.Count = 1;

            for (int point = 0; point < n; point++)
            {
                scratch.ResetFaces(checked(output.Count * 4));
                int badCount = 0;
                for (int t = 0; t < output.Count;)
                {
                    Tetrahedron tet = output.Tetrahedra[t];
                    Vector3 center = output.Circumcenters[t];
                    float radiusSq = Vector3.DistanceSquared(center, scratch.Points[tet.A]);
                    float pointSq = Vector3.DistanceSquared(center, scratch.Points[point]);
                    if (pointSq <= radiusSq * (1f + 1e-6f))
                    {
                        scratch.ToggleFace(tet.A, tet.B, tet.C);
                        scratch.ToggleFace(tet.A, tet.B, tet.D);
                        scratch.ToggleFace(tet.A, tet.C, tet.D);
                        scratch.ToggleFace(tet.B, tet.C, tet.D);
                        output.RemoveSwapBack(t);
                        badCount++;
                    }
                    else t++;
                }
                if (badCount == 0)
                    throw new InvalidOperationException("No Delaunay cavity found; check site precision and degeneracy");

                for (int f = 0; f < scratch.FaceCount; f++)
                {
                    if (!scratch.FaceActive[f]) continue;
                    FaceKey face = scratch.Faces[f];
                    var tet = new Tetrahedron(face.A, face.B, face.C, point);
                    if (!TryCircumcenter(scratch.Points, tet, out Vector3 center))
                        throw new InvalidOperationException("Degenerate tetrahedron; sites must be in general position");
                    output.Reserve(output.Count + 1);
                    output.Tetrahedra[output.Count] = tet;
                    output.Circumcenters[output.Count++] = center;
                }
            }

            for (int t = 0; t < output.Count;)
            {
                Tetrahedron tet = output.Tetrahedra[t];
                if (tet.A >= n || tet.B >= n || tet.C >= n || tet.D >= n)
                    output.RemoveSwapBack(t);
                else t++;
            }
            return output;
        }

        // Solve the three perpendicular-bisector equations in double precision.
        private static bool TryCircumcenter(Vector3[] points, Tetrahedron t, out Vector3 center)
        {
            Vector3 p = points[t.A], q = points[t.B], r = points[t.C], s = points[t.D];
            double ux = q.X - p.X, uy = q.Y - p.Y, uz = q.Z - p.Z;
            double vx = r.X - p.X, vy = r.Y - p.Y, vz = r.Z - p.Z;
            double wx = s.X - p.X, wy = s.Y - p.Y, wz = s.Z - p.Z;
            double vwx = vy * wz - vz * wy, vwy = vz * wx - vx * wz, vwz = vx * wy - vy * wx;
            double wux = wy * uz - wz * uy, wuy = wz * ux - wx * uz, wuz = wx * uy - wy * ux;
            double uvx = uy * vz - uz * vy, uvy = uz * vx - ux * vz, uvz = ux * vy - uy * vx;
            double det2 = 2 * (ux * vwx + uy * vwy + uz * vwz);
            double scale = Math.Max(Math.Sqrt(ux * ux + uy * uy + uz * uz),
                Math.Max(Math.Sqrt(vx * vx + vy * vy + vz * vz), Math.Sqrt(wx * wx + wy * wy + wz * wz)));
            if (Math.Abs(det2) <= 1e-12 * scale * scale * scale)
            {
                center = default;
                return false;
            }
            double u2 = ux * ux + uy * uy + uz * uz;
            double v2 = vx * vx + vy * vy + vz * vz;
            double w2 = wx * wx + wy * wy + wz * wz;
            center = new Vector3(
                (float)(p.X + (u2 * vwx + v2 * wux + w2 * uvx) / det2),
                (float)(p.Y + (u2 * vwy + v2 * wuy + w2 * uvy) / det2),
                (float)(p.Z + (u2 * vwz + v2 * wuz + w2 * uvz) / det2));
            return !float.IsNaN(center.X) && !float.IsInfinity(center.X) &&
                !float.IsNaN(center.Y) && !float.IsInfinity(center.Y) &&
                !float.IsNaN(center.Z) && !float.IsInfinity(center.Z);
        }
    }
}
