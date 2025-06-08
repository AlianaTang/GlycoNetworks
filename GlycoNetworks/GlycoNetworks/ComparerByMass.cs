// © 2024, 2025 Aliana Tang <alianatang17@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoNetworks
{
    internal class ComparerByMass : IComparer<Point>
    {
        public int Compare(Point A, Point B)
        {
            if (A.mass == B.mass)
            {
                return 0;
            }
            else if (A.mass < B.mass)
            {
                return 1;
            }
            return -1;
        }
    }
}
