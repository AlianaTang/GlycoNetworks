// © 2024 Aliana Tang <alianatang17@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoNet
{
    internal class ComparerByIntensity : IComparer<Point>
    {
        public int Compare(Point A, Point B)
        {
            if (A.intensity == B.intensity)
            {
                return 0;
            }
            else if (A.intensity < B.intensity)
            {
                return 1;
            }
            return -1;
        }
    }
}
