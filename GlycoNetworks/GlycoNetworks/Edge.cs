// © 2024, 2025 Aliana Tang <alianatang17@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoNetworks
{
    internal class Edge
    {
        public Glycopep source;
        public Glycopep target;
        public string glycDiff;
        public int increasing;
    }
}
