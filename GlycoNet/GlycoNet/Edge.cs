// © 2024 Aliana Tang <alianatang17@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoNet
{
    internal class Edge
    {
        public Glycopep source;
        public Glycopep target;
        public string glycDiff;
        public int increasing;
    }
}
