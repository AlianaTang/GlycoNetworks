// © 2024 Aliana Tang <alianatang17@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoNet
{
    internal static class Constants
    {
        public static double tolerance = 0.02;

        public static Dictionary<String, double> glycan_masses = new Dictionary<String, double>
        {
            ["Hex"] = 162.052823,
            ["HexNAc"] = 203.079373,
            ["Fuc"] = 146.057909,
            ["NeuAc"] = 291.095417,
            ["NeuGc"] = 307.090331,
            ["Pent"] = 150.05282342 - 18.010565,
            ["Phospho"] = 79.966331,
            ["Na"] = 21.981943,
            ["Acetyl"] = 42.010565,
            ["Hex-HexNAc"] = 162.052823 + 203.079373  // Special case of the most common disaccharide delta
        };

        public const double isotopeSpacing = 1.003;
        public const double adduct = 1.007825035 - 0.00054858;
    }
}
