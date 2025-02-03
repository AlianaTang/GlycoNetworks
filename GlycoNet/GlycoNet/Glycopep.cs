// © 2024 Aliana Tang <alianatang17@gmail.com>

namespace GlycoNet
{
    internal class Glycopep
    {
        // Retention time, Charge, Precursor Mz, Precursor Mass, Peptide Mass, Glycan Mass, Intensity, [Hexes, HexNacs, Fucs, NeuAcs, NeuGcs]*
        public double rt;
        public double charge; //change to int later
        public double precursorMz;
        public double precursorMass;
        public double peptideMass;
        public double glycanMass;
        public double intensity;
        public GlycanComposition comp;
        public bool compChanged = false;
    }
}