// © 2024 Aliana Tang <alianatang17@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace GlycoNet
{
    internal static class Utilities
    {
        public static List<MS2> readMgf(string inputfilename)
        {
            string[] lines = File.ReadAllLines(inputfilename);
            var spectra = new List<MS2>();

            int count = -1;

            for (int i = 0; i < lines.Length; i++)
            {
                if (lines[i] == "BEGIN IONS")
                {
                    spectra.Add(new MS2());
                    count++;
                }
                else if (lines[i].Split('=')[0] == "TITLE")
                {

                }
                else if (lines[i].Split('=')[0] == "RTINSECONDS")
                {
                    spectra[count].Rt = Convert.ToDouble(lines[i].Split('=')[1]);
                }
                else if (lines[i].Split('=')[0] == "PEPMASS")
                {
                    spectra[count].PrecursorMz = Convert.ToDouble(lines[i].Split('=')[1].Split(' ')[0]);
                }
                else if (lines[i].Split('=')[0] == "CHARGE")
                {
                    spectra[count].Charge = Int32.Parse(lines[i].Split('=')[1].Split('+')[0]);
                }
                else if (lines[i] == "END IONS")
                {

                }
                else if (Double.TryParse(lines[i].Split(' ')[0], out double mass) && Double.TryParse(lines[i].Split(' ')[1], out double intensity))
                {
                    spectra[count].AddData(mass, intensity);
                }
            }
            return spectra;
        }

        public static List<Glycopep> readCsv(string inputfilename)
        {
            string[] lines = File.ReadAllLines(inputfilename);
            List<Glycopep> infoTable = new List<Glycopep>(); //better name than infoTable
            CompareGlycan compareGlyc = new CompareGlycan();
            for (int i = 13; i < lines.Length; i++)
            {
                Glycopep glycopeptide = new Glycopep();
                glycopeptide.precursorMass = Convert.ToDouble(lines[i].Split("\t")[8]) - Constants.adduct;
                glycopeptide.comp = GlycanComposition.toComp(lines[i].Split("\t")[3]);

                if (glycopeptide.comp != null)
                {
                    glycopeptide.glycanMass = compareGlyc.calcMass(glycopeptide.comp);
                    glycopeptide.peptideMass = glycopeptide.precursorMass - glycopeptide.glycanMass;
                }
                infoTable.Add(glycopeptide);
            }
            Console.WriteLine(infoTable.Count);
            return infoTable;
        }

        public static List<GlycanComposition> readGlycanDatabase(string inputfilename = "")
        {
            string[] lines = File.ReadAllLines(inputfilename);
            List<GlycanComposition> comp_labels = new List<GlycanComposition>();

            for (int i = 0; i < lines.Length; i++)
            {
                string composition = lines[i].Split('%')[0];
                comp_labels.Add(GlycanComposition.toComp(composition));
                //Console.Write("comp:");
                //Console.WriteLine(GlycanComposition.toComp(composition).toString());       
            }
            return comp_labels;
        }

        public static double mzToMass(double mz, int charge_num)
        {
            return mz * charge_num - charge_num * Constants.adduct;
        }

        public static double massToMz(double mass, int charge_num)
        {
            return (mass + charge_num * Constants.adduct) / charge_num;
        }

        public static Point monoisotope(double peak, double charge, List<Point> data)
        {
            int i = 1;
            Point current_peak = Array.Find(data.ToArray(), ele => Math.Abs(ele.mass - (peak - (i - 1) * Constants.isotopeSpacing / charge)) <= Constants.tolerance);
            if (current_peak == null)
            {
                Point dun = new Point();
                dun.mass = -1;
                return dun;
            }
            while (true)
            {
                if (Array.Exists(data.ToArray(), ele => ele.mass >= peak - i * Constants.isotopeSpacing / charge - Constants.tolerance && ele.mass <= peak - i * Constants.isotopeSpacing / charge + Constants.tolerance))
                {
                    Point next_peak = Array.Find(data.ToArray(), ele => ele.mass >= peak - i * Constants.isotopeSpacing / charge - Constants.tolerance && ele.mass <= peak - i * Constants.isotopeSpacing / charge + Constants.tolerance);
                    if (current_peak.intensity / next_peak.intensity <= 10)
                    {
                        current_peak = next_peak;
                        i++;
                    }
                    else
                    {
                        return current_peak;
                    }
                }
                else
                {
                    return current_peak;
                }
            }
        }
    }
}
