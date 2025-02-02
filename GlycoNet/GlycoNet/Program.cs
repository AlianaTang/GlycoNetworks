// © 2024 Aliana Tang <alianatang17@gmail.com>

using GlycoNet;
using System.Collections.Generic;
using System.Diagnostics;

const double defaultMassTolerance = 0.02;

if (args.Length < 2 || args.Length > 4)
{
    Console.WriteLine("Usage: GlycoNet <Spectra file path> <Glycan database file path> [Mass tolerance] [Deltas]");
    Console.WriteLine("  <Spectra file path> is mandatory. Must be in MGF format, or ConvertToMgf.bat must be configured to convert to MGF format");
    Console.WriteLine("  <Glycan database file path> is mandatory. Must be a text file of glycan compositions");
    Console.WriteLine("  [Mass tolerance] is optional. Default is: " + defaultMassTolerance);
    Console.WriteLine("  [Deltas] is optional. List of glycan deltas to consider, separated by commas");
    Console.WriteLine("    Default is: all building blocks found in the glycan database, along with the disaccharide Hex-HexNAc");
    Console.WriteLine("    This program knows the following deltas: " + String.Join(",", Constants.glycan_masses.Keys));
    return;
}

string spectraFilePath = Path.GetFullPath(args[0]);

if (Path.GetExtension(spectraFilePath).ToLower() != ".mgf")
{
    var processStartInfo = new ProcessStartInfo(Path.Combine(AppContext.BaseDirectory, "ConvertToMgf.bat"), spectraFilePath);
    processStartInfo.UseShellExecute = false;
    processStartInfo.CreateNoWindow = true;
    processStartInfo.WorkingDirectory = Path.GetDirectoryName(spectraFilePath);
    Process? process = Process.Start(processStartInfo);
    process.WaitForExit();

    spectraFilePath = Path.Combine(Path.GetDirectoryName(spectraFilePath), Path.GetFileNameWithoutExtension(spectraFilePath) + ".mgf");
    if (!File.Exists(spectraFilePath))
    {
        Console.WriteLine("Error: Unable to convert " + args[0] + " to MGF format");
        return;
    }
}

string glycanDatabaseFilePath = args[1];

if (args.Length >= 3)
{
    if (!Double.TryParse(args[2], out Constants.tolerance))
    {
        Console.WriteLine("Error: Could not parse mass tolerance");
        return;
    }
}
else
{
    Constants.tolerance = defaultMassTolerance;
}

string[]? deltas = null;
if (args.Length >= 4)
{
    deltas = args[3].Split(",");
    foreach (string delta in deltas)
    {
        if (!Constants.glycan_masses.ContainsKey(delta))
        {
            Console.WriteLine("Error: Invalid delta");
            return;
        }
    }
}

string spectraFileName = Path.GetFileName(spectraFilePath);
string spectraDirectory = Path.GetDirectoryName(spectraFilePath);

Directory.CreateDirectory(Path.Combine(spectraDirectory, ("spectra - " + spectraFileName)));
Directory.CreateDirectory(Path.Combine(spectraDirectory, ("graphs - " + spectraFileName)));

List<GlycanComposition> compositions = Utilities.readGlycanDatabase(glycanDatabaseFilePath);
if (deltas == null)
{
    var deltasHashSet = new HashSet<string>();
    foreach (GlycanComposition composition in compositions)
    {
        foreach (string delta in composition.Keys)
        {
            deltasHashSet.Add(delta);
        }
    }
    if (deltasHashSet.Contains("Hex") && deltasHashSet.Contains("HexNAc"))
    {
        deltasHashSet.Add("Hex-HexNAc");
    }
    deltas = deltasHashSet.ToArray();
}
Console.WriteLine("Deltas: " + String.Join(",", deltas));
Console.WriteLine();

List<MS2> spectra = Utilities.readMgfAndFindGlycopeptides(Path.Combine(spectraDirectory, spectraFileName));
var glycanMasses = new List<double>();
var glycopepList = new List<Glycopep>();  // Retention time, Charge, Precursor Mz, Precursor Mass, Peptide Mass, Glycan Mass, Intensity, [Hexes, HexNacs, Fucs, NeuAcs, NeuGcs]*
var commonGlycans = new List<GlycInfo>();

var mapGlycantoPepList = new Dictionary<double, List<Glycopep>>(); //glycmass, list of pepmasses
var mapPeptoGlycanList = new Dictionary<double, List<Glycopep>>(); //used to be double, PepInfo2

var uniqueGlycansGlycopepList = new List<Glycopep>();

var graph2 = new Dictionary<Glycopep, List<Edge>>();
var graph3 = new Dictionary<Glycopep, List<Edge>>();

using (var gnuplotConvenienceFile = new StreamWriter(Path.Combine(spectraDirectory, ("spectra - " + spectraFileName + @"\all.txt"))))
{
    using (var spectraSummaryCsvFile = new StreamWriter(Path.Combine(spectraDirectory, (spectraFileName + " - spectra summary.csv"))))
    {
        spectraSummaryCsvFile.WriteLine("Retention Time (min),Charge,Precursor Mz,Precursor Mass,Peptide Mass,Glycan Mass,Intensity");

        for (int i = 0; i < spectra.Count; i++)
        {
            Point? pepMz = spectra[i].findGlycopep();  // Look for pep, pep+HexNAc, pep+2HexNAc, etc. signature
            if (pepMz != null && pepMz.mass > 0)
            {
                double pepMass = Utilities.mzToMass(pepMz.mass, spectra[i].returnBarePepCharge());
                double glycMass = Utilities.mzToMass(spectra[i].PrecursorMz, spectra[i].Charge) - pepMass;
                glycanMasses.Add(glycMass);

                spectraSummaryCsvFile.WriteLine((spectra[i].Rt / 60).ToString("F4") + "," + spectra[i].Charge + "," + spectra[i].PrecursorMz + ","
                    + Utilities.mzToMass(spectra[i].PrecursorMz, spectra[i].Charge) + "," + pepMass + ","
                    + glycMass + "," + pepMz.intensity);

                Glycopep glycopep = new Glycopep();
                glycopep.rt = Convert.ToDouble((spectra[i].Rt / 60).ToString("F4"));
                glycopep.charge = spectra[i].Charge;
                glycopep.precursorMz = spectra[i].PrecursorMz;
                glycopep.precursorMass = Utilities.mzToMass(spectra[i].PrecursorMz, spectra[i].Charge);
                glycopep.peptideMass = pepMass;
                glycopep.glycanMass = glycMass;
                glycopep.intensity = pepMz.intensity;
                glycopepList.Add(glycopep);

                string gnuplotSpectrumFileName = Path.Combine(spectraDirectory, @"spectra - " + spectraFileName + @"\rt=" + (spectra[i].Rt / 60).ToString("F4") + "_mz=" + spectra[i].PrecursorMz.ToString("F4") + ".txt");
                gnuplotConvenienceFile.WriteLine("load '" + Path.GetFileName(gnuplotSpectrumFileName) + "'; pause -1");

                using (var gnuplotSpectrumFile = new StreamWriter(gnuplotSpectrumFileName))
                {
                    string plotTitle = @"Retention time: " + spectra[i].Rt / 60 + @"\nCharge: " + spectra[i].Charge + @"\nPrecursor m/z: " + spectra[i].PrecursorMz;
                    gnuplotSpectrumFile.WriteLine("unset key\nset title \"" + plotTitle + "\"\nplot '-' with impulse linewidth 2 linetype - 1, '-' with impulse linewidth 2 linetype 7, '-' with impulse linewidth 2 linetype 17");

                    // Draw spectrum in black in gnuplot
                    for (int j = 0; j < spectra[i].findLen(); j++)
                    {
                        gnuplotSpectrumFile.WriteLine((spectra[i].GetMass(j)).ToString("F4") + " " + (spectra[i].GetIntensity(j)).ToString("F4"));
                    }

                    spectra[i].sortPeaksByIntensity();

                    // Make the pep, pep+HexNAc, pep+2HexNAc, pep+2HexNAc+Hex peaks red in gnuplot
                    Point plusHexnac = spectra[i].returnHexnac();
                    Point plus2Hexnac = spectra[i].return2Hexnac();
                    Point plusHex2Hexnac = spectra[i].returnHex2Hexnac();
                    gnuplotSpectrumFile.WriteLine("\n\ne\n\n" + pepMz.mass.ToString("F4") + " " + pepMz.intensity);
                    gnuplotSpectrumFile.WriteLine(plusHexnac.mass + " " + plusHexnac.intensity);
                    gnuplotSpectrumFile.WriteLine(plus2Hexnac.mass + " " + plus2Hexnac.intensity);
                    gnuplotSpectrumFile.WriteLine(plusHex2Hexnac.mass + " " + plusHex2Hexnac.intensity);
                    gnuplotSpectrumFile.WriteLine("# " + (spectra[i].Rt / 60).ToString("F4"));

                    // Make the 204 and 366 peaks purple in gnuplot
                    gnuplotSpectrumFile.WriteLine("\n\ne\n\n" + Constants.glycan_masses["HexNAc"] + Constants.adduct + " " + spectra[i].corrIntensity(Constants.glycan_masses["HexNAc"] + Constants.adduct));
                    gnuplotSpectrumFile.WriteLine(Constants.glycan_masses["HexNAc"] + Constants.glycan_masses["Hex"] + Constants.adduct + " " + spectra[i].corrIntensity(Constants.glycan_masses["HexNAc"] + Constants.glycan_masses["Hex"] + Constants.adduct));
                    gnuplotSpectrumFile.WriteLine("\n\ne\n\n");
                }
            }
        }
    }
}

CompareGlycan glycans = new CompareGlycan();

for (int i = 0; i < glycopepList.Count; i++)
{
    for (int j = 0; j < compositions.Count; j++)
    {
        if (Math.Abs(glycopepList[i].glycanMass - glycans.calcMass(compositions[j])) <= Constants.tolerance)
        {
            glycopepList[i].comp = compositions[j];
        }
    }
}

for (int i = 0; i < glycopepList.Count; i++)
{
    if (Array.Exists(mapGlycantoPepList.Keys.ToArray(), ele => Math.Abs(ele - glycopepList[i].glycanMass) <= Constants.tolerance))
    {
        double keyVal = Array.Find(mapGlycantoPepList.Keys.ToArray(), ele => Math.Abs(ele - glycopepList[i].glycanMass) <= Constants.tolerance);
        List<Glycopep> list = mapGlycantoPepList[keyVal];
        list.Add(glycopepList[i]);
        mapGlycantoPepList[keyVal] = list;
    } else
    {
        List<Glycopep> list = new List<Glycopep>(); 
        list.Add(glycopepList[i]);
        mapGlycantoPepList.Add(glycopepList[i].glycanMass, list);
    }
 
    if (Array.Exists(mapPeptoGlycanList.Keys.ToArray(), ele => Math.Abs(ele - glycopepList[i].peptideMass) <= Constants.tolerance))
    {
        double keyVal = Array.Find(mapPeptoGlycanList.Keys.ToArray(), ele => Math.Abs(ele - glycopepList[i].peptideMass) <= Constants.tolerance);
        List<Glycopep> list = mapPeptoGlycanList[keyVal];
        list.Add(glycopepList[i]);
        mapPeptoGlycanList[keyVal] = list;       
    }
    else
    {
        List<Glycopep> list = new List<Glycopep>();
        list.Add(glycopepList[i]);
        mapPeptoGlycanList.Add(glycopepList[i].peptideMass, list);
    }

    GlycInfo temp = new GlycInfo();
    temp.glycMass = glycopepList[i].glycanMass;
    temp.occurrences = 1;
    temp.rtOcc = glycopepList[i].rt;
    temp.comp = glycopepList[i].comp;
    temp.clusterSize = 0;
    commonGlycans.Add(temp);

    for (int j = i + 1; j < glycopepList.Count; j++)
    {
        if ((Math.Abs(glycopepList[i].glycanMass - glycopepList[j].glycanMass) <= Constants.tolerance) && (Math.Abs(glycopepList[i].peptideMass - glycopepList[j].peptideMass) <= Constants.tolerance))
        {
            commonGlycans[i].occurrences++;
            if (glycopepList[j].comp != null)
            {
                commonGlycans[i].comp = glycopepList[j].comp;
            }
            glycopepList.RemoveAt(j);
            j--;
        }
    }
}

for (int i = 0; i < mapGlycantoPepList.Count; i++)
{
    for (int j = 0; j < glycopepList.Count; j++)
    {
        if (mapGlycantoPepList.Keys.ToArray()[i] - Constants.tolerance <= glycopepList[j].glycanMass && mapGlycantoPepList.Keys.ToArray()[i] + Constants.tolerance >= glycopepList[j].glycanMass && glycopepList[j].comp != null)
        {
            double keyVal = mapGlycantoPepList.Keys.ToArray()[i];
            //mapGlycantoPepList[keyVal].comp = infoTable[j].comp;
        }
    }
}

for (int i = 0; i < mapPeptoGlycanList.Count; i++)
{
    double keyVal = mapPeptoGlycanList.Keys.ToArray()[i];
    for (int k = 0; k < mapPeptoGlycanList[keyVal].Count; k++)
    {
        //mapPeptoGlycanList[keyVal].comps.Add("");
        for (int j = 0; j < glycopepList.Count; j++)
        {
            if (Math.Abs(mapPeptoGlycanList[keyVal][k].glycanMass - glycopepList[j].glycanMass) <= Constants.tolerance && glycopepList[j].comp != null)
            {
                mapPeptoGlycanList[keyVal][k].comp = (glycopepList[j].comp);
            }
        }
    }
}

for (int i = 0; i < mapGlycantoPepList.Count; i++)
{
    double keyVal = mapGlycantoPepList.Keys.ToArray()[i];
    uniqueGlycansGlycopepList.Add(mapGlycantoPepList[keyVal][0]);
}

using (var glycanGraphFile = new StreamWriter(Path.Combine(spectraDirectory, @"glycan graph - " + spectraFileName + ".gml")))
{
    glycanGraphFile.WriteLine("Creator \"yFiles\"\nVersion 2.2\ngraph\n[ hierarchic  1\n  directed  1  ");

    for (int i = 0; i < uniqueGlycansGlycopepList.Count; i++)
    {
        string label = "";
        if (uniqueGlycansGlycopepList[i].comp != null)
        {
            label = uniqueGlycansGlycopepList[i].comp.toString() + ": ";
        }
        glycanGraphFile.WriteLine("  node");
        glycanGraphFile.WriteLine("  [ id  " + i.ToString());
        glycanGraphFile.WriteLine("    LabelGraphics");
        glycanGraphFile.WriteLine("    [ text \"" + label + uniqueGlycansGlycopepList[i].glycanMass + "\" ]");

        if (label != "")
        {
            glycanGraphFile.WriteLine("    graphics");
            glycanGraphFile.WriteLine("    [ fill	\"#4287F5\" ]");
        }
        glycanGraphFile.WriteLine("  ]");

        graph2.Add(uniqueGlycansGlycopepList[i], new List<Edge>());
    }

    for (int i = 0; i < uniqueGlycansGlycopepList.Count; i++)
    {
        for (int j = i + 1; j < uniqueGlycansGlycopepList.Count; j++)
        {
            int biggerIndex, smallerIndex;

            if (uniqueGlycansGlycopepList[i].glycanMass > uniqueGlycansGlycopepList[j].glycanMass)
            {
                biggerIndex = i;
                smallerIndex = j;
            }
            else
            {
                biggerIndex = j;
                smallerIndex = i;
            }

            foreach (string delta in deltas)
            {
                glycanGraphFile.Write(glycans.writeEdge(biggerIndex, smallerIndex, Constants.glycan_masses[delta], delta, uniqueGlycansGlycopepList, graph2));
            }
        }
    }

    glycanGraphFile.WriteLine("]");
}

for (int k = 0; k < mapPeptoGlycanList.Count; k++)
{
    double keyVal = mapPeptoGlycanList.Keys.ToArray()[k];
    using (var glycanGraphSpecificPeptideFile = new StreamWriter(Path.Combine(spectraDirectory, @"graphs - " + spectraFileName + @"\glycan graph pep - " + keyVal + " - " + spectraFileName + ".gml")))
    {
        int len = mapPeptoGlycanList[keyVal].Count;

        glycanGraphSpecificPeptideFile.WriteLine("Creator \"yFiles\"\nVersion 2.2\ngraph\n[ hierarchic  1\n  directed  1  ");

        for (int i = 0; i < len; i++)
        {
            string label = "";
            if (mapPeptoGlycanList[keyVal][i].comp != null)
            {
                label = mapPeptoGlycanList[keyVal][i].comp.toString() + ": ";
            }
            glycanGraphSpecificPeptideFile.WriteLine("  node");
            glycanGraphSpecificPeptideFile.WriteLine("  [ id  " + i.ToString());
            glycanGraphSpecificPeptideFile.WriteLine("    LabelGraphics");
            glycanGraphSpecificPeptideFile.WriteLine("    [ text \"" + label + mapPeptoGlycanList[keyVal][i].glycanMass + "\nrt: " + mapPeptoGlycanList[keyVal][i].rt + "\" ]");

            if (mapPeptoGlycanList[keyVal][i].comp != null)//(!label.Equals(""))
            {
                glycanGraphSpecificPeptideFile.WriteLine("    graphics");
                glycanGraphSpecificPeptideFile.WriteLine("    [ fill	\"#4287F5\" ]");
            }
            glycanGraphSpecificPeptideFile.WriteLine("  ]");

            //Glycopep key = new Glycopep();
            //key.glycanMass = Convert.ToDouble(mapPeptoGlycanList[keyVal][i].glycanMass);
            //key.peptideMass = keyVal;
            //key.comp = GlycanComposition.toComp(label);

            graph3.Add(mapPeptoGlycanList[keyVal][i], new List<Edge>());
        }
        
        for (int i = 0; i < len; i++)
        {
            for (int j = i + 1; j < len; j++)
            {
                int biggerIndex, smallerIndex;

                if (mapPeptoGlycanList[keyVal][i].glycanMass > mapPeptoGlycanList[keyVal][j].glycanMass)
                {
                    biggerIndex = i;
                    smallerIndex = j;
                }
                else
                {
                    biggerIndex = j;
                    smallerIndex = i;
                }

                List<Glycopep> list = mapPeptoGlycanList[keyVal];

                foreach (string delta in deltas)
                {
                    glycanGraphSpecificPeptideFile.Write(glycans.writeEdge(biggerIndex, smallerIndex, Constants.glycan_masses[delta], delta, list, graph3));
                }
            }
        }
    }
}

var uniqueAdditionalGlycans = new HashSet<string>();
using (var additionalGlycansFile = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName + " - additional glycans - details.csv")))
{
    additionalGlycansFile.WriteLine(@"Glycan composition,Glycan mass,Cluster size,Peptide mass,Retention time");
    bool continueIterating = true;
    while (continueIterating)
    {
        continueIterating = false;
        for (int i = 0; i < glycopepList.Count; i++)
        {
            glycans.visited.Clear();
            glycans.dfs2(graph3, glycopepList[i], Constants.tolerance, glycans.visited);
            if (glycopepList[i].comp != null && glycans.visited.Count >= 3)
            {
                foreach (Edge edge in graph3[glycopepList[i]])
                {
                    if (edge.target.comp == null)
                    {
                        continueIterating = true;
                        edge.target.comp = GlycanComposition.toComp(edge.source.comp.toString());
                        if (edge.glycDiff == "Hex-HexNAc")
                        {
                            // Hex-HexNAc disaccharide delta special case
                            if (!edge.target.comp.ContainsKey("Hex")) edge.target.comp.Add("Hex", 0);
                            if (!edge.target.comp.ContainsKey("HexNAc")) edge.target.comp.Add("HexNAc", 0);
                            edge.target.comp["Hex"] += edge.increasing;
                            edge.target.comp["HexNAc"] += edge.increasing;
                        }
                        else
                        {
                            // Monosaccharide delta
                            if (!edge.target.comp.ContainsKey(edge.glycDiff)) edge.target.comp.Add(edge.glycDiff, 0);
                            edge.target.comp[edge.glycDiff] += edge.increasing;
                        }
                        if (edge.target.comp.isValid() && !GlycanComposition.Contains(compositions, edge.target.comp))
                        {
                            string new_comp = edge.target.comp.toString();
                            additionalGlycansFile.WriteLine(new_comp + "," + edge.target.glycanMass + "," + glycans.visited.Count + "," + edge.target.peptideMass + "," + edge.target.rt);
                            uniqueAdditionalGlycans.Add(new_comp);
                        }
                    }
                }
            }
        }
    }
}

List<string> sortedUniqueAdditionalGlycans = uniqueAdditionalGlycans.ToList();
sortedUniqueAdditionalGlycans.Sort();
Console.WriteLine("Additional glycans:");
using (var uniqueAdditionalGlycansFile = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName + " - additional glycans - summary.txt")))
{
    foreach (string g in sortedUniqueAdditionalGlycans)
    {
        uniqueAdditionalGlycansFile.WriteLine(g);
        Console.WriteLine(g);
    }
}

/*
using (var file6 = new StreamWriter(Path.Combine(spectraDirectory, spectraFileName) + " - mostCommonGlycans.csv"))
{
    file6.WriteLine("Glycan Mass, # Occurences, Retention Times of Occurences, Composition, Cluster Size");
    for (int i = 0; i < commonGlycans.Count; i++)
    {
        file6.WriteLine(commonGlycans[i].glycMass + ", " + commonGlycans[i].occurrences + ", " + commonGlycans[i].rtOcc + ", " + commonGlycans[i].comp + ", " + commonGlycans[i].clusterSize);
    }
}
*/
