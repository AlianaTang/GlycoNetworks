// © 2024 Aliana Tang <alianatang17@gmail.com>

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GlycoNet
{
    internal class CompareGlycan
    {
        public List<string> glycanComps = new List<string>();
        public List<double> knownGlycans = new List<double>();

        //public Dictionary<double, List<double>> graph1 = new Dictionary<double, List<double>>(); // keys are the nodes/glycan masses; values are their neighbors/more glycan masses.
        public List<Glycopep> visited = new List<Glycopep>();

        //public Dictionary<double, List<double>> graph2 = new Dictionary<double, List<double>>();
        //public List<double> visited2 = new List<double>();

        //public Dictionary<Glycopep, List<Edge>> graph3 = new Dictionary<Glycopep, List<Edge>>();

        public List<String> compList = new List<String>();
        public double calcHex(string composition)
        {
            double num_hex = 0;
            if (composition.Contains("Hex("))
            {
                num_hex = Convert.ToDouble(composition.Split("Hex(")[1].Split(")")[0]);
            }
            return num_hex;
        }

        public double calcHexNac(string composition)
        {
            double num_hexNac = 0;
            if (composition.Contains("HexNAc("))
            {
                num_hexNac = Convert.ToDouble(composition.Split("HexNAc(")[1].Split(")")[0]);
            }
            return num_hexNac;
        }

        public double calcFuc(string composition)
        {
            double num_fuc = 0;
            if (composition.Contains("Fuc("))
            {
                num_fuc = Convert.ToDouble(composition.Split("Fuc(")[1].Split(")")[0]);

            }
            return num_fuc;
        }
        public double calcNeuAc(string composition)
        {
            double num_neuAc = 0;
            if (composition.Contains("NeuAc("))
            {
                num_neuAc = Convert.ToDouble(composition.Split("NeuAc(")[1].Split(")")[0]);
            }
            return num_neuAc;
        }

        public double calcNeuGc(string composition)
        {
            double num_neuGc = 0;
            if (composition.Contains("NeuGc("))
            {
                num_neuGc = Convert.ToDouble(composition.Split("NeuGc(")[1].Split(")")[0]);

            }
            return num_neuGc;
        }
        public double calcPentose(string composition)
        {
            double num_pentose = 0;
            if (composition.Contains("Pent("))
            {
                num_pentose = Convert.ToDouble(composition.Split("Pent(")[1].Split(")")[0]);

            }
            return num_pentose;
        }

        public double calcPhospho(string composition)
        {
            double num_phospho = 0;
            if (composition.Contains("Phospho("))
            {
                num_phospho = Convert.ToDouble(composition.Split("Phospho(")[1].Split(")")[0]);

            }
            return num_phospho;
        }

        public double calcMass(GlycanComposition composition)
        {
            double totalMass = 0;
            for (int i = 0; i < composition.Keys.Count; i++)
            {
                String key = composition.Keys.ToList()[i];
                int number = composition[key];
                double mass = Constants.glycan_masses[key];
                totalMass += number * mass;
            }
            
            return totalMass;
        }


        //graph has key may be discarded in favor of a built-in function that works better (but i can't find one)
        //or it may be moved to utilities? depends
        public bool graphHasKey(Dictionary<Glycopep, List<Edge>> graph, double keyGlycMass, double keyPepMass)
        {
            if (Array.Exists(graph.Keys.ToArray(), ele => Math.Abs(ele.glycanMass - keyGlycMass) <= Constants.tolerance && Math.Abs(ele.peptideMass - keyPepMass) <= Constants.tolerance))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public int graphFindKey(Dictionary<Glycopep, List<Edge>> graph, double keyGlycMass, double keyPepMass)
        {
            return (Array.FindIndex(graph.Keys.ToArray(), ele => Math.Abs(ele.glycanMass - keyGlycMass) <= Constants.tolerance && Math.Abs(ele.peptideMass - keyPepMass) <= Constants.tolerance));
        }

        public string writeEdge(int biggerIndex, int smallerIndex, double glycan_mass, string labelname, List<Glycopep> glycopeptides, Dictionary<Glycopep, List<Edge>> graph)
        {   
            string edge_code = "";
            double bigger = glycopeptides[biggerIndex].glycanMass;
            double smaller = glycopeptides[smallerIndex].glycanMass;
            double diff = bigger - smaller;
            
            if ((glycan_mass - Constants.tolerance < diff && diff < glycan_mass + Constants.tolerance))
            {
                /*
                if (!graph1.ContainsKey(smaller))
                {
                    List<double> temp = new List<double>();
                    graph1.Add(smaller, temp);
                }
                if (!graph1.ContainsKey(bigger))
                {
                    List<double> temp = new List<double>();
                    graph1.Add(bigger, temp);
                }
                */
                edge_code = "  edge\n  [ source  " + smallerIndex.ToString() + "\n    target  " + biggerIndex.ToString() + "\n    LabelGraphics" + 
                    "\n    [ text \"" + labelname + "\" ]\n  ]\n";
                /*
                List<double> temp1 = graph1[smaller];
                temp1.Add(Convert.ToDouble(bigger.ToString("F4")));
                graph1[smaller] = temp1;
                List<double> temp2 = graph1[bigger];
                temp2.Add(Convert.ToDouble(smaller.ToString("F4")));
                graph1[bigger] = temp2;
                */

                //

                //Console.WriteLine(graphFindKey(graph3, bigger, pep2));


                Glycopep glycopep1 = glycopeptides[biggerIndex];
                Glycopep glycopep2 = glycopeptides[smallerIndex];
                //Console.WriteLine(glycopep1.glycanMass);
                //Console.WriteLine(bigger);
                //Console.WriteLine();

                //Glycopep glycopep1 = graph3.Keys.ToArray()[graphFindKey(graph3, glycMasses[smallerIndex], pep1)];
                //Glycopep glycopep2 = graph3.Keys.ToArray()[graphFindKey(graph3, glycMasses[biggerIndex], pep2)];
                
                Edge edge1 = new Edge();
                edge1.source = glycopep1;
                edge1.target = glycopep2;
                edge1.glycDiff = labelname;
                edge1.increasing = -1;

                Edge edge2 = new Edge();
                edge2.source = glycopep2;
                edge2.target = glycopep1;
                edge2.glycDiff = labelname;
                edge2.increasing = 1;

                graph[glycopep1].Add(edge1);
                graph[glycopep2].Add(edge2);
            }
            /*
            else
            {
                if (!graph1.ContainsKey(smaller))
                {
                    List<double> temp = new List<double>();
                    graph1.Add(smaller, temp);
                }
                if (!graph1.ContainsKey(bigger))
                {
                    List<double> temp = new List<double>();
                    graph1.Add(bigger, temp);
                }
            }
            */
            return edge_code;
        }

        public void dfs(Dictionary<double, List<double>> graph, double node, double tolerance, List<double> visited)
        {
            visited.Add(node);
            foreach (double neighbor in graph[node])
            {
                if (!visited.Contains(neighbor))
                {
                    dfs(graph, neighbor, tolerance, visited);
                }
            }            
        }

        public void dfs2(Dictionary<Glycopep, List<Edge>> graph, Glycopep node, double tolerance, List<Glycopep> visited)
        {
            visited.Add(node);
            //Console.WriteLine(node.glycanMass);
            foreach (Edge neighbor in graph[node])
            {
                if (!visited.Contains(neighbor.target))
                {
                    //Console.WriteLine("Hello World");
                    dfs2(graph, neighbor.target, tolerance, visited);
                }
            }

        }
    }   
}
