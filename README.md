Glycosylation is critical to many biological processes, such as protein folding and intercellular communication. Mass spectrometry is widely used for studying proteins and glycoproteins. For glycoproteins, successful data interpretation of mass spectrometry data using database search requires a complete glycan database — a condition that is often not satisfied.

GlycoNet builds improved N-glycan databases by constructing a sample-specific glycan database based on the mass spectrometry data itself, rather than relying solely on existing glycan databases. GlycoNet reads in mass spectrometry data and an initial (seed) glycan database. Then it analyzes the data to augment the initial glycan database with additional glycans. Rather than analyzing individual spectra in isolation, GlycoNet constructs a network (graph) of putative glycans where each node represents a glycan and each edge represents a monosaccharide.

More details can be found in the paper "Network method for building a sample-specific glycan database for N-linked glycosylation from MS/MS data."

## Usage summary

GlycoNet is a command line tool.
```
Usage: GlycoNet <Spectra file path> <Glycan database file path> [Mass tolerance] [Deltas]
  <Spectra file path> is mandatory. Must be in MGF format, or ConvertToMgf.bat must be configured to convert to MGF format
  <Glycan database file path> is mandatory. Must be a text file of glycan compositions
  [Mass tolerance] is optional. Default is: 0.02
  [Deltas] is optional. List of glycan deltas to consider, separated by commas. You may also specify a delta mass
    Default is: all building blocks found in the glycan database, along with the disaccharide Hex-HexNAc
    This program knows the following deltas: Hex,HexNAc,Fuc,NeuAc,NeuGc,Pent,Phospho,Na,Acetyl,Hex-HexNAc
  ```

## Usage example

The `GlycanTestData` folder has example data for testing. This data is from dataset [PXD010308](https://www.ebi.ac.uk/pride/archive/projects/PXD010308) in the PRIDE data repository and is discussed in [this paper](https://doi.org/10.3389/fimmu.2018.02645).

Example command for running GlycoNet on this test dataset:
```
GlycoNet "HF01_20171017_QXH_LXY_2017YFF0205400_GIgM_F1_R1.pd.mgf" "N-glycan 182 human no multiple fucose.txt" 0.015 Hex,HexNAc,Fuc,NeuAc,Hex-HexNAc
```
The console output shows the additional glycans found by GlycoNet:
```
Deltas: Hex,HexNAc,Fuc,NeuAc,Hex-HexNAc

Number of spectra: 23286
Number of spectra with glycopeptide oxonium ion and Y0, Y1, Y2, etc.: 474 (2.04%)

Additional glycans:
HexNAc(5)Hex(7)NeuAc(1)
HexNAc(5)Hex(8)NeuAc(1)
HexNAc(5)Hex(9)
HexNAc(6)Hex(10)
HexNAc(6)Hex(10)Fuc(1)
HexNAc(6)Hex(10)Fuc(2)
HexNAc(6)Hex(10)NeuAc(1)
HexNAc(6)Hex(11)
HexNAc(6)Hex(11)NeuAc(1)
HexNAc(6)Hex(12)
HexNAc(6)Hex(8)
HexNAc(6)Hex(8)Fuc(1)
HexNAc(6)Hex(9)Fuc(1)
HexNAc(6)Hex(9)NeuAc(1)
HexNAc(6)Hex(9)NeuAc(1)Fuc(1)
HexNAc(7)Hex(11)
HexNAc(7)Hex(13)
```

GlycoNet also generates auxiliary files with more details:
* The additional glycans as a text file
* Details about the additional glycans
* Summary of all the spectra with putative glycans based on oxonium ions and Y<sub>0</sub>, Y<sub>1</sub>, Y<sub>2</sub>, etc.
* A `spectra` subfolder containing all the spectra with oxonium ions colored purple and Y<sub>0</sub>, Y<sub>1</sub>, Y<sub>2</sub>, etc. colored red (use [gnuplot](http://www.gnuplot.info/) to view)
* A `graphs` subfolder containing all the graphs, one for each peptide, in [GML](https://en.wikipedia.org/wiki/Graph_Modelling_Language) format (use [yEd](https://www.yworks.com/products/yed) with hierarchical layout, for example, to view)
