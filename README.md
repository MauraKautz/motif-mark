# Motif Mark
Maura Kautz, March 2025

## Overview

Motif mark is an object oriented program used to visualize motifs along genes. The tool takes a fasta file and a motif text file, then generates a .png figure illustrating where motifs were identified in each fasta record. An example figure is shown below:

***EXAMPLE IMAGE***

Boxes on the genes represent exons, while the different colored lines represent different motifs.

## Details

The input motif file should contain one motif per line that is less than or equal to 10 bases. Motifs can be ambiguous when they use IUPAC notation. IUPAC notation includes nucleobases represented by A, T, C, and G, as well as eleven ambiguity characters representing every possible combination of the four DNA bases. Ambiguous motifs are handled within this tool.

The output of this tool is a .png file containing the cairo figure drawn with the fasta and motif specifications. The output file will follow the same naming convention as the fasta file. For example Figure_1.fasta will output Figure_1.png.

## Installation and Usage

To use motif mark, the user will need ```motif-mark.py```, a fasta file, a motif text file, and the appropriate conda environment.

Conda Environment:

This program is set to run in an evironment called my_pycairo. The only package installed and necessary in this environment is pycairo.

Using:

```motif-mark.py -f <absolute path to fasta file> -m <absolute path to motif file>```

**-f**: The fasta file containing gene records. The absolute path to the fasta file is required for the script to run.

**-m**: The motif text file containing one motif per line. The absolute path to the motif file is required for the script to run.
