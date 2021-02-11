# motif-mark
Python code to visualize motif locations on gene sequences


# Inputs:
Two files.
1) FASTA file containing genes to be examined. Either DNA or RNA reads in conventional FASTA format. 
2) File containing known motifs. One motif per line, no header information. Motifs must be in IUPAC format.

# Output:
A single .svg file (named the same as the input FASTA file) containing color-coded motifs shown to scale on a gene line with exons represented as boxes. Legend for motif color is provided. Motifs are shown at 65% opacity so overlapping motifs will be visible.

# Example Usage:
python ./motifmarker.py -f [FASTA file] -m [motif file]
