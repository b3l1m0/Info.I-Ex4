# Info.I-Ex4
C++ scripts for Exercise 4 of the Informatik I Lecture of the first semester of my B.Sc. Chemical Engineering


fit.h, fit.cc & main.cpp: provided by the lecturer Prof. P. H. Hünenberger from ETHZ

ex4A: a program that extracts the amino-acid sequence of a protein from the corresponding PDB-file, and writes it into a corresponding sequence file listing the three-letter codes of the successive amino-acids

ex4B: modification of A so that it reads 2 PDB files, the first one the wild-type and the second one the mutant. Then the program determines directly the number of mutations and compares the corresponding three-dimensional structures and measures how much these differ by calculating the corresponding backbone atomic-positional root-mean-square
 deviations (RMSD)

ex4C: modification of B so that it can read 8 PDB files and compares each to each, sorts these comparisons in order of increasing RMSD 
