This program computes the double homology of the moment angle complex associated to a simplicial complex


It takes as an input the maximal simplices that define a simplicial complex K and outputs the double homology
of ZK with coefficients in Z/2. 

The output separates HH by homological degree and for every bidegree it outputs the rank and shows a basis of it.
Note: Due to the nature of the program, it can be easily modified to also print a basis for the double homology.

By "Homological degree" I mean the degree of the reduced homology in the hochster decomposition
namely, for the bidegree (-k,2l), the homological degree is l-k-1. 

For each simplex input the vertices it contains separated by spaces
For example, if you want to set K=C_4 (the four cycle), you'd write it as:

1 2
1 4
2 3
3 4


You can either choose to prepare the simplicial complex in advancein the text file called "Input.txt" or enter it 
during runtime


Make sure the vertices are indexed by a set of integers (preferrably [n]:={1,2,3,...,n})
and that they are written in increasing order for each simplex.