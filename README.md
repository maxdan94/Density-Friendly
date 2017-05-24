# Reference

Danisch, M., Chan, T. H. H., & Sozio, M. (2017, April).  
Large Scale Density-friendly Graph Decomposition via Convex Programming.  
In Proceedings of the 26th International Conference on World Wide Web (pp. 233-242). International World Wide Web Conferences Steering Committee.


# To compile:

"gcc simpleDF.c -o simpleDF -O3"  
"gcc approxDF.c -fopenmp -o approxDF -O3"  
"g++ exactDF.cpp -fopenmp -fpermissive -o exactDF -O3"


# To execute:

"./simpleDF edgelist.txt iter densest.txt decomp.txt".
- "edgelist.txt" should contain the graph: one edge on each line separated by a space.
- "iter" is the number of iterations, e.g. "100".
- "res.txt" contains the found subgraph: "upperbound density size node1 node2 node3...". "upperbound" is an upperbound on the maximum density
- "decomp.txt" contains on each line a node and it's density-score

If you also want the statistics (upperbound density and size of the found subgraph at each iteration) type:  
"./ds edgelist.txt iter densest.txt decomp.txt stat.txt".
- "stat.txt" will contain these statistics in the following format: "iteration density size upperbound".

It will be a bit slower.

"./approxDF nthreads iter net.txt rates.txt pavafit.txt cuts.txt"  
"./exactDF nthreads iter net.txt rates.txt pavafit.txt cuts.txt exact.txt"

- nthreads is the number of threads to use
- iter is the number of iterations over all edges to perform
- net.txt should contain the graph (one edge on each line: 2 unsigned separated by a space)
- rates.txt will contain the density value for each node
- pavafit.txt will contain the profile given by the PAVA fit, that is the week approximation of the density-friendly ("size density density-upperbound" on each line)
- cuts.txt will contain the profile given by correct cuts, that is the strong approximation of the density-friendly ("size density density-upperbound" on each line)
- exact.txt will contain the exact density-friendly ("size density" of each subgraph on each line, not necesarily in decreasing order of density

Some information will be printed in the terminal.

# Initial contributors

Maximilien Danisch  
Technical consultants: T-H. Hubert Chan and Mauro Sozio  
Avril 2017  
http://bit.ly/maxdan94  
maximilien.danisch@gmail.com
