## Create unit disk graphs

Use UnitDiskGraph(N, rho,1) to create a unit disk graph with N vertices, density rho, and periodic boundary condition. 

## Benchmarking instances

Data Files for benchmarking in the paper arXiv:1808.10816
1) Fig. 3a - UDGp-instances-for-phasediag.mat
2) Fig. 3b - UDGp2-instances-for-phasediag-line.mat

The variable xyCoords contains the x-y coordinates of vertices.
xyCoords{i,j,k} gives the coordinates for a graph with
	size N = NvertToTry(i)
	density rho = rhos(j)
	instance # k

Actual graph edge list and adjacency matrix can be generated using UnitDiskGraph.m
The benchmarking are done assuming the vertex are placed on a torus of size sqrt(N/rho) x sqrt(N/rho), i.e. square with periodic boundary condition.

For example, use

	[edges, A] = UnitDiskGraph(xyCoords{i,j,k}, rhos(j), 1)

to generate edge list 'edges' and adjacency matrix 'A'