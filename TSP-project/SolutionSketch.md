## Solution sketch

	1. Pick a first path
		a. Algorithm to use
	2. Pick a heuristic or combination thereof
		a. 2 opt
		b. How long? Until cant anymore?
	3. Repeat and choose best resulting path
		a. How many times?



### Options

	1. Pick a first path
		a. Nearest neighbour (Naive)
		b. Greedy 
		c. Insertion heuristics
		d. Christofides

	2. Heuristics
		a. 2-opt
		b. 3-opt
			i. remove 3 edges instead of 2
		c. k-opt
		d. Lin-Kerninghan
		e. Tabu-Search
		f. Simulated Annealing
		g. Genetic Algorithms
		h. Tour Data Structure


### Greedy first path
The Greedy heuristic gradually constructs a tour by repeatedly selecting the shortest edge and adding it to the tour as long as it doesn’t create a cycle with less than N edges, or increases the degree of any node to more than 2