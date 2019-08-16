# Centralized load-balancing strategies

## Main idea
The main idea of this strategy is that all the processors send their computational load to processor P<sub>0</sub>. 
P<sub>0</sub> then calculates the prefix sum of te load on elements ordered according to their global IDs. Then Proccesor zero
partitions the prefix sum array and broadcast the array to all the processors.  

## Element loads evaluation
The load on each element is
<a href="https://www.codecogs.com/eqnedit.php?latex=O(N^{4})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?O(N^{4})" title="O(N^{4})" /></a>.
Where N is `the polynomial order +1`/ `number of the grid point in one direction`.  

## Steps
### On the root processor
1. Organize the load as an one dimensional array of spectral elements. 
2. Compute the `prefix sum` of the 1D array.
3. Calculate the `Optimal Bottleneck` for each processor. 
4. Perform `H1` and `H2` Heuristics for partition. 
5. Reallocate the elements to each processor using `MPI_SCATTERV`

## Pseudo computational load
The program start with 4 element and they are evenly allocated to 4 processors. The elements achieve `hp-adaptivity`. 
The upper-bound of the element h-refinement is 4, i.e., each element can split up to four times. The h-refinement is represented by
the variable `H_DEPTH` and it is initialize with `1`. `p` stands for the polynomial order. The polynomial order start by `4` end with `12`. 

The code performs a random `hp-refinement` scheme, meaning that each element adapts randomly and there is no real numerical computation implemented inside the code. The program starts with `h-refinement`, each element has 30% chance to split. Once an element's `H_DEPTH` reaches 2, it has 30% chance to invoke `p-enrichment`. The load of each element is based on the `p-order`.



