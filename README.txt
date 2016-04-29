# Project 2: Pointerless Tree (Part I)
	Shubham Bansal	sb3766  
	Mark Rampton	mcr2176

## build instructions:

	make all


## how to run

	./main K P i1 i2 ... in

	K 			number of random keys to populate the tree with
	P 			number of probes to search the tree with
	i1 - in 	fanout of each level of the tree; root level is i1

## description
	
	The program accepts input from the user and builds and populates 
	a tree to match the specified structure. The program writes the
	results of the P number of probes against this tree. Each probe
	returns a range within the tree where the probed valued would appear. 

	Important details:

        Part1:
    		While the structure is a tree, it is a tree implemented with 
    		arrays -- one array per level. Pointers are not used for 
    		traversing the tree and we instead rely on the inherent structure
    		of one level's array to the previous. 

    		The tree is populated by filling in the leave nodes (lowest level array)
    		first. Parent nodes are populated one at a time as each lower level 
    		child is filled. If there are not enough keys to populate the entire tree,
    		MAXINT is inserted in place only if the node already has at least one
    		key value inserted.  If no key values exist in the node, each element is 
    		instead supplied with a 0 value.

        Part2: 
            Probing is carried out and the results are displayed inline with the results
            from probing in part1.  Timing is added to allow for analysis and comparison
            of the different probing techniques.  
            
            Building and executing is carried out the same way -- and doing so results
            in both parts being tested simultaneously. 
	  
