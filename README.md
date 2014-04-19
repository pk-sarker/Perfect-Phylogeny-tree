Implementation of Perfect Phylogeny Tree with C
-----------------
|    PROBLEM    |
-----------------

-------
Implement the algorithm for testing if a character-state matrix, with 2 states for each character,

admits a prefect phylogeny. When it does, display the phylogeny tree. Even when it doesn't, you can

display the phylogeny tree for the first i<=i<=m compatible characters. 

********************

-------------------------
|    Global Variables   |
-------------------------
 * N : Number of Taxon
 * M : Number of Characters in each Taxon
 * isPhylo : Flag which indicates whether it is a perfect phylogeny or not
 * MAT[][] : Character-State Matrix
 * MAT1[][] : Character-State Matrix where the columns are ordered in descending base on number of 1's per column 
 * O[][] : O[i] set
 * vs[][] : Vertex set those are already divided by characters or O[i].
 * vc : Number of vertexes in vs[][]
 * MO[][] : O[i] set created for newly formed character-state matrix M'
 * s1[] : Set of taxon that has 1's in a particular column
 * s0[] : Set of taxon that has 0's in a particular column
 *  



 --------------------
 |   Test Input     |
 -------------------- 
 Number of Taxon(n) : 5
 
 Number of Characters in each Taxon (m): 5
 
  
  Taxon 1 = 1 1 0 0 0
  
  Taxon 2 = 0 0 1 0 0
  
  Taxon 3 = 1 1 0 0 1
  
  Taxon 4 = 0 0 1 1 0
  
  Taxon 5 = 0 1 0 0 0
  

