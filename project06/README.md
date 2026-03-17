# Introduction
This project is an application of the Neighbor-Joining algorithm for analyzing distance-based phylogenetic relationships among species.  NJ works by identifying a pair of nodes that minimize a specific criterion (Qi,j) and joining them into a new intermediate node, reducing the taxa count by one for each iteration. The Q criterion is represented by the equation:
                            Qi,j = (n-2)Di,j - Si -Sj
where:
- Di,j = the distance between nodes i and j
- n = current number of nodes
- Si = sum of distances from all other nodes to i
- Sj = sum of distances from all other nodes to j

Neighbor-Joining (NJ) differs from UPGMA in that it produces an unrooted tree and does not assume a constant evolutionary rate across lineages. This makes NJ more biologically realistic for determining sequence-based relationships. However, the canonical NJ method has a complexity of O(n^3), making it highly sensitive to large datasets. One solution for this is optimizing the algorithm with dynamic programming, known as Dynamic Neighbor Joining (DNJ). With dynamic programming, NJ can be scaled to handle over 100,000 taxa, reducing the complexity for finding minimum taxa to O(dn) and resulting in a total time complexity of O(dn^2) (Claussen, 2023). This approach is mathematically guaranteed to result in the same unrooted tree, but in substantially less time by optimizing the selection of the best pair to join. Instead of recalculating the entire Q matrix in every step, DNJ maintains the minimum join criterion for each row in a vector (Bryant & Moultin, 2004). The algorithm updates the Q matrix in linear time when a new node is formed rather than full re-computation, providing substantial speedups for large datasets.

For even larger datasets, a heuristic version (HNJ) of the algorithm exists. HNJ updates the Q matrix differently, approximating the search by only checking the new node's distances. While this reduces the complexity to O(n^2), it may not always produce the exact same tree as NJ and DNJ, serving as a prime example of the sacrifices that must be made to accommodate larger datasets.

Local alignment with Smith-Waterman (SW) can be used to calculate distances for NJ because SW yields a pairwise similarity score that reflects how much two biological sequences resemble each other in their most conserved region. The key is that SW produces similarity, while NJ requires distance, so the similarity scores must be transformed in a way that preserves the mathematical assumptions of NJ. This is accomplished with the following equation:
                    d(i,j) = max(SW) - SW(i,j)
where: 
- d(i,j) = the distance between nodes i and j
- max(SW) = the maximum Smith-Waterman score
- SW(i,j) = the raw Smith-Waterman local alignment score between sequences i and j

Smith–Waterman scores are raw integers determined by sequence length, match/mismatch scoring, and gap penalties. They are not normalized, and they have no fixed upper bound. A pair of long, highly similar sequences might produce a score of 250, while a shorter pair might only reach 120 even if they are equally related biologically. Because the scale of SW scores varies across pairs, the only way to convert them into a distance that NJ can interpret is to use a global linear inversion. By subtracting each similarity score from the single largest score in the entire matrix, you create a distance matrix where more similar sequences have smaller distances, less similar sequences have larger distances, and all distances lie on a shared global scale. This transformation is monotonic, preserves ordering, and guarantees non‑negative distances, which are all requirements for NJ to behave correctly.

The alternative formula 𝑑 = 1 − 𝑠 is only valid when similarity values are already normalized to the interval [0,1] and share a fixed maximum of 1. Measures like cosine similarity or Jaccard similarity satisfy this requirement, but raw SW scores do not. Applying 1 − 𝑠 directly to SW scores produces negative distances, which violates NJ’s assumptions and breaks the algorithm. If SW scores were normalized first, they must be normalized globally using: sw' = sw(i,j) / sw(max) because normalizing each pair independently destroys comparability across the matrix. With global normalization, the transformation d = 1 - s' becomes mathematically equivalent to d(i,j) = sw(max) - sw(i,j), just scaled by a constant factor. NJ is invariant under positive scaling, so both produce the same topology. This is why the max‑similarity inversion is the standard and simplest method for converting raw alignment scores into NJ distances.

As such, we elected to apply the equation d(i,j) = SW(max) - SW(i,j) rather than distance = 1 - similarity in our program to ensure we obtained non-distorted distances that met all the requirements of the NJ algorithm.

 

# Pseudocode
Put pseudocode in this box:

```
Read fasta file:
1. initialize storage
2. Loop through file
3. Detect header ">"
4. Accumulate sequence lines
5. Save sequences
6. Return dictionary
Note: We considered returning an ordered dictionary because mapping is created the moment the FASTA file is read and if the order of sequences changes, the following downstream components change: the distance matrix layout, the Q‑matrix calculations, the pair chosen for joining at each iteration, the shape and branch lengths of the final tree, and the Newick output order. However, due to the small size of the dataset and that Python 3.7+ already preserves order, we elected to avoid the additional overhead of an ordered dictionary by using a regular dictionary. 

Smith_Waterman():
1. get lengths of both sequences
2. create scoring matrix filled with zeros
3. track maximum local alignment score
4. fill scoring matrix
    - determine match or mismatch score
    - calculate gap scores
    - smith-waterman recurrence
    - update max score
5. Return max score

build_distance_matrix():
1. get sequence IDs in a list to preserve order
2. get number of sequence
3. create empty square matrix filled with zeros
4. compare each pair of sequences
    - calculate similarity score using smith-waterman
    - convert similarity to distance using d(i,j) = SW(max) - SW(i,j)
    - fill both symmetric positions in matrix
5. return distance matrix, sequence ids

neighbor_joining()
1. make copies to original inputs are not modified
2. continue joining until only two nodes remain
    -compute total distance for each taxon
    -build Q matrix
    -find pair with minimum Q value
    -make sure i<j for easier removal later
    -calculate branch lengths from i to j to new node
    -prevent tiny negative values from floating point issues
    -create new joined label in Newick format
    -compute distances from new node to all remaining nodes
    -build reduced distance matrix
    -copy old distances among kept nodes
    -add distances from the new node
    -update labels
3. final join when only two nodes remain

plot tree()
- plots unrooted phylogenetic tree from newick string using matplotlib and biopython
1. Parse the Newick text into a tree structure
    tree_object ← PARSE_NEWICK(newick_string)
2. Create a drawing canvas with specified dimensions
    canvas ← CREATE_CANVAS(width = 8, height = 10)
3. Draw the tree on the canvas
    # Only show labels for terminal (leaf) nodes
    DRAW_TREE(
        tree = tree_object,
        canvas = canvas,
        show_internal_labels = FALSE,
        label_function = IF node IS LEAF THEN RETURN node.name ELSE RETURN NOTHING
    )

4. Reduce font size of all text labels
    FOR each label IN canvas.labels:
        SET_FONT_SIZE(label, size = 6)

5. Add margins around the drawing to reduce overlap
    SET_MARGINS(canvas, x_margin = 0.1, y_margin = 0.05)

6. Optimize layout to fit the figure area
    ADJUST_LAYOUT(canvas)

7. Display the final rendered tree
    SHOW(canvas)



```

# Successes
Combining last week's dynamic programming with neighbor joining to create a computationally efficient phylogenetic tree clarified a lot of misunderstanding for us regarding graph algorithm application.

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus


# References
Bryant, D. & Moulton, V. (2004). Neighbor-Net: An agglomerative method for the construction of phylogenetic networks. *Molecular Biology and Evolution*, 21(2): 255–265. https://doi.org/10.1093/molbev/msh018
Clausen, PTLC. (2023). Scaling neighbor joining to one million taxa with dynamic and heuristic neighbor joining. *Bioinformatics*,39(1): btac774. https://doi.org/10.1093/bioinformatics/btac774
