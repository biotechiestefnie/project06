# Introduction
Description of the project

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


Smith_Waterman():
1. get lengths of both sequences
2. create scoring matrix filled with zeros
3. track maximum local alignment score
4. fill scoring matrix
    - determine match or mismatch score
    - calculate gap scores
    - smith-waterman recurrence
    - update max score
5. find the best possible score for the shorter sequence
6. avoid division by zero ( max possible score must be different from 0)
7. return (max score / max possible score)

build_distance_matrix():
1. get sequence IDs in a list to preserve order
2. get number of sequence
3. create empty square matrix filled with zeros
4. compare each pair of sequences
    - calculate similarity score using smith-waterman
    - convert similarity to distance
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
```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus
