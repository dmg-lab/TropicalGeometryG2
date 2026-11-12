# Observations

## Variety

Automorphism group of $K$: G

On 51 vertices:

orbits(GV)
9-element Vector{Oscar.GSetByElements{PermGroup, Int64}}:
 - G-set of GV with seeds [1, 24, 39, 4, 50, 23, 38, 51, 49, 43, 48, 42] -> path of length 2
 - G-set of GV with seeds [2, 5, 16, 32, 3, 31, 26, 47, 25, 46, 27, 28] -> path of length 2
 - G-set of GV with seeds [6, 7, 12, 11, 41, 40] -> ordinary split
 - G-set of GV with seeds [8, 10, 30, 33, 9, 34] -> ordinary split, 30 (36|12457), 9 (14|23567), 10 (25|13467), 8 (367|1245), all in the link of 37 (snowflake), correspond to the terms $x_1 x4+x_2x_5+x_3x_6$ in the quadric and to the terms $e^{147} + e^{247}+e^{467}$ in the trilinear alternating form.
 - G-set of GV with seeds [13, 20, 19, 45, 44, 17] -> path of length 2
 - G-set of GV with seeds [14, 22, 21] -> path of length 2
 - G-set of GV with seeds [15, 29] -> 123|456: 1237|456 or 123|4567 correspond to the two terms $e^{123}+e^{456} of the trilinear alternating form. 
 - G-set of GV with seeds [18, 35, 36] -> path of length 2
 - G-set of GV with seeds [37] -> only fixed point : $e_{14}+e_{25}+e_{36}$. special snowflake (14|25|36)|7



-----
On 251 facets:
orbits(GF)
24-element Vector{Oscar.GSetByElements{PermGroup, Int64}}:
 - G-set of GF with 12 seeds [1, 50, 155, 2, 249, 49, 154, 248, 183, 193, 182, 192]
 - G-set of GF with 24 seeds [3, 38, 85, 4, 234, 37, 84, 10, 235, 213, 214, 9, 217, 212, 215, 175, 216, 226, 251, 174, 176, 227, 250, 177]
 - G-set of GF with 24 seeds [5, 48, 117, 6, 239, 47, 116, 8, 242, 198, 197, 7, 203, 196, 201, 115, 202, 178, 247, 114, 90, 179, 246, 91]
 - G-set of GF with 12 seeds [11, 54, 65, 12, 224, 53, 66, 222, 148, 225, 146, 223]
 - G-set of GF with 12 seeds [13, 46, 28, 129, 45, 130, 151, 221, 150, 218, 55, 56]
 - G-set of GF with 12 seeds [14, 40, 76, 15, 207, 39, 75, 206, 68, 167, 67, 166]
 - G-set of GF with 12 seeds [16, 52, 30, 17, 157, 51, 31, 34, 156, 209, 105, 33, 243, 204, 108, 169, 238, 88, 229, 164, 106, 89, 228, 109]
 - G-set of GF with 24 seeds [18, 80, 145, 19, 219, 79, 143, 25, 220, 211, 188, 24, 128, 210, 191, 172, 127, 232, 244, 170, 189, 233, 245, 190]
 - G-set of GF with 12 seeds [20, 112, 173, 21, 230, 110, 171, 231, 195, 113, 194, 111]
 - G-set of GF with 12 seeds [22, 144, 161, 23, 241, 142, 160, 240, 237, 200, 236, 199]
 - G-set of GF with 12 seeds [26, 36, 57, 125, 35, 126, 78, 181, 77, 180, 59, 63]
 - G-set of GF with 12 seeds [27, 44, 58, 131, 43, 132, 102, 185, 100, 184, 62, 61]
 - G-set of GF with 3 seeds [29, 147, 149]
 - G-set of GF with 3 seeds [32, 152, 153]
 - G-set of GF with 12 seeds [41, 119, 42, 118, 73, 208, 74, 205, 168, 104, 165, 107]
 - G-set of GF with 3 seeds [60, 69, 70]
 - G-set of GF with 3 seeds [64, 92, 96]
 - G-set of GF with 6 seeds [71, 103, 72, 101, 163, 162] -> non-trivalent
 - G-set of GF with 6 seeds [81, 98, 99, 82, 93, 95]
 - G-set of GF with 3 seeds [83, 94, 97]
 - G-set of GF with 6 seeds [86, 124, 87, 120, 187, 186]
 - G-set of GF with 6 seeds [121, 135, 134, 136, 133, 141] -> non-trivalent
 - G-set of GF with 6 seeds [122, 137, 123, 139, 138, 140]
 - G-set of GF with 2 seeds [158, 159] -> non-trivalent

 ## Prevariety
 -> indices refer to the rays modulo lineality indices for the fan of the variety (there is a permutation going on for the rays_modulo_lineality of the Prevariety) 
 ### 3 tetrahedra

 - [8,11,12,30] 
 - [7,10,33,41]
 - [6,9,34,40]
 ### 3 squares (the same three squares that are cut in half by the edges with triangle incidence 2 in the variety)
 - [13,14,17,18] 
 - [20,22,35,44]
 - [19,21,36,45]


## Toric degenerations




 ## Chevalley Polytope: the lecture hall simplex $(1,3,2,3,1)$


 Exhaustive UniTriSat enumeration gives 418 unimodular triangulations. All are regular and 180 are flag. The earlier “single triangulation” result omitted `find_all=true`.

-> switching to the face $x_1=1$ because it's a cone over it: LH(3,2,3,1)

Its projected $x_1=1$ facet also has 418 unimodular triangulations, all regular and 180 flag.

## Flag variety

$f$-vector: $(1, 344, 2515, 5709, 3879)$

Coarsest structure f-vector: $(141, 1352, 3584, 2714)$
Each maximal cell of the flag variety correspond to a maximal G2-tree decorated by a point in the relative interior of an edge.
Each ray of the flag variety corresponds to a minimal G2-tree whose one of its vertices is decorated. 
Still need to analyze the other cells.

