# Label-induced orbits of the G2 tree fan

## Result

We rooted every exact metric tree at taxon `7`, tested all `720` permutations of taxa `1,...,6`, and fixed taxon `7`. A permutation was retained only when it

1. permuted the 51 generator-ray metric signatures; and
2. its induced ray permutation preserved all 251 maximal triangles.

The ray test retained **12** permutations and the full fan-incidence test retained **12**. The final group is exactly the expected hexagonal Weyl group

```text
rotation   (1 5 3 4 2 6)
reflection (2 3)(5 6)
taxon 7    fixed
```

Thus the label-induced group has order **12** and is the dihedral Weyl group W(G2).

This is the permutation part visible after valuation. Algebraic Weyl representatives may also carry coordinate signs, but those signs disappear from the tropical tree data.

The resulting action has **11 ray orbits** and **32 triangle orbits**, with size distributions

```text
rays:      [1, 2, 3, 3, 3, 3, 6, 6, 6, 6, 12]
triangles: [2, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
```

## Why this differs from the old atlas orbits

The abstract simplicial automorphism group has order **24**, with 9 ray orbits and 24 triangle orbits. The label-induced Weyl group has index **2** in that group. It is the subgroup realized by permuting the six nonzero-weight taxa (with algebraic signs suppressed tropically). In particular, abstract triangle orbits of size 24 cannot be Weyl orbits because `|W(G2)|=12`.

The finer label-orbit partition splits the abstract ray orbits as follows:

| old abstract ray orbit | size | new label orbits |
|---:|---:|:---|
| 1 | 12 | `1` |
| 2 | 12 | `2, 3` |
| 3 | 6 | `4` |
| 4 | 6 | `5, 6` |
| 5 | 6 | `7` |
| 6 | 3 | `8` |
| 7 | 2 | `9` |
| 8 | 3 | `10` |
| 9 | 1 | `11` |

It splits the abstract triangle orbits as follows:

| old abstract orbit | size | new label orbits |
|---:|---:|:---|
| 1 | 12 | `1` |
| 2 | 24 | `2, 5` |
| 3 | 24 | `3, 4` |
| 4 | 12 | `6` |
| 5 | 12 | `7, 16` |
| 6 | 12 | `8` |
| 7 | 24 | `9, 18` |
| 8 | 24 | `10, 13` |
| 9 | 12 | `11` |
| 10 | 12 | `12` |
| 11 | 12 | `14, 20` |
| 12 | 12 | `15, 22` |
| 13 | 3 | `17` |
| 14 | 3 | `19` |
| 15 | 12 | `21` |
| 16 | 3 | `23` |
| 17 | 3 | `24` |
| 18 | 6 | `25` |
| 19 | 6 | `26` |
| 20 | 3 | `27` |
| 21 | 6 | `28` |
| 22 | 6 | `29, 31` |
| 23 | 6 | `30` |
| 24 | 2 | `32` |

## Ray orbits

| orbit | size | indices |
|---:|---:|:---|
| 1 | 12 | `1, 4, 23, 24, 38, 39, 42, 43, 48, 49, 50, 51` |
| 2 | 6 | `2, 25, 26, 28, 31, 32` |
| 3 | 6 | `3, 5, 16, 27, 46, 47` |
| 4 | 6 | `6, 7, 11, 12, 40, 41` |
| 5 | 3 | `8, 33, 34` |
| 6 | 3 | `9, 10, 30` |
| 7 | 6 | `13, 17, 19, 20, 44, 45` |
| 8 | 3 | `14, 21, 22` |
| 9 | 2 | `15, 29` |
| 10 | 3 | `18, 35, 36` |
| 11 | 1 | `37` |

## Triangle orbits

| orbit | size | indices |
|---:|---:|:---|
| 1 | 12 | `1, 2, 49, 50, 154, 155, 182, 183, 192, 193, 248, 249` |
| 2 | 12 | `3, 4, 174, 175, 176, 177, 212, 213, 214, 215, 234, 235` |
| 3 | 12 | `5, 6, 90, 91, 114, 115, 196, 197, 198, 201, 239, 242` |
| 4 | 12 | `7, 8, 47, 48, 116, 117, 178, 179, 202, 203, 246, 247` |
| 5 | 12 | `9, 10, 37, 38, 84, 85, 216, 217, 226, 227, 250, 251` |
| 6 | 12 | `11, 12, 53, 54, 65, 66, 146, 148, 222, 223, 224, 225` |
| 7 | 6 | `13, 56, 129, 130, 150, 151` |
| 8 | 12 | `14, 15, 39, 40, 67, 68, 75, 76, 166, 167, 206, 207` |
| 9 | 12 | `16, 17, 105, 106, 108, 109, 156, 157, 164, 169, 204, 209` |
| 10 | 12 | `18, 19, 170, 172, 188, 189, 190, 191, 210, 211, 219, 220` |
| 11 | 12 | `20, 21, 110, 111, 112, 113, 171, 173, 194, 195, 230, 231` |
| 12 | 12 | `22, 23, 142, 144, 160, 161, 199, 200, 236, 237, 240, 241` |
| 13 | 12 | `24, 25, 79, 80, 127, 128, 143, 145, 232, 233, 244, 245` |
| 14 | 6 | `26, 63, 77, 78, 125, 126` |
| 15 | 6 | `27, 61, 100, 102, 131, 132` |
| 16 | 6 | `28, 45, 46, 55, 218, 221` |
| 17 | 3 | `29, 147, 149` |
| 18 | 12 | `30, 31, 33, 34, 51, 52, 88, 89, 228, 229, 238, 243` |
| 19 | 3 | `32, 152, 153` |
| 20 | 6 | `35, 36, 57, 59, 180, 181` |
| 21 | 12 | `41, 42, 73, 74, 104, 107, 118, 119, 165, 168, 205, 208` |
| 22 | 6 | `43, 44, 58, 62, 184, 185` |
| 23 | 3 | `60, 69, 70` |
| 24 | 3 | `64, 92, 96` |
| 25 | 6 | `71, 72, 101, 103, 162, 163` |
| 26 | 6 | `81, 82, 93, 95, 98, 99` |
| 27 | 3 | `83, 94, 97` |
| 28 | 6 | `86, 87, 120, 124, 186, 187` |
| 29 | 3 | `121, 136, 141` |
| 30 | 6 | `122, 123, 137, 138, 139, 140` |
| 31 | 3 | `133, 134, 135` |
| 32 | 2 | `158, 159` |

## Files

- `label_permutation_group.csv`: every label permutation and its induced actions;
- `ray_label_orbits.csv`: exact rooted-metric ray orbits;
- `triangle_label_orbits.csv`: triangle orbits, ray triples, fiber types, and rooted topologies.

The paper ray catalog and the flag atlas read these files. Reproduce them from the repository root with

```sh
julia --project=julia julia/studies/g2/atlas/grassmannian/compute_label_orbits.jl
```

The maintained coarsest-fan symmetry artifacts are recomputed instead with

```sh
julia --project=julia julia/studies/g2/run_census.jl --target=base --stage=symmetry
```
