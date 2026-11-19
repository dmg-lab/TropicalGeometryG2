#!/usr/bin/env bash
###############################################################################
# Build the Zenodo tarball of the paper, TropicalGeometryG2_data.tar.gz, at the
# repository root. The script only copies stored files; it computes nothing.
# A missing file aborts it. make_paper_data.jl (in julia/) writes every file.
#
# The tarball holds the files of the table tab:data of the paper at its top
# level, and in data/ a copy of data/ in julia/, as make_paper_data.jl writes
# it. Unpacking data/ into julia/ is the same as running make_paper_data.jl;
# from julia/, with the tarball at the repository root:
#
#   tar -xzf ../TropicalGeometryG2_data.tar.gz --strip-components=1 TropicalGeometryG2_data/data
#
# Run from anywhere; from julia/:
#
#   bash ../production/make_zenodo_tarball.sh [TARBALL]
#
# (TARBALL defaults to TropicalGeometryG2_data.tar.gz at the repository root,
# ../ from julia/).
###############################################################################

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$ROOT/julia/data"
G2="$DATA/g2"
NAME="TropicalGeometryG2_data"
TARBALL="$(realpath -m "${1:-$ROOT/$NAME.tar.gz}")"

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
OUT="$STAGE/$NAME"

# toric_orbits CSV: the orbits, as three digits, of the lines of a stage-6 CSV whose
# column `toric` is true. A field in double quotes may contain commas.
toric_orbits() {
    julia --startup-file=no - "$1" <<'JULIA'
function fields(line)
    out, field, quoted, i = String[], IOBuffer(), false, 1
    while i <= lastindex(line)
        c = line[i]
        if quoted && c == '"' && i < lastindex(line) && line[nextind(line, i)] == '"'
            print(field, '"'); i = nextind(line, i)
        elseif c == '"'
            quoted = !quoted
        elseif c == ',' && !quoted
            push!(out, String(take!(field)))
        else
            print(field, c)
        end
        i = nextind(line, i)
    end
    push!(out, String(take!(field)))
    return out
end
lines = readlines(ARGS[1])
header = fields(lines[1])
orbit, toric = findfirst(==("orbit"), header), findfirst(==("toric"), header)
for line in lines[2:end]
    row = fields(line)
    length(row) == length(header) || error("$(ARGS[1]): a line has $(length(row)) fields, not $(length(header))")
    row[toric] == "true" && println(lpad(row[orbit], 3, '0'))
end
JULIA
}

# put SOURCE DESTINATION: copy one file into the tarball, relative to $OUT.
put() {
    [[ -f "$1" ]] || { echo "missing: $1" >&2; exit 1; }
    mkdir -p "$(dirname "$OUT/$2")"
    cp "$1" "$OUT/$2"
}

# 00: the ideals of the point space Q, the G2-Grassmannian G and the flag variety F.
put "$G2/pointspace/00_ideal.mrdi"       00_ideals/pointspace_qq.mrdi
put "$G2/pointspace_gf2/00_ideal.mrdi"   00_ideals/pointspace_gf2.mrdi
put "$G2/grassmannian/00_ideal.mrdi"     00_ideals/linespace_qq.mrdi
put "$G2/grassmannian_gf2/00_ideal.mrdi" 00_ideals/linespace_gf2.mrdi
put "$G2/flag/00_ideal.mrdi"             00_ideals/flagspace_qq.mrdi
put "$G2/flag_gf2/00_ideal.mrdi"         00_ideals/flagspace_gf2.mrdi

# The orbit correspondences between the orbit types of main.tex and the census
# orbits go to the root. The fans are the normalized coarsest fans of the
# census, whose numbering of rays and maximal cells the correspondences use.
put "$G2/dressian/dressian_orbits.csv"           dressian_orbits.csv
put "$G2/flag_dressian/flag_dressian_orbits.csv" flag_dressian_orbits.csv
put "$G2/grassmannian_orbits.csv"                grassmannian_orbits.csv
put "$G2/flag_orbits.csv"                        flag_orbits.csv

# 01: the G2-Dressian and the G2 flag Dressian.
put "$G2/dressian/03_normalized_coarsest_fan.mrdi"      01_dressians/dressian.mrdi
put "$G2/flag_dressian/03_normalized_coarsest_fan.mrdi" 01_dressians/flag_dressian.mrdi

# 02: trop(G) and trop(F) over QQ and over GF(2).
put "$G2/grassmannian/03_normalized_coarsest_fan.mrdi"     02_tropical_varieties/grassmannian_qq.mrdi
put "$G2/grassmannian_gf2/03_normalized_coarsest_fan.mrdi" 02_tropical_varieties/grassmannian_gf2.mrdi
put "$G2/flag/03_normalized_coarsest_fan.mrdi"             02_tropical_varieties/flag_qq.mrdi
put "$G2/flag_gf2/03_normalized_coarsest_fan.mrdi"         02_tropical_varieties/flag_gf2.mrdi

# 03: the Newton-Okounkov polytopes of the toric orbits of each census. The
# .poly files of stage 6 are Oscar mrdi files.
for target in grassmannian grassmannian_gf2 flag flag_gf2; do
    case "$target" in
        *_gf2) folder="$target" ;;
        *)     folder="${target}_qq" ;;
    esac
    orbits="$(toric_orbits "$G2/$target/06_canonical_interior_degenerations.csv")"
    [[ -n "$orbits" ]] || { echo "no toric orbits in $target" >&2; exit 1; }
    for orbit in $orbits; do
        put "$G2/$target/06_canonical_interior_degenerations_objects/orbit_${orbit}_newton_okounkov_body.poly" \
            "03_no_polytopes/$folder/orbit_${orbit}.mrdi"
    done
done

# 04: the tropicalized hexagon of fig:tropical_hexagon.
put "$G2/hexagon/tropicalized_hexagon.mrdi" 04_hexagons/tropicalized_hexagon.mrdi

# 05: the two genus-10 curves of fig:planar_curve and fig:nonplanar_curve, and
# the pieces of the curve of Method 1, computed from its shift q, from which
# the stage zenodo builds it.
put "$DATA/mukai_curves/genus_10_method_1.mrdi"        05_curves/genus_10_method_1.mrdi
put "$DATA/mukai_curves/genus_10_method_2.mrdi"        05_curves/genus_10_method_2.mrdi
put "$DATA/mukai_curves/genus_10_method_1_pieces.mrdi" 05_curves/genus_10_method_1_pieces.mrdi
ALIASES="$(find "$OUT" -type f | wc -l)"

# data/: the folder data/ of julia/ as make_paper_data.jl writes it, with the
# inputs that it reads and that git keeps. The list names what to copy, so that
# nothing else in that folder goes into the tarball.
DATA_FILES=(
    README.md
    fan_uniform_linear_space_16_21.mrdi
    mukai_curves/genus_10_method_1.mrdi
    mukai_curves/genus_10_method_1_pieces.mrdi
    mukai_curves/genus_10_method_1_k3.mrdi
    mukai_curves/genus_10_method_2.mrdi
    g2/SUMMARY.md
    g2/grassmannian_orbits.csv
    g2/flag_orbits.csv
    g2/maximal_cell_tree_metrics.csv
    g2/provenance/zenodo.toml
    g2/hexagon/tropicalized_hexagon.mrdi
)
DATA_FOLDERS=(dressian flag_dressian grassmannian grassmannian_gf2 flag flag_gf2 grassmannian_Slong
              pointspace pointspace_gf2)
# Files of these folders that no step of make_paper_data.jl writes: the
# report of the checks of the census over GF(2) (verify_gf2_census.jl).
DATA_EXCLUDED='^g2/grassmannian_gf2/07_verification\.toml$'
for file in "${DATA_FILES[@]}"; do
    put "$DATA/$file" "data/$file"
done
for folder in "${DATA_FOLDERS[@]}"; do
    [[ -d "$G2/$folder" ]] || { echo "missing: $G2/$folder" >&2; exit 1; }
    while IFS= read -r -d '' file; do
        relative="${file#"$DATA/"}"
        [[ "$relative" =~ $DATA_EXCLUDED ]] || put "$file" "data/$relative"
    done < <(find "$G2/$folder" -type f -print0 | sort -z)
done
# The copies in data/ without their private blocks (the lines from <!-- private -->
# to <!-- /private --> in Markdown, from # private to # /private in TOML), by
# strip_private.awk.
while IFS= read -r -d '' file; do
    status=0
    awk -f "$ROOT/production/strip_private.awk" "$file" > "$file.public" || status=$?
    case "$status" in
        0) mv "$file.public" "$file" ;;
        1) rm -f "$file.public" ;;
        *) echo "unbalanced private block in ${file#"$OUT/"}" >&2; exit 1 ;;
    esac
done < <(find "$OUT/data" -type f -print0)

cat > "$OUT/README.md" <<'EOF'
# Tropical Geometry of G2: data

The data of the article

> Michael Joswig, Kevin Kuehn, Bernd Sturmfels and Mathieu Vallée, *Tropical Geometry of G2*,
> arXiv:2611.24456 (2026), https://arxiv.org/abs/2611.24456

- Dataset: https://doi.org/10.5281/zenodo.24104003
- Code: https://github.com/dmg-lab/TropicalGeometryG2

This file summarizes Appendix A of the article. Equation, figure, table, theorem and section numbers are those of
arXiv:2611.24456v1. Q, G and F are the quadric G2/P1 in P^6, the G2-Grassmannian G2/P2 in P^20 and the flag
variety G2/B in P^6 x P^20.

## Contents

The 222 files at the top level are the objects of the article (Table 4). The `.mrdi` files are in the
serialization format of OSCAR and are read by its function `load`. They were written by OSCAR 1.8.2, which each file
names in its header.

| file                                                   | content                                               | article             |
|--------------------------------------------------------|-------------------------------------------------------|---------------------|
| `00_ideals/pointspace_{qq,gf2}.mrdi`                   | the ideal of Q: one quadric                           | (2)                 |
| `00_ideals/linespace_{qq,gf2}.mrdi`                    | the ideal of G: 7 linear forms, 35 quadrics           | (4), (5)            |
| `00_ideals/flagspace_{qq,gf2}.mrdi`                    | the ideal of F: 78 generators                         | Section 2           |
| `01_dressians/dressian.mrdi`                           | the G2-Dressian: 51 rays, 239 maximal cones           | Theorem 2.4         |
| `01_dressians/flag_dressian.mrdi`                      | the G2 flag Dressian: 141 rays, 2627 maximal cones    | Theorem 2.5         |
| `02_tropical_varieties/grassmannian_{qq,gf2}.mrdi`     | trop(G): 51, 54 rays; 248, 254 maximal cones          | Theorems 3.1, 4.1   |
| `02_tropical_varieties/flag_{qq,gf2}.mrdi`             | trop(F): 141, 156 rays; 2714, 2792 maximal cones      | Theorems 3.3, 4.2   |
| `03_no_polytopes/grassmannian_{qq,gf2}/orbit_NNN.mrdi` | 17, 20 Newton-Okounkov polytopes of G                 | Theorem 6.3, Table 2 |
| `03_no_polytopes/flag_{qq,gf2}/orbit_NNN.mrdi`         | 81, 83 Newton-Okounkov polytopes of F                 | Theorem 6.3         |
| `04_hexagons/tropicalized_hexagon.mrdi`                | the tropicalized hexagon of Figure 15                 | Example 9.3         |
| `05_curves/genus_10_method_1.mrdi`                     | the curve of genus 10 of Figure 16                    | Proposition 10.1    |
| `05_curves/genus_10_method_2.mrdi`                     | the curve of genus 10 of Figure 17                    | Proposition 10.1    |
| `05_curves/genus_10_method_1_pieces.mrdi`              | the polyhedra computed for Method 1                   | Proposition 10.1    |
| `grassmannian_orbits.csv`                              | the orbit types of trop(G)                            | Theorems 3.1, 4.1   |
| `flag_orbits.csv`                                      | the orbit types of trop(F)                            | Theorems 3.3, 4.2   |
| `dressian_orbits.csv`                                  | the orbit types of the G2-Dressian                    | Theorem 2.4         |
| `flag_dressian_orbits.csv`                             | the orbit types of the G2 flag Dressian               | Theorem 2.5         |
| `data/`                                                | the intermediate results, as the code writes them     | Appendix B          |
| `README.md`                                            | this file                                             |                     |

The folder `data/` is the folder `julia/data/` of the code repository after a full computation: every stage of
every computation, from which the files above are copied and from which the code checks the statements of the
article. Its file `README.md` describes it.

## Conventions

- The suffixes `qq` and `gf2` name the field QQ or GF(2). By Theorems 4.1, 4.2 and 6.3, the fans and polytopes
  over QQ hold for every characteristic other than 2.
- The coordinates of P^6 are x1, ..., x7, those of P^20 are p12, ..., p67 in lexicographic order, and those of
  P^6 x P^20 are x1, ..., x7, p12, ..., p67. The fans lie in R^21 and R^28.
- All tropical objects in the files use the min-convention of OSCAR. A point w in a file is -w in the
  max-convention of the article: the G2-tree of a ray r is -r, up to lineality, and the coordinates x_k give the
  d_k8 of the flag Dressian.
- The fans are `PolyhedralFan`s with the coarsest fan structure of Sections 2-4. Their rays are numbered as in
  `rays_modulo_lineality`, and their maximal cones as in `maximal_cones`.

## Orbit types

The computation numbers the rays and maximal cones of each fan, and their orbits under the Weyl group D12. The
four tables `*_orbits.csv` translate these orbits into the types of the article, one orbit per line. The members
of the orbits in `grassmannian_orbits.csv` and `flag_orbits.csv` are numbered as in the files of the same names in
`02_tropical_varieties/`, those of the other two tables as in `01_dressians/`. An orbit over QQ and an orbit over
GF(2) share a line when they consist of the same rays or cones. The columns:

- `kind`: `ray` or `cell`.
- `orbit_type`: the type `alpha_1`, ..., `epsilon` of a ray, as in Figure 7, or the type `1`, ..., `35` of a
  maximal cone, as in Figures 8 and 9. A ray or maximal cone of a flag fan has the type of the ray or maximal cone
  of trop(G) below it, and a ray over the lineality space has the type `star`. The tetrahedra of the G2-Dressian,
  and the cones of the flag Dressian over them, have the type `epsilon` (Figure 6).
- `apple` (flag fans): the edges of the G2-tree of the cell type that carry the apple, as in Figure 11, separated
  by `|`. `pendant k` is the pendant edge of the leaf k, and `inner a-b` the inner edge between the inner nodes
  a, b in {8, ..., 12}.
- `size`: the size of the orbit.
- `orbit` (Dressians), `qq_orbit`, `gf2_orbit`: for a cell, the number of its orbit, which is the number NNN in
  `03_no_polytopes/` for a toric cell type; for a ray, the smallest index of a ray in the orbit. Empty if the orbit
  does not occur over that field.
- `members`, `qq_members`, `gf2_members`: the indices of the rays or maximal cones of the orbit.

For example, the line

    kind,orbit_type,size,qq_orbit,gf2_orbit,qq_members,gf2_members
    ray,alpha_2,3,8,1,8;33;34,1;19;22

of `grassmannian_orbits.csv` says that the rays of type alpha_2 form an orbit of size 3: the rays 8, 33 and 34 of
`grassmannian_qq.mrdi`, and in characteristic 2 the rays 1, 19 and 22 of `grassmannian_gf2.mrdi`.

## Newton-Okounkov polytopes

`03_no_polytopes/grassmannian_qq/orbit_NNN.mrdi` is the polytope P_C of Section 6 for a maximal cone C in the
orbit NNN of `02_tropical_varieties/grassmannian_qq.mrdi`, and likewise for the other three folders. There is one
file for each orbit of toric degenerations in Theorem 6.3. The polytopes are lattice polytopes of dimension 5 in
R^5 for G, and of dimension 6 in R^6 for F, where P_C is the Minkowski sum of the polytopes of the coordinates
x1, ..., x7 and of the coordinates p12, ..., p67.

## Hexagon and curves

`tropicalized_hexagon.mrdi` is a named tuple with the fields `points`, the six vectors val(x^(1)), ..., val(x^(6))
in QQ^7, and `lines`, the six vectors val(p^(12)), val(p^(23)), ..., val(p^(61)) in QQ^21, each shifted to the
minimum 0. For such a vector v, the vector a^(i) or d^(ij) of Example 9.3 is max(v) 1 - v.

`genus_10_method_1.mrdi` and `genus_10_method_2.mrdi` hold the tropical curves of genus 10 obtained with
Methods 1 and 2 of Section 10, as one-dimensional `PolyhedralComplex`es in R^21 with one point over each point of
R^21/1. They have 449 and 447 maximal cells, of which 239 and 210 are unbounded, and their bounded cells form
graphs of genus 10 on 201 and 228 vertices. The first file is a named tuple with the field `complex`; the second
also has the field `multiplicities`, the multiplicities 1, 2, 4 or 5 of the maximal cells in the order of
`maximal_polyhedra`. `genus_10_method_1_pieces.mrdi` is the vector of the 454 one-dimensional polyhedra computed
for Method 1 from its shift q; they make up the 449 maximal cells of `genus_10_method_1`, five of which consist of
two collinear pieces.

## Example

In Julia with OSCAR, from this folder:

```julia
using Oscar
G = load("02_tropical_varieties/grassmannian_qq.mrdi")   # trop(G), a PolyhedralFan in R^21
lineality_dim(G), n_maximal_cones(G)                     # (3, 248)
r = rays_modulo_lineality(G).rays_modulo_lineality[8]    # a ray of type alpha_2; its G2-tree is -r
P = load("03_no_polytopes/grassmannian_qq/orbit_001.mrdi")
dim(P), normalized_volume(P)                             # (5, 18)
H = load("04_hexagons/tropicalized_hexagon.mrdi")
maximum(H.points[1]) .- H.points[1]                      # a^(1) = (11, 12, 4, 15, 0, 14, 13)
C = load("05_curves/genus_10_method_2.mrdi")
n_maximal_polyhedra(C.complex)                           # 447
```

## Reproducing and checking

The code repository computes all of `data/` with `make_paper_data.jl` (about 1 h 40 min on a laptop, at most
6 GB per step) and checks every computational statement of the article from it with `paper_results.jl` (about
7 minutes). To check the article from this tarball instead, put it at the root of a clone of the repository and,
from its folder `julia/`, run

    tar -xzf ../TropicalGeometryG2_data.tar.gz --strip-components=1 TropicalGeometryG2_data/data
    julia --project=. paper_results.jl

## Citation

Please cite the article (arXiv:2611.24456) and, for the data, this dataset (doi:10.5281/zenodo.24104003).
EOF

tar --sort=name --owner=0 --group=0 --numeric-owner --mtime="$(date -u +%F)" \
    -czf "$TARBALL" -C "$STAGE" "$NAME"

echo "wrote $TARBALL: $((ALIASES + 1)) files at the top level, $(find "$OUT/data" -type f | wc -l) in data/," \
     "$(du -h "$TARBALL" | cut -f1)"
