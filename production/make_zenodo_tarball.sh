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
# Data for the paper on the tropical G2-Grassmannian and Mukai curves

The `.mrdi` files are written by OSCAR (1.8.2 for most of them; the version
is in the header of each file) and are read with `Oscar.load`.
Tropical objects use the min-convention. Q, G and F are the point space
(the quadric in P^6), the G2-Grassmannian G2/P2 in P^20 and the flag variety
G2/B in P^6 x P^20. The suffixes `qq` and `gf2` name the field of coefficients:
QQ or GF(2). Over QQ the fans are those of every characteristic other than 2.
Figures and tables of the paper are cited by their labels in main.tex.

- `*_orbits.csv`: the orbit correspondences between the orbit types of the
  paper and the census orbits of the fans; see below.
- `00_ideals/`: the ideals of Q, G and F over QQ and GF(2).
- `01_dressians/`: the G2-Dressian and the G2 flag Dressian, as coarsest fans.
  They are the same over every field.
- `02_tropical_varieties/`: trop(G) and trop(F) over QQ and GF(2), as coarsest fans.
- `03_no_polytopes/`: the Newton-Okounkov polytopes of the toric degenerations
  (`tab:g2-grassmannian-newton-okounkov`), one per orbit of prime cones.
  `orbit_NNN` is the census orbit NNN of the maximal cells of the fan of the
  same name in `02_tropical_varieties/`.
- `04_hexagons/tropicalized_hexagon.mrdi`: the tropicalized hexagon of
  `fig:tropical_hexagon`, a named tuple with `points` (the apples a^(1..6) in
  trop(Q), in R^7) and `lines` (the trees d^(12), d^(23), ..., d^(61) in
  trop(G), in R^21, with the coordinates p_ij in lexicographic order).
- `05_curves/`: the two genus-10 curves as polyhedral complexes in R^21:
  Method 1 (`fig:planar_curve`) as `(complex,)`, and Method 2
  (`fig:nonplanar_curve`) as `(complex, multiplicities)` with one
  multiplicity per maximal cell. `genus_10_method_1_pieces.mrdi` is the
  vector of the 454 polyhedra computed for Method 1 from its shift q, which
  `genus_10_method_1` is assembled from.
- `data/`: the folder `data/` that `make_paper_data.jl` writes in the folder
  `julia/` of the code repository: every stage of every computation, from
  which the files above are copied and from which `paper_results.jl` checks
  the statements of the paper. See `data/README.md`.

## Using the data with the code

The code is at https://github.com/dmg-lab/TropicalGeometryG2. Instead of running
`make_paper_data.jl` (about 1.5 hours), put the tarball at the root of a clone
and unpack `data/` into its folder `julia/`; from `julia/`:

    tar -xzf ../TropicalGeometryG2_data.tar.gz --strip-components=1 TropicalGeometryG2_data/data
    julia --project=. paper_results.jl

## Orbit correspondences

The census numbers the rays and maximal cells of each fan, and their orbits
under the Weyl group D12, by machine. The four `*_orbits.csv` files translate
these orbits to the orbit types of the paper. One line is one orbit; the
orbits of the two fields share a line when they consist of the same rays or
cells. The member indices refer to these fans:

| file                       | fans                                               |
|----------------------------|----------------------------------------------------|
| `dressian_orbits.csv`      | `01_dressians/dressian.mrdi`                       |
| `flag_dressian_orbits.csv` | `01_dressians/flag_dressian.mrdi`                  |
| `grassmannian_orbits.csv`  | `02_tropical_varieties/grassmannian_{qq,gf2}.mrdi` |
| `flag_orbits.csv`          | `02_tropical_varieties/flag_{qq,gf2}.mrdi`         |

For example, `grassmannian_orbits.csv` has the lines

    kind,orbit_type,size,qq_orbit,gf2_orbit,qq_members,gf2_members
    ray,alpha_2,3,8,1,8;33;34,1;19;22

The rays of orbit type alpha_2 form an orbit of size 3. They are the rays
8, 33 and 34 of `grassmannian_qq.mrdi` and the rays 1, 19 and 22 of
`grassmannian_gf2.mrdi`.

The columns:

- `kind`: `ray` or `cell`.
- `orbit_type`: the orbit type in the paper. For rays, alpha_1, ..., epsilon
  of `fig:atlas_rays`. For cells, the types 1-35 of
  `fig:atlas_toric_cells_prime` and `fig:atlas_toric_cells_nontoric`. For
  the flag fans, a ray gets the orbit type of the Grassmannian ray below it,
  or `star` over the lineality space. In the flag Dressian, the cells over the
  tetrahedra of `fig:tetrahedra` have orbit type epsilon.
- `apple` (flag fans): the edges of the tree of the cell that carry the apple
  (`fig:decorated_trees`). `pendant k` is the pendant edge of leaf k.
  `inner a-b` is the inner edge between the nodes a and b of the tree.
- `size`: the orbit size.
- `orbit` (Dressians), `qq_orbit`, `gf2_orbit`: the census orbit. For a cell,
  the orbit number of the census; in `qq_orbit` and `gf2_orbit` of a toric
  cell it is the NNN of `03_no_polytopes/`. For a ray, the smallest ray index
  in the orbit. Empty if the orbit does not occur over that field.
- `members`, `qq_members`, `gf2_members`: the indices of the rays or maximal
  cells in the orbit. Rays are numbered as in `rays_modulo_lineality`, and
  maximal cells as in `maximal_cones` of the fan.
EOF

tar --sort=name --owner=0 --group=0 --numeric-owner --mtime="$(date -u +%F)" \
    -czf "$TARBALL" -C "$STAGE" "$NAME"

echo "wrote $TARBALL: $((ALIASES + 1)) files at the top level, $(find "$OUT/data" -type f | wc -l) in data/," \
     "$(du -h "$TARBALL" | cut -f1)"
