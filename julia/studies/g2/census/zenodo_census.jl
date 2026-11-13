###############################################################################
# Stage `zenodo` of the G2 census, the target `zenodo` of run_census.jl: the
# files of the Zenodo tarball (../production/make_zenodo_tarball.sh) that no
# other stage writes. It reads the artifacts of every target of run_census.jl
# and of run_census_gf2.jl, so it runs after both censuses.
#
#   ideals   data/g2/pointspace{,_gf2}/00_ideal.mrdi: the quadric Q, the one
#            generator of the flag ideal in x1, ..., x7 alone, over QQ and GF(2)
#   hexagon  data/g2/hexagon/tropicalized_hexagon.mrdi: the tropicalized hexagon
#            of the Example in Section 9 of main.tex (min-convention), completed
#            from its flags 1 and 4 by tropical_complete_hexagon (about 70 s,
#            2.5 GB) and checked against the table of main.tex
#   tables   the orbit correspondences between the labels of main.tex and the
#            census orbits: data/g2/{grassmannian,flag}_orbits.csv,
#            data/g2/dressian/dressian_orbits.csv and
#            data/g2/flag_dressian/flag_dressian_orbits.csv
#   curve of Method 1, curve of Method 2
#            data/mukai_curves/genus_10_method_{1,2}.mrdi: the two genus-10
#            curves of Section 10 as polyhedral complexes. Method 1 is cut
#            from the raw fan by B + q, from its shift q alone
#            (curves_g10/overnight_20261028/wpart_2722/translation_1.mrdi), by
#            studies/g2/curves/method_1_curve.jl (about 1 min with julia -t 8,
#            the stage then peaks at 3.2 GB), which also writes its 454 pieces
#            to genus_10_method_1_pieces.mrdi; Method 2 is the stored curve of
#            curves_g10/overnight_20261028/experiment_100128; both complexes
#            by studies/g2/curves/curve_complex.jl
#
# A part whose files all exist is skipped unless --force, even when its inputs
# are missing, and its stored files are checked (the hexagon against main.tex,
# the curves for genus 10). The provenance, data/g2/provenance/zenodo.toml, is
# written when a part ran.
###############################################################################

module G2Zenodo

using Oscar
using TropicalGeometryG2

export zenodo_target_definition

module HexagonCompletion
include(joinpath(@__DIR__, "..", "..", "hexagon", "tropical_complete_hexagon.jl"))
end

module CurveComplex
include(joinpath(@__DIR__, "..", "curves", "curve_complex.jl"))
end

module Method1Curve
include(joinpath(@__DIR__, "..", "curves", "method_1_curve.jl"))
end

using .CurveComplex: write_curve, bounded_betti_number

census_path(target, name) = data_path("g2", target, name)

const ZENODO_DIRECTORY = data_path("g2")
const CURVES_G10 = project_path("curves_g10", "overnight_20261028")
const RAY_LABEL_ORBITS_CSV = study_path(:g2, "atlas", "grassmannian", "ray_label_orbits.csv")
const RAY_ORBIT_LABELS_CSV = study_path(:g2, "atlas", "grassmannian", "ray_orbit_labels.csv")
const FLAG_ATLAS_CSV = study_path(:g2, "atlas", "flag", "local_flag_cell_atlas.csv")
const CELL_METRICS_CSV = data_path("g2", "maximal_cell_tree_metrics.csv")

const IDEAL_SOURCES = (
    (census_path("flag", "00_ideal.mrdi"), census_path("pointspace", "00_ideal.mrdi")),
    (census_path("flag_gf2", "00_ideal.mrdi"), census_path("pointspace_gf2", "00_ideal.mrdi")),
)
const HEXAGON_OUTPUT = census_path("hexagon", "tropicalized_hexagon.mrdi")
const TABLE_TARGETS = ("grassmannian", "grassmannian_gf2", "dressian", "flag", "flag_gf2", "flag_dressian")
const TABLE_OUTPUTS = (
    grassmannian=data_path("g2", "grassmannian_orbits.csv"),
    dressian=census_path("dressian", "dressian_orbits.csv"),
    flag=data_path("g2", "flag_orbits.csv"),
    flag_dressian=census_path("flag_dressian", "flag_dressian_orbits.csv"),
)
# Method 1 is recomputed from its shift q, its pieces kept; Method 2 is the
# stored curve of its night-search run, with one multiplicity per piece.
const METHOD_1_INPUTS = [Method1Curve.METHOD_1_SHIFT; collect(Method1Curve.METHOD_1_FILES)]
const METHOD_1_PIECES = data_path("mukai_curves", "genus_10_method_1_pieces.mrdi")
const METHOD_1_OUTPUT = data_path("mukai_curves", "genus_10_method_1.mrdi")
const METHOD_2_SOURCE = joinpath(CURVES_G10, "experiment_100128", "curve.mrdi")
const METHOD_2_OUTPUT = data_path("mukai_curves", "genus_10_method_2.mrdi")

table_inputs() = vcat(
    [census_path(target, name) for target in TABLE_TARGETS
     for name in ("03_normalized_coarsest_fan.mrdi", "05_symmetry_ray_actions.mrdi", "05_symmetry_cell_orbits.csv")],
    [census_path(target, "02_coarsest_fan.mrdi") for target in ("flag", "flag_gf2")],
    [RAY_LABEL_ORBITS_CSV, RAY_ORBIT_LABELS_CSV, CELL_METRICS_CSV, FLAG_ATLAS_CSV],
)

"""
The parts of the stage, each with the files it reads and writes. Every output
of the stage is an output of exactly one part.
"""
zenodo_parts() = (
    (name="ideals", run=export_ideals, check=nothing,
     inputs=[source for (source, _) in IDEAL_SOURCES], outputs=[output for (_, output) in IDEAL_SOURCES]),
    (name="hexagon", run=export_hexagon, check=check_stored_hexagon,
     inputs=[normpath(HexagonCompletion.RAW_FAN_PATH)], outputs=[HEXAGON_OUTPUT]),
    (name="tables", run=export_tables, check=nothing,
     inputs=table_inputs(), outputs=collect(values(TABLE_OUTPUTS))),
    (name="curve of Method 1", run=export_method_1, check=check_stored_method_1,
     inputs=METHOD_1_INPUTS, outputs=[METHOD_1_PIECES, METHOD_1_OUTPUT]),
    (name="curve of Method 2", run=export_method_2, check=() -> check_stored_curve(METHOD_2_OUTPUT),
     inputs=[METHOD_2_SOURCE], outputs=[METHOD_2_OUTPUT]),
)

zenodo_outputs() = reduce(vcat, [part.outputs for part in zenodo_parts()])

###############################################################################
# CSV
###############################################################################

semicolon_integers(field) = isempty(field) ? Int[] : parse.(Int, split(field, ';'))

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            println(io, join((string(row[column]) for column in header), ','))
        end
    end
    println("  wrote $path ($(length(rows)) rows)")
end

###############################################################################
# Ideals of the point space Q
###############################################################################

"""
The ideal of the quadric Q in P^6: the one generator of the flag ideal in the
point coordinates x1, ..., x7 alone, in a ring of these seven variables.
"""
function pointspace_ideal(flag_ideal)
    ring = base_ring(flag_ideal)
    x = gens(ring)[1:7]
    string.(x) == ["x$i" for i in 1:7] || error("the flag ring does not start with x1, ..., x7")
    point_generators = [f for f in gens(flag_ideal) if all(v -> v in x, vars(f))]
    length(point_generators) == 1 || error("expected one generator of the flag ideal in x1, ..., x7 alone")
    point_ring, y = polynomial_ring(coefficient_ring(ring), ["x$i" for i in 1:7])
    quadric = evaluate(only(point_generators), vcat(y, zeros(point_ring, ngens(ring) - 7)))
    total_degree(quadric) == 2 || error("the generator in x1, ..., x7 is not a quadric")
    return ideal(point_ring, [quadric])
end

function export_ideals()
    for (source, output) in IDEAL_SOURCES
        mkpath(dirname(output))
        save(output, pointspace_ideal(load(source)))
        println("  wrote $output")
    end
end

###############################################################################
# The tropicalized hexagon
###############################################################################

# The apples a^(k) in trop(Q) and the trees d^(ij) in trop(G) of the Example in
# Section 9 of main.tex, in the min-convention of the comments there; the
# d_ij are in lexicographic order.
const HEXAGON_POINTS_MIN = [
    [4, 3, 11, 0, 15, 1, 2],
    [6, 3, 3, 0, 5, 1, 2],
    [3, 0, 2, 1, 4, 5, 5],
    [6, 3, 10, 0, 3, 4, 15],
    [6, 3, 3, 0, 3, 4, 6],
    [3, 0, 2, 1, 4, 5, 3],
]
const HEXAGON_LINES_MIN = [
    [4, 4, 1, 6, 2, 3, 3, 8, 5, 13, 6, 0, 13, 1, 2, 2, 6, 3, 3, 4, 4],
    [6, 6, 3, 8, 4, 5, 3, 0, 5, 1, 2, 2, 9, 3, 4, 4, 2, 3, 5, 6, 6],
    [9, 8, 3, 6, 7, 11, 5, 0, 3, 4, 8, 2, 5, 6, 15, 15, 8, 5, 11, 8, 9],
    [7, 6, 4, 7, 8, 9, 3, 13, 18, 9, 6, 0, 3, 4, 13, 13, 6, 3, 9, 6, 7],
    [7, 6, 3, 6, 7, 9, 3, 0, 3, 4, 6, 2, 5, 6, 6, 6, 6, 3, 9, 6, 7],
    [4, 6, 3, 8, 4, 5, 5, 0, 7, 1, 2, 2, 15, 3, 4, 4, 2, 5, 5, 6, 4],
]
# The same vectors as displayed in main.tex, in the max-convention.
const HEXAGON_POINTS_MAX = [
    [11, 12, 4, 15, 0, 14, 13],
    [0, 3, 3, 6, 1, 5, 4],
    [2, 5, 3, 4, 1, 0, 0],
    [9, 12, 5, 15, 12, 11, 0],
    [0, 3, 3, 6, 3, 2, 0],
    [2, 5, 3, 4, 1, 0, 2],
]
const HEXAGON_LINES_MAX = [
    [9, 9, 12, 7, 11, 10, 10, 5, 8, 0, 7, 13, 0, 12, 11, 11, 7, 10, 10, 9, 9],
    [3, 3, 6, 1, 5, 4, 6, 9, 4, 8, 7, 7, 0, 6, 5, 5, 7, 6, 4, 3, 3],
    [6, 7, 12, 9, 8, 4, 10, 15, 12, 11, 7, 13, 10, 9, 0, 0, 7, 10, 4, 7, 6],
    [11, 12, 14, 11, 10, 9, 15, 5, 0, 9, 12, 18, 15, 14, 5, 5, 12, 15, 9, 12, 11],
    [2, 3, 6, 3, 2, 0, 6, 9, 6, 5, 3, 7, 4, 3, 3, 3, 3, 6, 0, 3, 2],
    [11, 9, 12, 7, 11, 10, 10, 15, 8, 14, 13, 13, 0, 12, 11, 11, 13, 10, 10, 9, 11],
]

"""Check that the two conventions of main.tex differ by a sign and a shift."""
function check_hexagon_conventions()
    for (minimum_vectors, maximum_vectors) in ((HEXAGON_POINTS_MIN, HEXAGON_POINTS_MAX),
                                               (HEXAGON_LINES_MIN, HEXAGON_LINES_MAX))
        for (v, w) in zip(minimum_vectors, maximum_vectors)
            length(unique(v .+ w)) == 1 || error("the conventions of main.tex disagree on $v")
        end
    end
end

"""Check that a hexagon is the one of main.tex, every vector with minimum 0."""
function check_hexagon(points, lines)
    points == [QQ.(v) for v in HEXAGON_POINTS_MIN] && lines == [QQ.(v) for v in HEXAGON_LINES_MIN] ||
        error("the hexagon differs from the Example of Section 9: points $points, lines $lines")
end

function export_hexagon()
    check_hexagon_conventions()
    flag(k) = (HEXAGON_POINTS_MIN[k], HEXAGON_LINES_MIN[k])
    seconds = @elapsed completion = HexagonCompletion.tropical_complete_hexagon(flag(1), flag(4))
    all(isempty, completion.no_representative) ||
        error("coordinates of the hexagon have no tropically faithful representative: $(completion.no_representative)")
    check_hexagon(completion.points, completion.lines)
    println("  completed the hexagon from flags 1 and 4 in $(round(seconds; digits=1)) s, as in main.tex")
    mkpath(dirname(HEXAGON_OUTPUT))
    save(HEXAGON_OUTPUT, (points=completion.points, lines=completion.lines))
    println("  wrote $HEXAGON_OUTPUT")
end

function check_stored_hexagon()
    check_hexagon_conventions()
    hexagon = load(HEXAGON_OUTPUT)
    check_hexagon(hexagon.points, hexagon.lines)
end

###############################################################################
# The genus-10 curves
###############################################################################

function export_method_1()
    shift = Method1Curve.METHOD_1_SHIFT
    seconds = @elapsed curve = Method1Curve.method_1_pieces(Vector{Int}(load(shift)))
    println("  cut the raw fan by B + q, q of $shift, in $(round(seconds; digits=1)) s on ",
            "$(Threads.nthreads()) threads: $(length(curve.pieces)) pieces, balanced, of degree 18 in every coordinate")
    mkpath(dirname(METHOD_1_PIECES))
    save(METHOD_1_PIECES, curve.pieces)
    println("  wrote $METHOD_1_PIECES")
    # The complex alone, as section_A.jl of paper_results reads it: the
    # multiplicities of the pieces are not saved.
    write_curve(METHOD_1_OUTPUT, curve.pieces, nothing)
    check_stored_method_1()
end

function export_method_2()
    println("  loading $METHOD_2_SOURCE")
    stored = load(METHOD_2_SOURCE)
    write_curve(METHOD_2_OUTPUT, stored.polyhedra, stored.multiplicities)
    check_stored_curve(METHOD_2_OUTPUT)
end

"""Print the counts of the appendix of main.tex for the stored curve `output`, and check its genus."""
function check_stored_curve(output)
    complex = load(output).complex
    cells = collect(maximal_polyhedra(complex))
    bounded = filter(is_bounded, cells)
    vertices = unique([collect(v) for P in bounded for v in Oscar.vertices(P)])
    genus = bounded_betti_number(complex)
    genus == 10 || error("the curve of $output has genus $genus, not 10")
    println("  $(basename(output)): $(length(cells)) maximal cells, $(length(cells) - length(bounded)) ",
            "unbounded; the bounded cells form a graph of genus $genus on $(length(vertices)) vertices")
end

"""Check the stored curve of Method 1 and that its stored pieces are one-dimensional."""
function check_stored_method_1()
    pieces = load(METHOD_1_PIECES)
    all(P -> dim(P) == 1, pieces) || error("a piece of $METHOD_1_PIECES is not one-dimensional")
    println("  $(basename(METHOD_1_PIECES)): $(length(pieces)) one-dimensional pieces")
    check_stored_curve(METHOD_1_OUTPUT)
end

###############################################################################
# Orbit correspondences
###############################################################################

# The labels of main.tex of the ray orbits of ray_label_orbits.csv, by its
# orbit number, and of orbit 12, the rays that occur only over GF(2), and
# their order, that of the ray catalog: ray_orbit_labels.csv, which
# to_tex/build_styled_ray_orbit_tables.jl reads as well.
const RAY_LABEL_ROWS = sort!(csv_rows(RAY_ORBIT_LABELS_CSV); by=row -> parse(Int, row["orbit"]))
const RAY_LABELS = Tuple(row["label"] for row in RAY_LABEL_ROWS)
const RAY_LABEL_ORDER = Tuple(row["label"] for row in sort(RAY_LABEL_ROWS; by=row -> parse(Int, row["catalog_position"])))

"""
The coarsest fan of a census target: its rays modulo lineality, its lineality
and its maximal cells, in the numbering of the census, with the orbits of its
rays and maximal cells under the Weyl group.
"""
function census_fan(target)
    fan = load(census_path(target, "03_normalized_coarsest_fan.mrdi"))
    data = rays_modulo_lineality(fan)
    rays = [QQ.(collect(ray)) for ray in data.rays_modulo_lineality]
    actions = Matrix{Int}(load(census_path(target, "05_symmetry_ray_actions.mrdi")))
    size(actions, 2) == length(rays) || error("the ray actions of $target do not act on its rays")
    ray_orbits = unique!([sort!(unique(actions[:, ray])) for ray in eachindex(rays)])
    sort!(reduce(vcat, ray_orbits)) == collect(eachindex(rays)) ||
        error("the ray orbits of $target do not partition its rays")
    cell_orbits = [(orbit=parse(Int, row["orbit"]), members=semicolon_integers(row["cells"]))
                   for row in csv_rows(census_path(target, "05_symmetry_cell_orbits.csv"))]
    return (target=target, rays=rays, lineality=matrix(QQ, data.lineality_basis),
            cells=TropicalGeometryG2.fan_cells(fan), ray_orbits=sort!(ray_orbits; by=first),
            cell_orbits=sort!(cell_orbits; by=row -> row.orbit))
end

"""The projection modulo the lineality of `reference`, shared by `others`."""
function common_quotient(reference, others...)
    for other in others
        rank(reference.lineality) == rank(other.lineality) ==
            rank(vcat(reference.lineality, other.lineality)) ||
            error("$(other.target) has another lineality space than $(reference.target)")
    end
    return Matrix{QQFieldElem}(kernel(reference.lineality; side=:right))
end

ray_signatures(fan, quotient) = [TropicalGeometryG2.quotient_ray_signature(ray, quotient) for ray in fan.rays]

"""Every maximal cell of `fan` as the set of the signatures of its rays."""
cell_keys(fan, signatures) = [Set(signatures[cell]) for cell in fan.cells]

"""The orbit numbers of the maximal cells of `fan`, by cell."""
function cell_orbit_of(fan)
    orbit = zeros(Int, length(fan.cells))
    for row in fan.cell_orbits, cell in row.members
        orbit[cell] = row.orbit
    end
    all(>(0), orbit) || error("a maximal cell of $(fan.target) lies in no orbit")
    return orbit
end

orbit_row(fan, number) = fan.cell_orbits[findfirst(row -> row.orbit == number, fan.cell_orbits)]

members_field(members) = join(sort(members), ';')

"""
The labels of main.tex of the rays of the Grassmannian over QQ and over GF(2),
as a dictionary from the signature of a ray to its label.
"""
function grassmannian_ray_labels(qq, gf2, quotient)
    qq_signatures = ray_signatures(qq, quotient)
    labels = Dict{Any,String}()
    for row in csv_rows(RAY_LABEL_ORBITS_CSV)
        rays = sort!(semicolon_integers(row["rays"]))
        rays in qq.ray_orbits || error("orbit $(row["orbit"]) of ray_label_orbits.csv is no census orbit")
        for ray in rays
            labels[qq_signatures[ray]] = RAY_LABELS[parse(Int, row["orbit"])]
        end
    end
    length(labels) == length(qq.rays) || error("ray_label_orbits.csv does not cover the rays over QQ")
    gf2_only = [signature for signature in ray_signatures(gf2, quotient) if !haskey(labels, signature)]
    length(gf2_only) == 3 || error("expected the three rays $(RAY_LABELS[12]) over GF(2), orbit 12")
    foreach(signature -> labels[signature] = RAY_LABELS[12], gf2_only)
    return labels
end

"""The rows of the ray orbits of `fans`, matched by signature, with their labels."""
function ray_rows(fans, signatures, label_of)
    rows = Dict{Any,Dict{String,Any}}()
    for (field, fan) in fans
        for orbit in fan.ray_orbits
            key = Set(signatures[field][orbit])
            row = get!(rows, key) do
                Dict{String,Any}("kind" => "ray", "orbit_type" => label_of(first(key)),
                                 "size" => length(orbit))
            end
            row["$(field)_orbit"] = minimum(orbit)
            row["$(field)_members"] = members_field(orbit)
        end
    end
    return collect(values(rows))
end

function fill_missing!(rows, header)
    for row in rows, column in header
        haskey(row, column) || (row[column] = "")
    end
    return rows
end

label_rank(label) = something(findfirst(==(label), RAY_LABEL_ORDER), length(RAY_LABEL_ORDER) + 1)

"""
Write the four orbit tables, to the files of `outputs` (by default those of the
stage, `TABLE_OUTPUTS`).
"""
function export_tables(outputs=TABLE_OUTPUTS)
    grassmannian = (qq=census_fan("grassmannian"), gf2=census_fan("grassmannian_gf2"))
    dressian = census_fan("dressian")
    flag = (qq=census_fan("flag"), gf2=census_fan("flag_gf2"))
    flag_dressian = census_fan("flag_dressian")

    # The Grassmannian over QQ and over GF(2).
    quotient = common_quotient(grassmannian.qq, grassmannian.gf2, dressian)
    labels = grassmannian_ray_labels(grassmannian.qq, grassmannian.gf2, quotient)
    signatures = Dict(field => ray_signatures(fan, quotient) for (field, fan) in pairs(grassmannian))
    header = ["kind", "orbit_type", "size", "qq_orbit", "gf2_orbit", "qq_members", "gf2_members"]
    rows = ray_rows(pairs(grassmannian), signatures, signature -> labels[signature])
    sort!(rows; by=row -> label_rank(row["orbit_type"]))
    qq_number = Dict{Int,String}()
    for metric in csv_rows(CELL_METRICS_CSV)
        row = Dict{String,Any}("kind" => "cell", "orbit_type" => metric["number"], "size" => metric["orbit_size"])
        for field in (:qq, :gf2)
            orbit = metric["$(field)_orbit"]
            isempty(orbit) && continue
            members = orbit_row(grassmannian[field], parse(Int, orbit)).members
            length(members) == parse(Int, metric["orbit_size"]) ||
                error("cell $(metric["number"]) has another orbit size over $field")
            row["$(field)_orbit"] = orbit
            row["$(field)_members"] = members_field(members)
        end
        isempty(metric["qq_orbit"]) || (qq_number[parse(Int, metric["qq_orbit"])] = metric["number"])
        push!(rows, row)
    end
    write_csv(outputs.grassmannian, header, fill_missing!(rows, header))

    # The Dressian, by the rays and maximal cells over QQ.
    dressian_signatures = ray_signatures(dressian, quotient)
    qq_cell = Dict(key => cell for (cell, key) in
                   enumerate(cell_keys(grassmannian.qq, signatures[:qq])))
    qq_orbit = cell_orbit_of(grassmannian.qq)
    dressian_keys = cell_keys(dressian, dressian_signatures)
    header = ["kind", "orbit_type", "size", "orbit", "members"]
    rows = ray_rows([(:dressian, dressian)], Dict(:dressian => dressian_signatures),
                    signature -> labels[signature])
    sort!(rows; by=row -> label_rank(row["orbit_type"]))
    for row in rows
        row["orbit"], row["members"] = row["dressian_orbit"], row["dressian_members"]
    end
    for orbit in dressian.cell_orbits
        key = dressian_keys[first(orbit.members)]
        label = haskey(qq_cell, key) ? qq_number[qq_orbit[qq_cell[key]]] : "epsilon"
        label == "epsilon" && length(key) != 4 && error("Dressian orbit $(orbit.orbit) is no tetrahedron")
        push!(rows, Dict{String,Any}("kind" => "cell", "orbit_type" => label, "size" => length(orbit.members),
                                     "orbit" => orbit.orbit, "members" => members_field(orbit.members)))
    end
    write_csv(outputs.dressian, header, rows)

    # The flag variety over QQ and over GF(2): a flag ray is labelled by the
    # Grassmannian ray below it, or star over the lineality space.
    flag_quotient = common_quotient(flag.qq, flag.gf2, flag_dressian)
    flag_signatures = Dict(field => ray_signatures(fan, flag_quotient) for (field, fan) in pairs(flag))
    flag_label(ray) = all(iszero, transpose(ray[8:28]) * quotient) ? "star" :
        labels[TropicalGeometryG2.quotient_ray_signature(ray[8:28], quotient)]
    flag_ray_label = Dict{Any,String}()
    for (field, fan) in pairs(flag), (ray, signature) in zip(fan.rays, flag_signatures[field])
        flag_ray_label[signature] = flag_label(ray)
    end
    header = ["kind", "orbit_type", "apple", "size", "qq_orbit", "gf2_orbit", "qq_members", "gf2_members"]
    rows = ray_rows(pairs(flag), flag_signatures, signature -> flag_ray_label[signature])
    sort!(rows; by=row -> (label_rank(row["orbit_type"]), get(row, "qq_orbit", typemax(Int))))
    # Maximal cells: the edge-decorated cells of the flag atlas, by the rays of
    # the atlas fan of each field, 02_coarsest_fan.mrdi.
    atlas_fans = Dict(field => let fan = load(census_path(flag[field].target, "02_coarsest_fan.mrdi"))
                          [TropicalGeometryG2.quotient_ray_signature(QQ.(collect(ray)), flag_quotient)
                           for ray in rays_modulo_lineality(fan).rays_modulo_lineality]
                      end for field in (:qq, :gf2))
    census_cell = Dict(field => Dict(key => cell for (cell, key) in
                                     enumerate(cell_keys(flag[field], flag_signatures[field])))
                       for field in (:qq, :gf2))
    census_orbit = Dict(field => cell_orbit_of(flag[field]) for field in (:qq, :gf2))
    cell_rows = Dict{Any,Dict{String,Any}}()
    for atlas in csv_rows(FLAG_ATLAS_CSV)
        atlas["location_kind"] == "edge" || continue
        orbits = Dict{Symbol,Int}()
        for field in (:qq, :gf2)
            string(field) in split(atlas["fields"], '+') || continue
            key = Set(atlas_fans[field][semicolon_integers(atlas["$(field)_flag_rays"])])
            orbits[field] = census_orbit[field][census_cell[field][key]]
        end
        isempty(orbits) && continue
        location = atlas["location"]
        m = match(r"^\(:ray, (\d+)\)$", location)
        apple = m !== nothing ? "pendant $(m[1])" :
            "inner " * join(match(r"^\(:edge, (\d+), (\d+)\)$", location).captures, '-')
        key = (atlas["label"], get(orbits, :qq, 0), get(orbits, :gf2, 0))
        row = get!(cell_rows, key) do
            sizes = unique([length(orbit_row(flag[field], orbit).members) for (field, orbit) in orbits])
            length(sizes) == 1 || error("atlas label $(atlas["label"]) has orbits of different sizes")
            row = Dict{String,Any}("kind" => "cell", "orbit_type" => atlas["label"], "size" => only(sizes),
                                   "apples" => String[])
            for (field, orbit) in orbits
                row["$(field)_orbit"] = orbit
                row["$(field)_members"] = members_field(orbit_row(flag[field], orbit).members)
            end
            row
        end
        push!(row["apples"], apple)
    end
    for field in (:qq, :gf2)
        listed = sort!([orbit for (_, qq, gf2) in keys(cell_rows) for orbit in (field == :qq ? qq : gf2) if orbit > 0])
        listed == [row.orbit for row in flag[field].cell_orbits] ||
            error("the flag atlas does not list every orbit over $field once")
    end
    for row in values(cell_rows)
        row["apple"] = join(sort!(unique(pop!(row, "apples"))), '|')
    end
    cell_list = sort!(collect(values(cell_rows));
                      by=row -> (parse(Int, row["orbit_type"]), get(row, "qq_orbit", typemax(Int)),
                                 get(row, "gf2_orbit", typemax(Int))))
    flag_rows = fill_missing!(vcat(rows, cell_list), header)
    write_csv(outputs.flag, header, flag_rows)

    # The flag Dressian, by the rays and maximal cells over QQ; the cells that
    # are not cells of trop(F) lie over the tetrahedra epsilon.
    flag_dressian_signatures = ray_signatures(flag_dressian, flag_quotient)
    qq_flag_row = Dict(row["qq_orbit"] => row for row in cell_list if row["qq_orbit"] != "")
    header = ["kind", "orbit_type", "apple", "size", "orbit", "members"]
    rows = ray_rows([(:flag_dressian, flag_dressian)], Dict(:flag_dressian => flag_dressian_signatures),
                    signature -> flag_ray_label[signature])
    sort!(rows; by=row -> (label_rank(row["orbit_type"]), row["flag_dressian_orbit"]))
    for row in rows
        row["orbit"], row["members"], row["apple"] = row["flag_dressian_orbit"], row["flag_dressian_members"], ""
    end
    flag_dressian_keys = cell_keys(flag_dressian, flag_dressian_signatures)
    for orbit in flag_dressian.cell_orbits
        key = flag_dressian_keys[first(orbit.members)]
        cell = get(census_cell[:qq], key, 0)
        source = cell > 0 ? qq_flag_row[census_orbit[:qq][cell]] : nothing
        push!(rows, Dict{String,Any}(
            "kind" => "cell", "orbit_type" => isnothing(source) ? "epsilon" : source["orbit_type"],
            "apple" => isnothing(source) ? "" : source["apple"], "size" => length(orbit.members),
            "orbit" => orbit.orbit, "members" => members_field(orbit.members)))
    end
    write_csv(outputs.flag_dressian, header, rows)
end

###############################################################################
# The stage
###############################################################################

function stage_zenodo(target, options)
    parts = zenodo_parts()
    provenance = stage_provenance_path(ZENODO_DIRECTORY, zenodo)
    ran = false
    for part in parts
        # Stored outputs are reused before any input is looked for.
        if all(isfile, part.outputs) && !options.force
            println("  $(part.name): using existing canonical artifact(s)")
            isnothing(part.check) || part.check()
            continue
        end
        missing_inputs = filter(!isfile, part.inputs)
        if !isempty(missing_inputs)
            hint = any(contains("_gf2"), missing_inputs) ?
                "; run studies/g2/run_census_gf2.jl first" : ""
            error("the part $(part.name) needs the missing files\n  " * join(missing_inputs, "\n  ") * hint)
        end
        println("  $(part.name)")
        part.run()
        ran = true
    end
    # Reused outputs keep their provenance: their inputs need not be present.
    ran || !isfile(provenance) || return nothing
    write_provenance(
        provenance;
        study="G2", target=target.name, stage=zenodo,
        inputs=unique(reduce(vcat, [part.inputs for part in parts])),
        outputs=zenodo_outputs(),
        configuration=census_configuration(options),
        metadata=Dict{String,Any}(
            "tarball" => "production/make_zenodo_tarball.sh",
            "pointspace_ideal" => "the one generator of the flag ideal in x1, ..., x7 alone",
            "hexagon" => "tropical_complete_hexagon from flags 1 and 4, equal to the Example of Section 9",
            "curves" => "Method 1 cut from the raw fan by B + q, q of wpart 2722, by curves/method_1_curve.jl; " *
                        "Method 2 the stored curve of experiment 100128",
        ),
    )
end

"""The target `zenodo` of the G2 census: the stage `zenodo` alone."""
zenodo_target_definition() = StudyTarget(
    name="zenodo",
    stages=Set([zenodo]),
    callbacks=Dict{CensusStage,Function}(zenodo => stage_zenodo),
    artifacts=Dict(zenodo => zenodo_outputs()),
)

end # module G2Zenodo
