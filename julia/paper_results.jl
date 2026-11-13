#=
paper_results.jl

The computational statements of the paper "Tropical Geometry of G2" (Joswig, Kuehn, Sturmfels, Vallée), in the order
in which they appear, computed from the data stored in data/ and from the certificates and reports stored next to
the scripts. For each statement the script prints what it computes and compares it with the paper:

    [ ok   ]  computed here, equal to the paper
    [FAIL  ]  computed here, different from the paper (the paper's value follows)
    [warn  ]  computed here, different from the paper, left as it is by the authors (not a failure)
    [stored]  read from a stored certificate or report, written by the script named after it
    [cited ]  a statement of the paper that cites the literature
    [skip  ]  not run (--quick, or an input that only a longer script writes)

Tables that the paper describes in words only (for instance the toric degenerations of the flag variety) are printed
as well. A summary ends the output; the exit status is 1 when a statement fails.

The helpers are in this file; the statements of each section are in paper_results/section_N.jl (a function
section_N()), loaded only when the section is run.

Usage, from julia/:

    julia --project=. paper_results.jl [--sections=1,2,...,10,A] [--quick]

--quick skips the slower checks: in_{w0}(J) recomputed (Section 3); the ±1 coefficients of the stored Gröbner bases
and the quadrilaterals that are two Gröbner cones (Section 4); in_C(G) at a second point of every prime cone and the
unimodular equivalence of the Newton-Okounkov polytopes over QQ and GF(2) (Section 6); the 23 fiber complexes of
Table tab:local_flag_complex (Section 7); the random lifts over QQ(t) and the completion of the hexagon (Section 9);
the Betti table (Section 10). The full run takes about 7 minutes and 4.5 GB, --quick about 5 minutes and 3.6 GB.
=#

using TropicalGeometryG2
using Oscar
using Printf
using Random
import Combinatorics
using TOML
import TropicalGeometryG2: cell_faces, fan_cells, load_symmetry_artifacts

const QUICK = "--quick" in ARGS
const SECTION_OPTION = [String(split(a, '='; limit=2)[2]) for a in ARGS if startswith(a, "--sections=")]
const ALL_SECTIONS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "A"]
const SECTIONS = isempty(SECTION_OPTION) ? ALL_SECTIONS : String.(split(only(SECTION_OPTION), ','))
for option in ARGS
    option == "--quick" || startswith(option, "--sections=") ||
        error("unknown option $option; the options are --sections=1,2,...,10,A and --quick")
end
all(in(ALL_SECTIONS), SECTIONS) ||
    error("unknown section in --sections: $(join(setdiff(SECTIONS, ALL_SECTIONS), ", ")); the sections are $(join(ALL_SECTIONS, ","))")
const REPOSITORY = normpath(joinpath(@__DIR__, ".."))
repository(parts...) = joinpath(REPOSITORY, parts...)
g2(parts...) = data_path("g2", parts...)

###############################################################################
# Output and bookkeeping
###############################################################################

const COUNTS = Dict(tag => 0 for tag in ("ok", "FAIL", "warn", "stored", "cited", "skip"))
const FAILED = String[]
const WARNED = String[]

show_value(x::Union{Tuple,AbstractVector}) = "(" * join(show_value.(x), ", ") * ")"
show_value(x) = string(x)

function heading(number, title)
    println("\n", "="^118, "\n", number, ". ", title, "\n", "="^118)
end

subheading(title) = println("\n  -- ", title)

"Print a statement computed here and compare it with the value of the paper."
function check(statement, paper, computed)
    ok = isequal(paper, computed)
    tag = ok ? "ok" : "FAIL"
    COUNTS[tag] += 1
    ok || push!(FAILED, statement)
    println("  [", rpad(tag, 6), "] ", statement, ": ", show_value(computed),
            ok ? "" : "   <- the paper says $(show_value(paper))")
    flush(stdout)
    return ok
end

"""
As `check`, for a difference from the paper that the authors decided to leave: printed as a warning, which does not
count as a failure and does not change the exit status.
"""
function warned(statement, paper, computed)
    ok = isequal(paper, computed)
    tag = ok ? "ok" : "warn"
    COUNTS[tag] += 1
    ok || push!(WARNED, statement)
    println("  [", rpad(tag, 6), "] ", statement, ": ", show_value(computed),
            ok ? "" : "   <- the paper says $(show_value(paper)) (left as it is)")
    flush(stdout)
    return ok
end

"A statement that holds when `condition` does, printed with what was computed."
claim(statement, condition::Bool, detail="") = check(statement * (isempty(detail) ? "" : " ($detail)"), true, condition)

"""
Print a value read from a stored certificate or report. With `paper`, compare it with the value of the paper: a
difference is a FAIL, as for `check`.
"""
function stored(statement, value; paper=nothing)
    if paper !== nothing && !isequal(paper, value)
        COUNTS["FAIL"] += 1
        push!(FAILED, statement)
        println("  [FAIL  ] ", statement, " (stored): ", show_value(value), "   <- the paper says $(show_value(paper))")
    else
        COUNTS["stored"] += 1
        println("  [stored] ", statement, ": ", paper === nothing ? value : show_value(value))
    end
    flush(stdout)
end

function cited(statement)
    COUNTS["cited"] += 1
    println("  [cited ] ", statement)
end

function skipped(statement, reason)
    COUNTS["skip"] += 1
    println("  [skip  ] ", statement, " (", reason, ")")
end

"""
Whether the repository has the figure and table files `files` of the paper (paths from the repository root, the parent
of julia/, such as "figures/ray_tree.tex"); if not, print `statement`, the checks that read them, as skipped. The
public repository has no figures and no tables.
"""
function drawings(statement, files...)
    absent = [file for file in files if !isfile(repository(file))]
    isempty(absent) || skipped(statement, "$(join(absent, ", ")): the figures and tables are not in this repository")
    return isempty(absent)
end

"The text of a figure file of the paper, figures/`name`.tex."
figure_text(name) = read(repository("figures", "$name.tex"), String)

info(text) = (println("           ", text); flush(stdout))

"Run one section (`f` loads and runs it), timing it, and report an error without stopping the others."
function run_section(f, number, title)
    string(number) in SECTIONS || return
    heading(number, title)
    seconds = @elapsed try
        Base.invokelatest(f)
    catch error
        COUNTS["FAIL"] += 1
        push!(FAILED, "section $number raised an error")
        println("  [FAIL  ] the section raised an error: ", sprint(showerror, error))
        for frame in stacktrace(catch_backtrace())[1:min(end, 6)]
            println("           ", frame)
        end
    end
    @printf("\n  (section %s: %.0f s; peak memory of the run so far %.1f GB)\n", number, seconds, Sys.maxrss() / 2^30)
    flush(stdout)
end

###############################################################################
# Stored data
###############################################################################

"The rows of a CSV file as dictionaries (`csv_rows`); fields in double quotes may contain commas, and lines starting with # are skipped."
rows(path) = csv_rows(path; comments=true)

const TARGET_DIRECTORIES = Dict(
    "dressian" => g2("dressian"), "flag_dressian" => g2("flag_dressian"),
    "grassmannian" => g2("grassmannian"), "grassmannian_gf2" => g2("grassmannian_gf2"),
    "flag" => g2("flag"), "flag_gf2" => g2("flag_gf2"),
    "slong" => g2("grassmannian_Slong"), "short_long" => g2("grassmannian_Slong", "short_long"),
)

const SUMMARIES = Dict{String,Any}()

"The statements of `census_summary` about a census target, computed from its stored fan and symmetry data."
function summary(target)
    get!(SUMMARIES, target) do
        directory = TARGET_DIRECTORIES[target]
        stage1 = joinpath(directory, isfile(joinpath(directory, "01_raw_fan.mrdi")) ? "01_raw_fan.mrdi" :
                                     "01_prevariety_fan.mrdi")
        census_summary(directory; stage1_fan=stage1)
    end
end

const FANS = Dict{String,Any}()
"The normalized coarsest fan of a census target (the fans of the Zenodo tarball)."
fan(target) = get!(() -> load(joinpath(TARGET_DIRECTORIES[target], "03_normalized_coarsest_fan.mrdi")), FANS, target)

"The rays of a fan modulo lineality, as rational vectors, and a basis of its lineality space."
function fan_rays(F)
    data = rays_modulo_lineality(F)
    lineality = Matrix{QQFieldElem}(matrix(QQ, data.lineality_basis))
    return [QQ.(collect(r)) for r in data.rays_modulo_lineality], [lineality[k, :] for k in 1:size(lineality, 1)]
end

"""
All faces of a fan modulo lineality, by orbits under its symmetry group: a dictionary from one face of each orbit (a
sorted list of rays) to its dimension and the size of its orbit. As in `census_summary`, only one maximal cell per
orbit is examined.
"""
face_orbits(target) = first(TropicalGeometryG2.face_orbits(fan(target),
    load_symmetry_artifacts(joinpath(TARGET_DIRECTORIES[target], "05_symmetry"))))

"The number of faces of dimension `d` with `n` vertices, from `face_orbits`."
face_count(orbits, d, n) = sum((size for (face, (k, size)) in orbits if k == d && length(face) == n); init=0)

"The paper's f-vector of a census target: the faces modulo lineality, counted by dimension."
f_vector_of(target) = summary(target)["complex"]["f_vector"]

reduced_homology(target) = summary(target)["reduced_homology"]["description"]

maximal_types(target) = Dict((t["dimension"], t["type"]) => (t["cells"], t["orbits"]) for t in summary(target)["maximal_cell_types"])

"""
Locating points in a fan: for each maximal cone its equations and facet inequalities, so that `locate(L, v)` returns the
maximal cones containing `v` in their relative interior.
"""
function locator_entry(C)
    equations = Matrix{QQFieldElem}(linear_equation_matrix(linear_span(C)))
    inequalities = Matrix{QQFieldElem}(linear_inequality_matrix(facets(C)))
    approximate = Float64.(vcat(equations, inequalities))
    return (equations=equations, facets=inequalities, approximate=approximate,
            scale=vec(sum(abs, approximate; dims=2)), equation_count=size(equations, 1))
end
locator(F) = [locator_entry(C) for C in maximal_cones(F)]

"""
The cones of `L` that may contain `v` (closed): a floating-point test with a tolerance far above the rounding error, so it
never drops a cone that contains `v`. The exact test runs on these only, which saves most of the rational arithmetic.
"""
function candidate_cones(L, v)
    x = Float64.(collect(v))
    size_of_x = maximum(abs, x; init=0.0)
    possible(c) = !hasproperty(c, :approximate) || begin
        y = c.approximate * x
        all(i -> (i <= c.equation_count ? abs(y[i]) : y[i]) <= 1e-6 * c.scale[i] * size_of_x + 1e-9, eachindex(y))
    end
    return [k for (k, c) in enumerate(L) if possible(c)]
end
locate(L, v) = filter(k -> all(iszero, L[k].equations * v) && all(<(0), L[k].facets * v), candidate_cones(L, v))
locate_closed(L, v) = filter(k -> all(iszero, L[k].equations * v) && all(<=(0), L[k].facets * v), candidate_cones(L, v))

"The orbit types of the paper of the rays and maximal cells, by census index, from the orbit tables of the stage zenodo."
function orbit_types(file, field)
    types = Dict{Tuple{String,Int},String}()
    for row in rows(g2(file))
        members = row["$(field)_members"]
        isempty(members) && continue
        for m in split(members, ';')
            types[(row["kind"], parse(Int, m))] = row["orbit_type"]
        end
    end
    return types
end

const PAIRS = [(i, j) for i in 1:7 for j in i+1:7]
const PAIR_INDEX = Dict(p => k for (k, p) in enumerate(PAIRS))
pair_name(k) = "p$(PAIRS[k][1])$(PAIRS[k][2])"

rank_of(vectors) = isempty(vectors) ? 0 : rank(matrix(QQ, reduce(vcat, permutedims.(vectors))))

"A name of a module defined by `run_in_module`, read in the latest world."
latest(m, name) = Base.invokelatest(getproperty, m, name)
"Call a function of such a module."
call(m, name, args...) = Base.invokelatest(latest(m, name), args...)

"Run a script of the repository in a fresh module, so that its top-level names do not clash; returns the module."
function run_in_module(path; name=:Script)
    m = Module(name)
    Core.eval(m, :(include(x) = Base.include($m, x)))
    Base.include(m, path)
    return m
end

"The toric orbits of a stage-6 census, with the paper's type of each, and its Newton-Okounkov polytope."
function toric_orbits(target, field)
    census = rows(g2(target, "06_canonical_interior_degenerations.csv"))
    type_of = Dict{String,String}()
    file = startswith(target, "grassmannian") ? "grassmannian_orbits.csv" : "flag_orbits.csv"
    apple = Dict{String,String}()
    for r in rows(g2(file))
        r["kind"] == "cell" && !isempty(r["$(field)_orbit"]) || continue
        type_of[r["$(field)_orbit"]] = r["orbit_type"]
        apple[r["$(field)_orbit"]] = get(r, "apple", "")
    end
    return [(row=r, type=type_of[r["orbit"]], apple=apple[r["orbit"]]) for r in census if r["toric"] == "true"]
end

body(target, row) = load(joinpath(g2(target), row["newton_okounkov_body"]))

const PROJECTIONS = Dict{String,Bool}()

"""
Whether the exact projection certificate 07_projection_certificate.toml of a flag census target (flag or flag_gf2) is
valid: the images in trop(G) of the same field of its rays and maximal cells, and their D12-transport. Verified once per
target, for Sections 3, 4 and A; a verification that raises an error is invalid.
"""
function projection_certificate_valid(target)
    get!(PROJECTIONS, target) do
        grassmannian = replace(target, "flag" => "grassmannian")
        try
            verify_ray_projection_certificate(TOML.parsefile(g2(target, "07_projection_certificate.toml")),
                fan(target), load(g2(target, "03_normalized_rays.mrdi")),
                fan(grassmannian), load(g2(grassmannian, "03_normalized_rays.mrdi")),
                load_symmetry_artifacts(g2(target, "05_symmetry")), load_symmetry_artifacts(g2(grassmannian, "05_symmetry"));
                projection_columns=8:28, source_coordinate_indices=p -> tropical_coordinate_permutation(p, 28),
                target_coordinate_indices=coordinate_permutation)
        catch error
            info("the verification raised: " * sprint(showerror, error))
            false
        end
    end
end

"The curve of a stored complex: vertices, bounded edges and rays."
function curve_graph(complex)
    vertex = Dict{Vector{QQFieldElem},Int}()
    index(v) = get!(vertex, v, length(vertex) + 1)
    edges, rays = Tuple{Int,Int}[], Tuple{Int,Vector{QQFieldElem}}[]
    for P in maximal_polyhedra(complex)
        V = [QQ.(collect(v)) for v in vertices(P)]
        if is_bounded(P)
            push!(edges, (index(V[1]), index(V[2])))
        else
            push!(rays, (index(V[1]), QQ.(collect(first(Oscar.rays(P))))))
        end
    end
    points = Vector{Vector{QQFieldElem}}(undef, length(vertex))
    foreach(((v, k),) -> points[k] = v, vertex)
    return (; points, edges, rays)
end

"""
The tropical hexagon displayed in Section 9, in the max convention: the apples a^(1), ..., a^(6) and the trees d^(12),
..., d^(61). Sections 9 and A compare the stored hexagon with them.
"""
const HEXAGON_POINTS = [[11, 12, 4, 15, 0, 14, 13], [0, 3, 3, 6, 1, 5, 4], [2, 5, 3, 4, 1, 0, 0],
                        [9, 12, 5, 15, 12, 11, 0], [0, 3, 3, 6, 3, 2, 0], [2, 5, 3, 4, 1, 0, 2]]
const HEXAGON_LINES = [[9, 9, 12, 7, 11, 10, 10, 5, 8, 0, 7, 13, 0, 12, 11, 11, 7, 10, 10, 9, 9],
                       [3, 3, 6, 1, 5, 4, 6, 9, 4, 8, 7, 7, 0, 6, 5, 5, 7, 6, 4, 3, 3],
                       [6, 7, 12, 9, 8, 4, 10, 15, 12, 11, 7, 13, 10, 9, 0, 0, 7, 10, 4, 7, 6],
                       [11, 12, 14, 11, 10, 9, 15, 5, 0, 9, 12, 18, 15, 14, 5, 5, 12, 15, 9, 12, 11],
                       [2, 3, 6, 3, 2, 0, 6, 9, 6, 5, 3, 7, 4, 3, 3, 3, 3, 6, 0, 3, 2],
                       [11, 9, 12, 7, 11, 10, 10, 15, 8, 14, 13, 13, 0, 12, 11, 11, 13, 10, 10, 9, 11]]

###############################################################################
# G2-trees: the rows of 04_trees.csv and the trees drawn in the figures
###############################################################################

"The side of a split without taxon 7, sorted."
side_without_7(A) = sort!(collect(7 in A ? setdiff(Set(1:7), A) : A))

"The internal splits of a row of a 04_trees.csv: the side without taxon 7 (sorted) => length."
function tree_splits(r)
    lengths = Dict{Vector{Int},Int}()
    for s in split(r["internal_splits"], ';'; keepempty=false)
        side, len = split(s, ':')
        lengths[side_without_7(Set(c - '0' for c in first(split(side, '|'))))] = parse(Int, len)
    end
    return lengths
end

"The tree of a row of a 04_trees.csv: the pendant lengths of the taxa 1-7, and the internal splits of `tree_splits`."
tree_row(r) = (pendant=parse.(Int, split(r["pendant_lengths"], ';')), lengths=tree_splits(r))

"The tree metric d_ij, in the order of PAIRS, of the pendant lengths of 1-7 and the internal splits (side => length)."
tree_metric(pendant, lengths) =
    [pendant[i] + pendant[j] + sum((l for (S, l) in lengths if (i in S) != (j in S)); init=0) for (i, j) in PAIRS]

"The tree metric d_ij, in the order of PAIRS, of a row of a 04_trees.csv."
tree_metric(r::AbstractDict) = tree_metric(tree_row(r)...)

"""
The G2-trees drawn in a TikZ text, in their order, as name => block: a tree starts at its vertex v1 (v1-v7 are the taxa)
and its block runs to the start of the next. The name is the first capture of the last match of `label` before the tree
(nothing if there is none, or no `label`).
"""
function tree_blocks(tex, label=nothing)
    marks = label === nothing ? Tuple{Int,String}[] :
            [(m.offset, String(something(m.captures...))) for m in eachmatch(label, tex)]
    starts = [m.offset for m in eachmatch(r"\\coordinate \(v1\) at", tex)]
    return map(eachindex(starts)) do n
        mark = findlast(m -> first(m) < starts[n], marks)
        (mark === nothing ? nothing : marks[mark][2]) =>
            SubString(tex, starts[n], n < length(starts) ? prevind(tex, starts[n+1]) : lastindex(tex))
    end
end

"""
The G2-tree drawn in a block of `tree_blocks`, or nothing if it draws no edge: its internal splits (sides without 7,
sorted) and, when weights are drawn, the pendant lengths of 1-7 and the internal splits => length (else both nothing).
Each weight goes to its nearest edge, since some figures draw it beside its edge, or with `nearest=false` to the edges it
lies on, away from their ends. Every edge must get exactly one weight.
"""
function drawn_tree(block; nearest=true)
    point = Dict(parse(Int, m[1]) => parse.(Float64, (m[2], m[3]))
                 for m in eachmatch(r"\\coordinate \(v(\d+)\) at \(([-0-9.]+),([-0-9.]+)\)", block))
    edges = [(String(m[1]), parse(Int, m[2]), parse(Int, m[3]))
             for m in eachmatch(r"\\draw\[g2tree (internal|pendant)[^\]]*\] \(v(\d+)\)--\(v(\d+)\)", block)]
    isempty(edges) && return nothing
    # The position t of the foot of p along the edge a-b (0 at a, 1 at b), and the distance from p to the segment.
    function foot(p, (_, a, b))
        (x, y), (u, v) = point[a], point[b]
        t = ((p[1] - x) * (u - x) + (p[2] - y) * (v - y)) / ((u - x)^2 + (v - y)^2)
        s = clamp(t, 0.0, 1.0)
        return t, hypot(p[1] - x - s * (u - x), p[2] - y - s * (v - y))
    end
    on = [Int[] for _ in edges]
    for m in eachmatch(r"\\node\[g2tree weight[^\]]*\] at \(([-0-9.]+),([-0-9.]+)\) \{(\d+)\}", block)
        p, w = parse.(Float64, (m[1], m[2])), parse(Int, m[3])
        if nearest
            push!(on[argmin(k -> last(foot(p, edges[k])), eachindex(edges))], w)
        else
            for k in eachindex(edges)
                t, distance = foot(p, edges[k])
                0.1 < t < 0.9 && distance < 1e-3 && push!(on[k], w)
            end
        end
    end
    # The split of the edge a-b: the taxa reached from a without crossing it, or the others if 7 is among them.
    function split_of(a, b)
        seen, todo = Set([a]), [a]
        while !isempty(todo)
            v = pop!(todo)
            for (_, x, y) in edges, (s, w) in ((x, y), (y, x))
                s == v && !(w in seen) && Set([s, w]) != Set([a, b]) && (push!(seen, w); push!(todo, w))
            end
        end
        return side_without_7(Set(v for v in seen if v <= 7))
    end
    internal = [(k, a, b) for (k, (kind, a, b)) in enumerate(edges) if kind == "internal"]
    splits = sort([split_of(a, b) for (_, a, b) in internal])
    # One length per split: an internal edge drawn in two segments would lose one.
    allunique(splits) || error("a drawn tree with an internal edge drawn in two segments")
    all(isempty, on) && return (; splits, pendant=nothing, lengths=nothing)
    all(ws -> length(ws) == 1, on) || error("a drawn tree with an edge without exactly one weight")
    pendant = [only(on[only(k for (k, (kind, a, b)) in enumerate(edges) if kind == "pendant" && i in (a, b))]) for i in 1:7]
    return (; splits, pendant, lengths=Dict(split_of(a, b) => only(on[k]) for (k, a, b) in internal))
end

###############################################################################
# The trilinear form Φ, the quadric q and the linear forms
###############################################################################

const PHI = Dict((1, 2, 3) => 1, (4, 5, 6) => 1, (1, 4, 7) => 1, (2, 5, 7) => 1, (3, 6, 7) => 1)

"Φ(e_i, e_j, e_k) for the alternating form Φ = e123 + e456 + e147 + e257 + e367."
function phi(i, j, k)
    (i == j || j == k || i == k) && return 0
    v = [i, j, k]
    order = sortperm(v)
    inversions = count(a < b && order[a] > order[b] for a in 1:3 for b in 1:3)
    return get(PHI, Tuple(sort(v)), 0) * (isodd(inversions) ? -1 : 1)
end

"The dimension of the Lie algebra of the stabilizer of Φ in gl_7 (and, with `trace_zero`, in sl_7) over the field K."
function stabilizer_dimension(K; trace_zero=false)
    # X e_j = sum_m X[m, j] e_m acts on Φ as a derivation; the coefficient of e_ijk of X.Φ is
    # sum_m X[m,i] Φ(m,j,k) + X[m,j] Φ(i,m,k) + X[m,k] Φ(i,j,m).
    unknown(m, i) = 7 * (i - 1) + m
    equations = Vector{Vector{elem_type(K)}}()
    for (i, j, k) in ((i, j, k) for i in 1:7 for j in i+1:7 for k in j+1:7)
        row = [zero(K) for _ in 1:49]
        for m in 1:7
            row[unknown(m, i)] += phi(m, j, k)
            row[unknown(m, j)] += phi(i, m, k)
            row[unknown(m, k)] += phi(i, j, m)
        end
        push!(equations, row)
    end
    trace_zero && push!(equations, [K(i == m ? 1 : 0) for i in 1:7 for m in 1:7])
    A = matrix(K, reduce(vcat, permutedims.(equations)))
    return 49 - rank(A), A
end

"The polar form of q = x1x4 + x2x5 + x3x6 - x7^2 as a symmetric matrix."
const POLAR = [(i <= 3 && j == i + 3) || (j <= 3 && i == j + 3) ? 1 : (i == j == 7 ? -2 : 0) for i in 1:7, j in 1:7]

"The seven linear forms L_i(p) = sum_{j<k} Φ(e_i, e_j, e_k) p_jk, the coefficients of Φ(x, ., .) as a linear form in x."
phi_linear_forms() = [[QQ(phi(i, j, k)) for (j, k) in PAIRS] for i in 1:7]

###############################################################################
# The paper, section by section
###############################################################################

# The sections live in paper_results/section_N.jl, one function section_N() each; only the requested ones are
# loaded, so that a section can be run on its own.
const SECTION_TITLES = [("1", "Introduction"), ("2", "Trees and Dressians"), ("3", "Grassmannian and Flag Variety"),
                        ("4", "Characteristic Two"), ("5", "The G2-forest"), ("6", "Toric Degenerations"),
                        ("7", "The G2-Orchard"), ("8", "Tropicalizing Homogeneous Spaces"), ("9", "Hexagons"),
                        ("10", "Curves of Genus Ten"), ("A", "Data")]

println("The computational statements of \"Tropical Geometry of G2\", from the stored data", QUICK ? " (--quick)" : "")
total = @elapsed for (number, title) in SECTION_TITLES
    number in SECTIONS || continue
    run_section(number, title) do
        include(joinpath(@__DIR__, "paper_results", "section_$number.jl"))
        Base.invokelatest(Base.invokelatest(getfield, Main, Symbol("section_$number")))
    end
end
println("\n", "="^118)
@printf("%d ok, %d failed, %d warnings, %d stored, %d cited, %d skipped (%.0f s)\n", COUNTS["ok"], COUNTS["FAIL"],
        COUNTS["warn"], COUNTS["stored"], COUNTS["cited"], COUNTS["skip"], total)
isempty(FAILED) || (println("Failed:"); foreach(f -> println("  ", f), FAILED))
isempty(WARNED) || (println("Warnings (left as they are by the authors):"); foreach(f -> println("  ", f), WARNED))
exit(isempty(FAILED) ? 0 : 1)
