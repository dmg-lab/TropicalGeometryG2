#=
paper_results.jl

The computational statements of the paper "Tropical Geometry of G2" (Joswig, Kuehn, Sturmfels, Vallée), in the order
in which they appear, computed from the data stored in julia/data/ and from the certificates and reports stored next to
the scripts. For each statement the script prints what it computes and compares it with the paper:

    [ ok   ]  computed here, equal to the paper
    [FAIL  ]  computed here, different from the paper (the paper's value follows)
    [warn  ]  computed here, different from the paper, left as it is by the authors (not a failure)
    [stored]  read from a stored certificate or report, written by the script named after it
    [cited ]  a statement of the paper that cites the literature
    [skip  ]  not run (--quick, or an input that only a longer script writes)

Tables that the paper describes in words only (for instance the toric degenerations of the flag variety) are printed
as well. A summary ends the output; the exit status is 1 when a statement fails.

The helpers are in this file; the statements of each section are in julia/paper_results/section_N.jl (a function
section_N()), loaded only when the section is run.

Usage, from the repository root:

    julia --project=julia julia/paper_results.jl [--sections=1,2,...,10,A] [--quick]

--quick skips the slower checks: in_{w0}(J) recomputed (Section 3); the ±1 coefficients of the stored Gröbner bases
and the quadrilaterals that are two Gröbner cones (Section 4); in_C(G) at a second point of every prime cone and the
unimodular equivalence of the Newton-Okounkov polytopes over QQ and GF(2) (Section 6); the 23 fiber complexes of
Table tab:local_flag_complex (Section 7); the random lifts over QQ(t) and the completion of the hexagon (Section 9);
the Betti table (Section 10). The full run takes about 7 minutes and 4.5 GB, --quick about 5 minutes and 3.6 GB.
=#

using TropicalCensus
using Oscar
using Printf
using Random
import Combinatorics
using TOML
import TropicalCensus: csv_fields, cell_faces, fan_cells, load_symmetry_artifacts

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
studies(parts...) = repository("julia", "studies", parts...)

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

"The rows of a CSV file as dictionaries; fields in double quotes may contain commas, and lines from # on are skipped."
function rows(path)
    lines = filter(line -> !isempty(line) && !startswith(line, '#'), readlines(path))
    header = csv_fields(first(lines))
    return map(lines[2:end]) do line
        fields = csv_fields(line)
        length(fields) == length(header) ||
            error("$path: a line with $(length(fields)) fields, the header has $(length(header)): $line")
        Dict(zip(header, fields))
    end
end

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
function face_orbits(target)
    directory = TARGET_DIRECTORIES[target]
    F = fan(target)
    symmetry = load_symmetry_artifacts(joinpath(directory, "05_symmetry"))
    rays, lineality = fan_rays(F)
    cells = fan_cells(F)
    images(face) = unique(sort!([action[r] for r in face]) for action in symmetry.ray_actions)
    orbits = Dict{Vector{Int},Tuple{Int,Int}}()
    for orbit in symmetry.cell_orbits
        faces, _ = cell_faces(cells[first(orbit)], rays, matrix(QQ, reduce(vcat, permutedims.(lineality))))
        for (level, list) in enumerate(faces), face in list
            face_images = images(face)
            get!(orbits, minimum(face_images), (level - 1, length(face_images)))
        end
    end
    return orbits
end

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

"The coefficient vector of a linear form in the variables of its ring."
linear_vector(f, n) = [coeff(f, gen(parent(f), k)) for k in 1:n]

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

"The Pfaffian of a skew-symmetric matrix of even size, by expansion along the first row."
function pfaffian_of(M)
    n = size(M, 1)
    n == 2 && return M[1, 2]
    total = zero(M[1, 1])
    for j in 2:n
        iszero(M[1, j]) && continue
        rest = [k for k in 2:n if k != j]
        total += (iseven(j) ? 1 : -1) * M[1, j] * pfaffian_of(M[rest, rest])
    end
    return total
end

"The seven linear forms L_i(p) = sum_{j<k} Φ(e_i, e_j, e_k) p_jk, the coefficients of Φ(x, ., .) as a linear form in x."
phi_linear_forms() = [[QQ(phi(i, j, k)) for (j, k) in PAIRS] for i in 1:7]

###############################################################################
# The paper, section by section
###############################################################################

# The sections live in julia/paper_results/section_N.jl, one function section_N() each; only the requested ones are
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
