###############################################################################
# sample_genus_10: fast sampling of the curves of
# hyperplane_intersection_census.jl, in search of a curve of genus 10 whose
# core is trivalent.
#
# Experiment k draws the 5 hyperplanes of experiment k of the census
# (Xoshiro(k), shifts in -2^20:2^20, supports of sizes 2:21) and makes the
# first 4 cuts only, with the floating-point filter :float of
# intersect_complexes.jl: no exact LP but for the pairs HiGHS leaves
# undecided, and nothing stored but a row of DIR/samples/task_T.csv. The
# curve is on target when it is balanced, of genus 10, and the valences of its
# core (legs dropped, leaves removed, vertices of valence 2 smoothed) are all
# 3; the rest of the curve may have vertices of any valence. A curve of genus
# 10 with a trivalent core, balanced or not, goes to DIR/candidates.csv and is
# then verified: experiment k of
# the census is run in DIR/verified/experiment_k with the solver certified,
# all 5 cuts and every object stored, and its row goes to DIR/verified.csv.
# Its target column requires the degree 18 as well.
#
#   julia --project=. -t 2 studies/g2/curves/sample_genus_10.jl \
#     --first=K --experiments=N [--stride=S] [--task=T] [--output=DIR] [--solver=float] [--max-rss=M]
#
# runs experiments K, K+S, ..., N of them (S defaults to 1), but stops after
# the first experiment that brings the peak RSS above M MiB (6000 by
# default): polymake memory is not returned, and on hydra processes of 10G
# were killed after 1 to 5 experiments, whereas no experiment has yet
# needed 9 GB in a fresh process (2026-10-27).
# studies/g2/curves/sample_genus_10.slurm runs it as a Slurm array on
# hydra. DIR defaults to results/curves/genus_10_sampling.
###############################################################################

include(joinpath(@__DIR__, "hyperplane_intersection_census.jl"))

const SAMPLE_HEADER = "experiment,status,supports,cells,genus,core_vertices,core_edges,core_valences," *
                      "multiplicities,balanced,target,verified,seconds"

"The cuts 1 to 4 of experiment `k` by `solver`, and its row of the samples."
function sample_experiment(k, grassmannian, weights; solver)
    hyperplanes = draw_hyperplanes(k, ambient_dim(grassmannian))
    row = Dict{String,Any}("experiment" => k, "supports" => join((count(!iszero, h.normal) for h in hyperplanes[1:CURVE_STAGE]), ';'))
    complex, multiplicities, cells = grassmannian, weights, Int[]
    for s in 1:CURVE_STAGE
        _, stage, _ = cut_by(complex, multiplicities, hyperplanes[s].complex, s, dim(grassmannian) - s;
                             solver, label="Experiment $k")
        complex, multiplicities = stage.complex, stage.multiplicities
        push!(cells, length(multiplicities))
        row["cells"] = join(cells, ';')
        stage.transverse || (row["status"] = "cut $s not transverse"; return row)
    end
    L = lineality_dim(complex)
    L == 1 || (row["status"] = "lineality of dimension $L"; return row)
    list = modulo_lineality(complex)
    G = decorated_graph(list, nothing)
    C = desubdivide(core(G))
    core_valences = graph_valences(C)
    g = genus(G)
    balanced = is_balanced(list, multiplicities)
    trivalent = !isempty(core_valences) && all(==(3), core_valences)
    merge!(row, Dict("status" => "curve", "genus" => g, "core_vertices" => length(C.points), "core_edges" => length(C.edges),
                     "core_valences" => histogram(core_valences), "multiplicities" => histogram(multiplicities),
                     "balanced" => balanced, "target" => g == TARGET_GENUS && trivalent && balanced,
                     "candidate" => g == TARGET_GENUS && trivalent))
    return row
end

"Run experiment `k` of the census, certified and stored in `output/verified`, and append its row to `verified.csv`."
function verify_experiment(output, k, grassmannian, weights, cells, labels)
    directory = joinpath(output, "verified", "experiment_$k")
    start = time()
    row = is_summarized(directory) ? summary_row(directory) :
          run_census_experiment(directory, k, grassmannian, weights, cells, labels; shift=1 << 20, solver=:certified,
                                coordinates=1:ambient_dim(grassmannian), support_sizes=2:ambient_dim(grassmannian))
    haskey(row, "seconds") || (row["seconds"] = round(Int, time() - start))
    write(joinpath(directory, "summary.csv"), CENSUS_HEADER * '\n' * census_line(row) * '\n')
    path = joinpath(output, "verified.csv")
    # Tasks of the array may verify at the same time: one line per write.
    open(io -> write(io, (isfile(path) && filesize(path) > 0 ? "" : CENSUS_HEADER * '\n') * census_line(row) * '\n'), path, "a")
    println("Verified experiment $k: ", census_line(row))
    return row
end

function parse_sampling_options(args)
    options = Dict("first" => "1", "experiments" => "1", "stride" => "1", "task" => "0", "solver" => "float",
                   "output" => project_path("results", "curves", "genus_10_sampling"), "max-rss" => "6000")
    for arg in args
        m = match(r"^--([a-z-]+)=(.+)$", arg)
        !isnothing(m) && haskey(options, m[1]) || error("unknown argument $arg; see the header of $(@__FILE__)")
        options[m[1]] = m[2]
    end
    return (first=parse(Int, options["first"]), experiments=parse(Int, options["experiments"]),
            stride=parse(Int, options["stride"]), task=options["task"], solver=Symbol(options["solver"]),
            output=abspath(options["output"]), max_rss=parse(Int, options["max-rss"]))
end

if abspath(PROGRAM_FILE) == @__FILE__
    options = parse_sampling_options(ARGS)
    grassmannian, weights = g2_grassmannian()
    cells = collect(maximal_cones(load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))))
    labels = cell_labels()
    mkpath(joinpath(options.output, "samples"))
    samples = joinpath(options.output, "samples", "task_$(options.task).csv")
    isfile(samples) || write(samples, SAMPLE_HEADER * '\n')
    done = Set(parse(Int, first(split(line, ','))) for line in readlines(samples)[2:end])
    for k in options.first .+ options.stride .* (0:options.experiments-1)
        k in done && continue
        start = time()
        row = try
            sample_experiment(k, grassmannian, weights; options.solver)
        catch error
            error isa InterruptException && rethrow()
            @error "Experiment $k failed" exception=(error, catch_backtrace())
            Dict{String,Any}("experiment" => k, "status" => "failed: " * replace(sprint(showerror, error), r"[,\n]" => ' ')[1:min(end, 80)])
        end
        row["seconds"] = round(Int, time() - start)
        line = join((get(row, column, "") for column in split(SAMPLE_HEADER, ',')), ',')
        if get(row, "candidate", false)
            # Recorded before the verification, which may not finish.
            open(io -> println(io, line), joinpath(options.output, "candidates.csv"), "a")
            verified = try
                verify_experiment(options.output, k, grassmannian, weights, cells, labels)
            catch error
                error isa InterruptException && rethrow()
                @error "Verification of experiment $k failed" exception=(error, catch_backtrace())
                nothing
            end
            row["verified"] = isnothing(verified) ? "failed" : get(verified, "target", "")
            line = join((get(row, column, "") for column in split(SAMPLE_HEADER, ',')), ',')
        end
        open(io -> println(io, line), samples, "a")
        println("Sample $k: ", line, " (max RSS $(Sys.maxrss() >> 20) MiB)")
        flush(stdout)
        GC.gc(true)
        if Sys.maxrss() >> 20 > options.max_rss
            println("Peak RSS above $(options.max_rss) MiB: the rest of the experiments is left out")
            break
        end
    end
end
