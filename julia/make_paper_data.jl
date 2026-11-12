#=
make_paper_data.jl

Writes the data of the paper "Tropical Geometry of G2" (Joswig, Kuehn, Sturmfels, Vallée) into julia/data/: the files
of the folder data/ of the Zenodo tarball, with the same names, from the ideals, the shift vector q of the first curve
of genus ten, the run of the second, and the inputs kept in git: among them the tables that give the census orbits
the names of the paper (julia/studies/g2/atlas/, julia/data/g2/maximal_cell_tree_metrics.csv). It runs the census,
in order, one Julia process per step, so that the memory of a step is returned before the next one. A step whose
outputs are stored is skipped, so the script continues where it stopped.

Usage, from the repository root:

    julia --project=julia julia/make_paper_data.jl              # every step whose outputs are missing
    julia --project=julia julia/make_paper_data.jl --list       # the steps, their time, and what is stored
    julia --project=julia julia/make_paper_data.jl --only=qq_base,qq_dressian
    julia --project=julia julia/make_paper_data.jl --from=gf2_base

    --force       run the selected steps even if their outputs are stored (the census gets --force)
    --dry-run     print the commands only
    --threads=N   Julia threads of the step zenodo, the only threaded one (default: min(8, number of CPU threads));
                  the other steps run on one thread, since polymake must not be called from two threads at once

A full run takes about 1.5 hours on a laptop with 16 GB (Intel i7-12700H); no step needs more than 6 GB. It needs
gfan 0.8 (julia/scripts/install_gfan.sh). Besides julia/data/, the steps write checkpoints to julia/work/. Then
julia/paper_results.jl checks every statement of the paper, and production/make_zenodo_tarball.sh packs julia/data/
into the Zenodo tarball. The label tables are computed by the scripts of julia/studies/g2/atlas/ (see the README
there), and the checks of the census over GF(2) by julia/studies/g2/supplemental/verify_gf2_census.jl.
=#

const JULIA_DIRECTORY = @__DIR__
const OPTIONS = filter(startswith("--"), ARGS)
const FLAGS = ("list", "force", "dry-run", "help")
const VALUE_OPTIONS = ("only", "from", "threads")
for a in OPTIONS
    a in ("--" .* FLAGS) || any(name -> startswith(a, "--$name=") && length(a) > length(name) + 3, VALUE_OPTIONS) ||
        error("unknown option $a; the options are --list, --force, --dry-run, --only=STEPS, --from=STEP, --threads=N")
end
length(OPTIONS) == length(ARGS) || error("make_paper_data.jl takes options only, not $(setdiff(ARGS, OPTIONS))")
flag(name) = "--$name" in OPTIONS
option(name) = (values = [split(a, '='; limit=2)[2] for a in OPTIONS if startswith(a, "--$name=")];
                isempty(values) ? nothing : String(only(values)))

data(parts...) = joinpath(JULIA_DIRECTORY, "data", parts...)
g2(parts...) = data("g2", parts...)
provenance(target, stage) = g2(target, "provenance", "$stage.toml")

const CENSUS = "studies/g2/run_census.jl"
const CENSUS_GF2 = "studies/g2/run_census_gf2.jl"
const CENSUS_SLONG = "studies/g2/run_census_slong.jl"
const FLAG_STAGES = ("tropicalize", "coarsen", "normalize", "trees", "symmetry", "degenerations", "projections",
                     "summarize")
const FLAG_DRESSIAN_STAGES = ("tropicalize", "coarsen", "normalize", "trees", "symmetry", "summarize")
const GF2_BASE_STAGES = ("tropicalize", "coarsen", "normalize", "trees", "symmetry", "degenerations", "summarize")
const SHORT_LONG_STAGES = ("tropicalize", "coarsen", "normalize", "symmetry")
const SLONG_STAGES = ("tropicalize", "coarsen", "normalize", "symmetry", "degenerations", "projections", "summarize")

"""
A step: the script (relative to julia/) with its arguments, the files whose presence means it is done, and its time in
minutes on the laptop. The census runs with one Gfan worker, and the large targets one stage per process: two workers,
or all the stages of a large target in one process, need more than 6 GB, against about 3 GB for each stage alone.
"""
step(name, script, arguments, outputs, minutes, purpose; threaded=false) =
    (; name, script, arguments, outputs, minutes, purpose, threaded)

const STEPS = [
    step("qq_base", CENSUS, ["--target=base", "--from=tropicalize", "--through=summarize", "--workers=1"],
         [provenance("grassmannian", "summarize")], 5, "trop(G) over QQ (Oscar's tropical_variety) and its stages"),
    step("qq_dressian", CENSUS, ["--target=dressian", "--from=tropicalize", "--through=summarize", "--workers=1"],
         [provenance("dressian", "summarize")], 14, "the G2-Dressian (Oscar's tropical_prevariety, 13 min)"),
    [step("qq_flag_$stage", CENSUS, ["--target=flag", "--stage=$stage", "--workers=1"],
          [provenance("flag", stage)], stage == "tropicalize" ? 4 : 1, "trop(F) over QQ, stage $stage")
     for stage in FLAG_STAGES]...,
    [step("qq_flag_dressian_$stage", CENSUS, ["--target=flag_dressian", "--stage=$stage", "--workers=1"],
          [provenance("flag_dressian", stage)], 1, "the G2 flag Dressian (gfan), stage $stage")
     for stage in FLAG_DRESSIAN_STAGES]...,
    [step("gf2_base_$stage", CENSUS_GF2, ["--target=base", "--stage=$stage", "--workers=1", "--bits=0"],
          [provenance("grassmannian_gf2", stage)], stage == "tropicalize" ? 3 : 1,
          "trop(G) over GF(2), stage $stage" * (stage == "tropicalize" ? " (with its tropical basis)" : ""))
     for stage in GF2_BASE_STAGES]...,
    [step("gf2_flag_$stage", CENSUS_GF2, ["--target=flag", "--stage=$stage", "--workers=1", "--bits=0"],
          [provenance("flag_gf2", stage)], stage == "tropicalize" ? 17 : 2, "trop(F) over GF(2), stage $stage")
     for stage in FLAG_STAGES]...,
    [step("short_long_$stage", CENSUS_SLONG, ["--target=short_long", "--stage=$stage", "--workers=1"],
          [provenance("grassmannian_Slong/short_long", stage)], stage == "tropicalize" ? 4 : 1,
          "Section 8: the short and long roots, stage $stage")
     for stage in SHORT_LONG_STAGES]...,
    [step("slong_$stage", CENSUS_SLONG, ["--target=slong", "--stage=$stage", "--workers=1"],
          [provenance("grassmannian_Slong", stage)], stage == "tropicalize" ? 6 : 1,
          "Section 8: the long roots, stage $stage")
     for stage in SLONG_STAGES]...,
    step("zenodo", CENSUS, ["--target=zenodo", "--stage=zenodo"],
         [provenance(".", "zenodo"), data("mukai_curves", "genus_10_method_1_pieces.mrdi")], 5,
         "point-space ideals, orbit tables, tropicalized hexagon, curve of Method 1 from q, curve of Method 2";
         threaded=true),
]

const NAMES = [s.name for s in STEPS]
done(s) = all(isfile, s.outputs)

function selected_steps()
    only_names = option("only")
    from = option("from")
    isnothing(only_names) || isnothing(from) || error("pass --only or --from, not both")
    if !isnothing(only_names)
        names = String.(split(only_names, ','))
        unknown = setdiff(names, NAMES)
        isempty(unknown) || error("unknown steps $(join(unknown, ", ")); the steps are $(join(NAMES, ", "))")
        return filter(s -> s.name in names, STEPS)
    end
    isnothing(from) && return STEPS
    first_step = findfirst(==(from), NAMES)
    isnothing(first_step) && error("unknown step $from; the steps are $(join(NAMES, ", "))")
    return STEPS[first_step:end]
end

function list_steps(steps)
    println(rpad("step", 32), rpad("min", 5), rpad("stored", 8), "purpose")
    for s in steps
        println(rpad(s.name, 32), rpad(s.minutes, 5), rpad(done(s) ? "yes" : "no", 8), s.purpose)
    end
    println("about ", sum(s.minutes for s in steps if !done(s); init=0), " minutes for the missing steps")
end

"The arguments of a step, and the command that runs it."
function command(s, threads, force)
    arguments = copy(s.arguments)
    force && s.script in (CENSUS, CENSUS_GF2, CENSUS_SLONG) && push!(arguments, "--force")
    threads = s.threaded ? threads : "1"
    return arguments, `$(Base.julia_cmd()) --startup-file=no --project=$JULIA_DIRECTORY --threads=$threads $(joinpath(JULIA_DIRECTORY, s.script)) $arguments`
end

function main()
    steps = selected_steps()
    (flag("list") || flag("help")) && return list_steps(steps)
    force = flag("force")
    threads = something(option("threads"), string(min(8, Sys.CPU_THREADS)))
    pending = filter(s -> force || !done(s), steps)
    needs_gfan = any(s -> s.script in (CENSUS, CENSUS_GF2, CENSUS_SLONG) && s.name != "zenodo", pending)
    needs_gfan && !flag("dry-run") && Sys.which("gfan_tropicalprevariety") === nothing &&
        error("gfan_tropicalprevariety of gfan 0.8 is not on the PATH; the census needs it (julia/scripts/install_gfan.sh)")
    total = 0.0
    for (k, s) in enumerate(steps)
        header = "[$k/$(length(steps))] $(s.name): $(s.purpose)"
        if done(s) && !force
            println(header, ": stored")
            continue
        end
        arguments, cmd = command(s, threads, force)
        println("\n", header, "\n  julia --project=julia", s.threaded ? " --threads=$threads" : "", " julia/", s.script,
                isempty(arguments) ? "" : " " * join(arguments, ' '))
        flush(stdout)
        flag("dry-run") && continue
        foreach(mkpath ∘ dirname, s.outputs)
        seconds = @elapsed run(Cmd(cmd; dir=dirname(JULIA_DIRECTORY)))
        total += seconds
        all(isfile, s.outputs) || error("step $(s.name) finished without writing $(join(filter(!isfile, s.outputs), ", "))")
        println("  $(s.name) done in $(round(seconds / 60; digits=1)) min")
        flush(stdout)
    end
    missing_steps = [s.name for s in STEPS if !done(s)]
    println("\n", isempty(missing_steps) ? "julia/data/ is complete" : "missing: $(join(missing_steps, ", "))",
            " ($(round(total / 60; digits=1)) min in this run)")
end

main()
