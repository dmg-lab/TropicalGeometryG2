###############################################################################
# The tree metrics d_ij of the 35 orbits of maximal cells of the tropical
# G2-Grassmannian in characteristic 0 and in characteristic 2: one tree per
# orbit, that of the two tables of maximal trees of main.tex, numbered as there.
#
# Run from julia/ with `julia --project=.`; the script writes
# data/g2/maximal_cell_tree_metrics.csv, or --output=FILE, in about half
# a minute.
###############################################################################

# The numbers of the orbits and the trees of their support cells.
include(joinpath(@__DIR__, "..", "flag", "cell_orbit_trees.jl"))

const MAXIMAL_CELL_TREE_METRICS_CSV =
    TropicalGeometryG2.data_path("g2", "maximal_cell_tree_metrics.csv")
const MAXIMAL_CELL_TREE_PAIRS = [(i, j) for i in 1:6 for j in (i + 1):7]
# The six equalities of eq:dressian in main.tex.
const MAXIMAL_CELL_TREE_EQUALITIES = (
    ((1, 7), (5, 6)), ((2, 7), (4, 6)), ((3, 7), (4, 5)),
    ((4, 7), (2, 3)), ((5, 7), (1, 3)), ((6, 7), (1, 2)),
)

"""The integral tree distances d_ij, i < j, of a seven-taxon tree."""
function maximal_cell_tree_distances(tree)
    d = [tree_distance(tree, i, j) for (i, j) in MAXIMAL_CELL_TREE_PAIRS]
    all(value -> denominator(value) == 1, d) || error("a tree metric is not integral")
    return Int.(numerator.(d))
end

"""Check that `d` satisfies eq:dressian of main.tex."""
function check_g2_tree_metric(d, number)
    D(i, j) = d[findfirst(==((min(i, j), max(i, j))), MAXIMAL_CELL_TREE_PAIRS)]
    for (left, right) in MAXIMAL_CELL_TREE_EQUALITIES
        D(left...) == D(right...) || error("orbit $number has d$left != d$right")
    end
    m = (D(1, 4), D(2, 5), D(3, 6))
    count(==(maximum(m)), m) >= 2 || error("orbit $number attains max(d14, d25, d36) once")
end

"""
The 35 orbits of maximal cells of both fields, by the numbers of the tables of
maximal trees: orbit size, census orbit of each field, support cell, and the
d_ij of its tree.
"""
function maximal_cell_tree_metric_rows()
    rows = Dict{Int,Any}()
    for toric in (true, false), row in toric_tree_table_rows(; toric=toric)
        number = toric_table_label_number(row.label)
        d = maximal_cell_tree_distances(row.tree)
        if haskey(rows, number)
            old = rows[number]
            old.d == d && old.support == row.support || error("orbit $number has two trees")
            rows[number] = merge(old, (orbits=merge(old.orbits, Dict(row.orbits)),))
        else
            rows[number] = (number=number, size=row.orbit_size,
                            orbits=Dict(row.orbits), support=row.support, d=d)
        end
    end
    sort(collect(keys(rows))) == 1:length(rows) || error("the orbits are not numbered 1, 2, ...")
    for row in values(rows)
        check_g2_tree_metric(row.d, row.number)
    end
    for field in TORIC_TABLE_FIELDS
        cells = sum(row.size for row in values(rows) if haskey(row.orbits, field))
        cells == G2_TREE_FAN_SIZES[field].cells ||
            error("the orbit sizes of $field do not sum to its number of maximal cells")
    end
    return [rows[number] for number in 1:length(rows)]
end

function write_maximal_cell_tree_metrics(path=MAXIMAL_CELL_TREE_METRICS_CSV)
    rows = maximal_cell_tree_metric_rows()
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "number,orbit_size,qq_orbit,gf2_orbit,support,",
                join(("d$i$j" for (i, j) in MAXIMAL_CELL_TREE_PAIRS), ","))
        for row in rows
            println(io, join((row.number, row.size,
                              get(row.orbits, :qq, ""), get(row.orbits, :gf2, ""),
                              join(row.support, ';'), row.d...), ","))
        end
    end
    println("Wrote ", path)
    return nothing
end

function export_maximal_cell_tree_metrics_main(args=ARGS)
    outputs = [split(arg, '='; limit=2)[2] for arg in args if startswith(arg, "--output=")]
    length(outputs) == length(args) <= 1 ||
        throw(ArgumentError("the only option is --output=FILE"))
    write_maximal_cell_tree_metrics(isempty(outputs) ? MAXIMAL_CELL_TREE_METRICS_CSV : only(outputs))
end

if abspath(PROGRAM_FILE) == @__FILE__
    export_maximal_cell_tree_metrics_main()
end
