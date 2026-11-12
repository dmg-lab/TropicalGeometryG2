using Oscar
using TropicalCensus: data_path, project_path

"""
An edge from `tail` to `head`, or a leg from `tail` to infinity when `head` is
0. It runs through the polyhedra `pieces` of the intersection, in order from
`tail`, and `cells[k]` is the maximal cell of the coarsest tropical
G2-Grassmannian containing `pieces[k]`.
"""
const DecoratedEdge = @NamedTuple{tail::Int, head::Int, cells::Vector{Int}, pieces::Vector{Int}}

"The vertices, in the coordinates of `data/g2`, and the decorated edges and legs; loops and multiple edges are allowed."
const DecoratedGraph = @NamedTuple{points::Vector{Vector{QQFieldElem}}, edges::Vector{DecoratedEdge}}

is_leg(e::DecoratedEdge) = iszero(e.head)
reversed(e::DecoratedEdge) = DecoratedEdge((e.head, e.tail, reverse(e.cells), reverse(e.pieces)))

const PLUCKER_PAIRS = [(i, j) for i in 1:7 for j in i+1:7]

"The permutation of the coordinates `p_ij` of `QQ^21` sending `p_ij` to `p_π(i)π(j)`."
function plucker_relabelling(π)
    σ = [findfirst(==(minmax(π[i], π[j])), PLUCKER_PAIRS) for (i, j) in PLUCKER_PAIRS]
    return x -> (y = similar(x); y[σ] = x; y)
end

"""
The labels of `data/g2` of the labels 1..7 of the archives. `curves_g10/test0`
was computed with `Mukai_Grassmannian_g_10_split.mrdi` of commit e0154e1; the
runners since commit 40196cb load `data/g2`. The labels (2,1,3,7,5,4,6) also send the split
complex onto `data/g2`; they differ by the symmetry (1 2)(4 5), which changes
the cells but not their orbits.
"""
const RELABELLINGS = ["data/g2" => [1, 2, 3, 4, 5, 6, 7], "split" => [1, 2, 3, 7, 4, 5, 6]]

"The indices of the maximal cones `cells` containing the polyhedron `P`, once its coordinates are relabelled."
containing_cells(P::Polyhedron, cells::Vector{<:Cone}; relabel=identity) =
    findall(C -> all(v -> relabel(collect(v)) in C, vertices(P)) && all(r -> relabel(collect(r)) in C, rays(P)), cells)

"""
The graph of the segments and rays `list`, relabelled to the coordinates of
`cells`: segments are edges and rays legs, each decorated by the maximal cone
containing it. The vertices are numbered as in `build_intersection_graph`, the
endpoints of the segments first. Empty polyhedra are skipped: the solver
`certified` keeps the pairs it cannot certify, and their exact intersections
may be empty.
"""
function decorated_graph(list::Vector{Polyhedron{QQFieldElem}}, cells::Vector{<:Cone}; relabel=identity)
    points, index = Vector{QQFieldElem}[], Dict{Vector{QQFieldElem},Int}()
    vertex!(v) = get!(() -> (push!(points, v); length(points)), index, v)
    edges = DecoratedEdge[]
    pieces = findall(is_feasible, list)
    for k in [filter(k -> is_bounded(list[k]), pieces); filter(k -> !is_bounded(list[k]), pieces)]
        P = list[k]
        @req dim(P) == 1 && lineality_dim(P) == 0 "polyhedron $k is neither a segment nor a ray"
        containing = containing_cells(P, cells; relabel)
        @req length(containing) == 1 "polyhedron $k lies in the maximal cells $containing"
        ends = [vertex!(relabel(collect(v))) for v in vertices(P)]
        push!(edges, DecoratedEdge((ends[1], is_bounded(P) ? ends[2] : 0, containing, [k])))
    end
    return DecoratedGraph((points, edges))
end

"""
Replace every vertex of degree 2 and its two edges by one edge between its
neighbours, or by one leg when one of the two is a leg. A vertex is kept when
its two edges are one loop or two legs, and with `within_cells` also when they
do not lie in one maximal cell.
"""
function desubdivide(G::DecoratedGraph; within_cells::Bool=false)
    edges = Union{DecoratedEdge,Nothing}[G.edges...]
    incident = [Int[] for _ in G.points]
    for (k, e) in enumerate(G.edges)
        push!(incident[e.tail], k)
        is_leg(e) || push!(incident[e.head], k)
    end
    removed = falses(length(G.points))
    for v in eachindex(G.points)
        length(incident[v]) == 2 || continue
        i, j = incident[v]
        (i == j || is_leg(edges[i]) && is_leg(edges[j])) && continue
        within_cells && !allequal([edges[i].cells; edges[j].cells]) && continue
        is_leg(edges[i]) && ((i, j) = (j, i))
        a = edges[i].head == v ? edges[i] : reversed(edges[i])
        b = edges[j].tail == v ? edges[j] : reversed(edges[j])
        push!(edges, DecoratedEdge((a.tail, b.head, [a.cells; b.cells], [a.pieces; b.pieces])))
        edges[i] = edges[j] = nothing
        for (x, old) in ((a.tail, i), (b.head, j))
            iszero(x) || replace!(incident[x], old => length(edges))
        end
        empty!(incident[v])
        removed[v] = true
    end
    kept = findall(!, removed)
    number = zeros(Int, length(G.points))
    number[kept] = eachindex(kept)
    renumbered(e) = DecoratedEdge((number[e.tail], is_leg(e) ? 0 : number[e.head], e.cells, e.pieces))
    return DecoratedGraph((G.points[kept], [renumbered(e) for e in edges if !isnothing(e)]))
end

"The core of `G`, which has the same genus: the legs are dropped, then the leaves removed until none is left."
function core(G::DecoratedGraph)
    edges = filter(!is_leg, G.edges)
    while true
        valence = zeros(Int, length(G.points))
        foreach(e -> (valence[e.tail] += 1; valence[e.head] += 1), edges)
        any(==(1), valence) || break
        filter!(e -> valence[e.tail] > 1 && valence[e.head] > 1, edges)
    end
    kept = sort!(unique!([v for e in edges for v in (e.tail, e.head)]))
    number = zeros(Int, length(G.points))
    number[kept] = eachindex(kept)
    return DecoratedGraph((G.points[kept], [DecoratedEdge((number[e.tail], number[e.head], e.cells, e.pieces))
                                            for e in edges]))
end

"`E - V + C`, the legs not counted."
function genus(G::DecoratedGraph)
    parent = collect(eachindex(G.points))
    root(x) = parent[x] == x ? x : (parent[x] = root(parent[x]))
    edges = filter(!is_leg, G.edges)
    foreach(e -> parent[root(e.tail)] = root(e.head), edges)
    return length(edges) - length(G.points) + count(x -> root(x) == x, eachindex(G.points))
end

"The census orbit of every maximal cell of the G2-Grassmannian, and its number in the tables of maximal trees of `main.tex`."
function cell_labels()
    orbit = Dict{Int,Int}()
    for line in readlines(data_path("g2", "grassmannian", "05_symmetry_cell_orbits.csv"))[2:end]
        fields = split(line, ',')
        foreach(c -> orbit[parse(Int, c)] = parse(Int, fields[1]), split(strip(fields[4], '"'), ';'))
    end
    number = Dict{Int,Int}()
    for line in readlines(data_path("g2", "maximal_cell_tree_metrics.csv"))[2:end]
        fields = split(line, ',')
        isempty(fields[3]) || (number[parse(Int, fields[3])] = parse(Int, fields[1]))
    end
    cells = 1:maximum(keys(orbit))
    return [orbit[c] for c in cells], [number[orbit[c]] for c in cells]
end

"One row per edge or leg (`head` 0): its cells, their orbits and tree numbers, and its pieces, from `tail`."
function write_decorated_graph(path, G::DecoratedGraph, orbit, number)
    open(path, "w") do io
        println(io, "edge,tail,head,cells,orbits,tree_numbers,pieces")
        for (k, e) in enumerate(G.edges)
            println(io, "$k,$(e.tail),$(e.head),", join(e.cells, ';'), ',', join(orbit[e.cells], ';'), ',',
                    join(number[e.cells], ';'), ',', join(e.pieces, ';'))
        end
    end
end

"""
Save the NamedTuple `graphs` of `DecoratedGraph`s with Oscar, their edges as
columns: Oscar cannot save an empty vector of edges.
"""
function save_decorated_graphs(path, graphs)
    columns(G) = (points=G.points, tail=[e.tail for e in G.edges], head=[e.head for e in G.edges],
                  cells=[e.cells for e in G.edges], pieces=[e.pieces for e in G.edges])
    save(path, map(columns, graphs))
end

"The graphs saved by `save_decorated_graphs`, with their types: empty vectors come back untyped."
load_decorated_graphs(path) =
    map(load(path)) do G
        DecoratedGraph((Vector{Vector{QQFieldElem}}(G.points),
                        DecoratedEdge[DecoratedEdge(e) for e in zip(G.tail, G.head, G.cells, G.pieces)]))
    end

"""
Save the polygons `x18` of k3_surface_polygons.jl with Oscar, flat: one
vector per field of the polygons, which Oscar saves, unlike their pairs.
"""
function save_k3_polygons(path, x18)
    P, Key = x18.polygons, Tuple{Symbol,Int,Vector{Int}}
    save(path, (cells=[c for ((c, _), _) in P], subsets=[S for ((_, S), _) in P], bounded=[p.bounded for (_, p) in P],
                verts=[p.verts for (_, p) in P], keys=[Key[k for k in p.keys] for (_, p) in P],
                curve_keys=Dict{Int,Key}(x18.curve_keys), cells_of_wall=x18.cells_of_wall, t=x18.t))
end

"The polygons saved by `save_k3_polygons`, as k3_surface_polygons.jl returns them."
function load_k3_polygons(path)
    x = load(path)
    polygons = [(c, S) => (; keys, bounded, verts)
                for (c, S, bounded, verts, keys) in zip(x.cells, x.subsets, x.bounded, x.verts, x.keys)]
    return (; polygons, curve_keys=Dict{Int,Any}(x.curve_keys), cells_of_wall=x.cells_of_wall, t=x.t)
end

function describe(io::IO, name, G::DecoratedGraph)
    valence = zeros(Int, length(G.points))
    foreach(e -> (valence[e.tail] += 1; is_leg(e) || (valence[e.head] += 1)), G.edges)
    distinct_cells = [length(unique(e.cells)) for e in G.edges]
    println(io, "$name: $(length(G.points)) vertices, $(count(!is_leg, G.edges)) edges, ",
            "$(count(is_leg, G.edges)) legs, genus $(genus(G))")
    println(io, "  valences: ", join(("$d^$(count(==(d), valence))" for d in sort!(unique(valence))), ", "))
    println(io, "  edges and legs through k distinct cells: ",
            join(("k=$k: $(count(==(k), distinct_cells))" for k in sort!(unique(distinct_cells))), ", "))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (2, 3) || error("usage: decorated_curve_graph.jl ARCHIVE_DIR EXPERIMENT [OUTPUT_DIR]")
    archive, experiment = abspath(ARGS[1]), parse(Int, ARGS[2])
    output = length(ARGS) == 3 ? abspath(ARGS[3]) : project_path("results", "curves", "genus_10", basename(archive))
    mkpath(output)

    cells = collect(maximal_cones(load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))))
    list = load(joinpath(archive, "list_polyhedra_$(experiment).mrdi"))
    piece = list[findfirst(is_feasible, list)]
    found = findfirst(((_, π),) -> !isempty(containing_cells(piece, cells; relabel=plucker_relabelling(π))),
                      RELABELLINGS)
    @req !isnothing(found) "the first nonempty polyhedron lies in no maximal cell, whatever the coordinates"
    println("Coordinates of the archive: $(RELABELLINGS[found])")
    G = decorated_graph(list, cells; relabel=plucker_relabelling(last(RELABELLINGS[found])))
    graph_file = joinpath(archive, "graph_$(experiment).mrdi")
    if isfile(graph_file)
        stored = Set(minmax(src(e), dst(e)) for e in edges(load(graph_file)))
        @req stored == Set(minmax(e.tail, e.head) for e in G.edges if !is_leg(e)) "the edges differ from $graph_file"
    end

    orbit, number = cell_labels()
    graphs = (subdivided=G, desubdivided=desubdivide(G), reduced=desubdivide(G; within_cells=true),
              core=desubdivide(core(G)))
    foreach(name -> describe(stdout, name, graphs[name]), keys(graphs))
    save_decorated_graphs(joinpath(output, "decorated_graph_$(experiment).mrdi"), graphs)
    for name in (:desubdivided, :reduced, :core)
        path = joinpath(output, "decorated_graph_$(experiment)_$(name).csv")
        write_decorated_graph(path, graphs[name], orbit, number)
        println("Wrote $path")
    end
end
