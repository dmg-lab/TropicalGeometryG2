# Writes hexagon_example.tex and hexagon_examples.md: tropical hexagons drawn as six decorated G2-trees in the style of
# the paper, with the cones of their lines and the positions of their points.
#
#   Examples 1-5: the input flags are two random points of the tropical flag variety, random_tropical_flag(coeff=10)
#                 after Random.seed!(101), ..., Random.seed!(105).
#   Example 6:    the input flags are the t-adic valuations of two random flags over Q(t), random_flag with
#                 random_lift_element after Random.seed!(1).
#
# All hexagons are completed by tropical_complete_hexagon and checked with is_tropical_hexagon (tests/test_tropical_hexagon.jl).
#
# Run from the julia folder of the repository:   julia --project=. studies/hexagon/figures/hexagon_example_figure.jl

include(joinpath(@__DIR__, "..", "tests", "test_tropical_hexagon.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "to_tex", "render_maximal_cell_tree.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "to_tex", "toric_degeneration_table.jl"))   # CSV reader, atlas labels

const REPO = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const OUTPUT = joinpath(@__DIR__, "hexagon_example.tex")
const DATA_OUTPUT = joinpath(@__DIR__, "hexagon_examples.md")
const NEW_SEEDS = 101:105

# ---------------------------------------------------------------------------------------------
# Trees and positions of the points (max convention: tree distances d_ij = -val(p_ij), d_i8 = -val(x_i))
# ---------------------------------------------------------------------------------------------

"Leaves reachable from `start` without passing `blocked`."
function side_leaves(T, start, blocked)
    leaves = Int[]
    stack = [(start, blocked)]
    while !isempty(stack)
        v, parent = pop!(stack)
        v <= T.ntaxa && push!(leaves, v)
        for (w, _) in T.adjacency[v]
            w == parent || push!(stack, (w, v))
        end
    end
    return sort!(leaves)
end

"Distance between two vertices of the tree."
function vertex_distance(T, a, b)
    a == b && return QQ(0)
    stack = [(a, 0, QQ(0))]
    while !isempty(stack)
        v, parent, d = pop!(stack)
        for (w, l) in T.adjacency[v]
            w == parent && continue
            w == b && return d + l
            push!(stack, (w, v, d + l))
        end
    end
    error("disconnected tree")
end

"""
The position of the point with valuations `x` on the tree `T` of a line. A point at distance `t` from `u` on the edge
`(u, v)` of length `l` has distance `d(i,u) + t` to the leaves `i` beyond `u` and `d(i,v) + l - t` to the others, and these
must equal `-x_i + c`. Returns the edges that contain the point, as named tuples; a point at a vertex lies on several.
"""
function point_positions(T, x)
    found = []
    for u in sort!(collect(keys(T.adjacency))), (v, l) in T.adjacency[u]
        u < v || continue
        A, B = side_leaves(T, u, v), side_leaves(T, v, u)
        alpha = unique([QQ(x[i]) + vertex_distance(T, i, u) for i in A])
        beta = unique([QQ(x[i]) + vertex_distance(T, i, v) + l for i in B])
        length(alpha) == 1 && length(beta) == 1 || continue
        t = (beta[1] - alpha[1]) / 2
        if u <= 7 || v <= 7
            leaf, inner = u <= 7 ? (u, v) : (v, u)
            s = u <= 7 ? l - t : t
            s >= 0 && push!(found, (kind=:leaf, leaf=leaf, inner=inner, distance=s))
        elseif 0 <= t <= l
            push!(found, (kind=:edge, u=u, v=v, from_u=t, length=l, side_u=A, side_v=B))
        end
    end
    isempty(found) && error("the point does not lie on the tree")
    return found
end

"One position: on a vertex if the point lies on several edges, else the edge."
function point_position(T, x)
    found = point_positions(T, x)
    length(found) == 1 && return only(found)
    vertices = [f.kind == :leaf ? f.inner : (f.from_u == 0 ? f.u : f.v) for f in found]
    length(unique(vertices)) == 1 || error("ambiguous position")
    return (kind=:vertex, vertex=only(unique(vertices)), edges=found)
end

# ---------------------------------------------------------------------------------------------
# Cones of the lines in the coarsest fan of the tropical G2-Grassmannian, with the labels of the paper
# ---------------------------------------------------------------------------------------------

const SOURCE = load_g2_maximal_cell_source()
const LINEALITY = balanced_lineality_basis()
ray_vector(r) = [SOURCE.balanced_rays[r, j] for j in 1:21]
const _CONES = Dict{Vector{Int},Any}()
"The cone on the given rays plus the lineality space (cached: building cones is slow)."
function cone_on(rays)
    get!(_CONES, sort(rays)) do
        R = isempty(rays) ? zeros(QQ, 0, 21) : reduce(vcat, [permutedims(ray_vector(r)) for r in rays])
        positive_hull(R, Matrix{QQFieldElem}(LINEALITY))
    end
end
csv(path) = (lines = readlines(joinpath(REPO, path)); header = split(lines[1], ',');
             [Dict(zip(header, toric_table_csv_fields(l))) for l in lines[2:end]])

const CELL_ORBIT = Dict(parse(Int, c) => parse(Int, row["orbit"])
                        for row in csv("julia/data/g2/grassmannian/06_canonical_interior_degenerations.csv")
                        for c in split(row["cells"], ';'))
const ORBIT_SIZE = Dict(parse(Int, row["orbit"]) => parse(Int, row["orbit_size"])
                        for row in csv("julia/data/g2/grassmannian/06_canonical_interior_degenerations.csv"))
const PAPER_TREE = Dict(parse(Int, row["qq_orbit"]) => toric_table_label_number(row["label"])
                        for row in csv("julia/studies/g2/atlas/flag/local_flag_cell_atlas.csv") if !isempty(row["qq_orbit"]))
const RAY_ORBIT = Dict(parse(Int, r) => parse(Int, row["orbit"])
                       for row in csv("julia/studies/g2/atlas/grassmannian/ray_label_orbits.csv") for r in split(row["rays"], ';'))
# The article labels of the ray orbits, as in build_styled_ray_orbit_tables.jl.
const RAY_LABEL = (raw"\beta_1", raw"\beta_2", raw"\alpha_1", raw"\alpha_3", raw"\alpha_2", raw"\alpha_4",
                   raw"\gamma_1", raw"\gamma_2", raw"\alpha_5", raw"\gamma_3", raw"\delta", raw"\epsilon")
ray_label(r) = RAY_LABEL[RAY_ORBIT[r]]
paper_tree(cell) = PAPER_TREE[CELL_ORBIT[cell]]

"""
The maximal cells (of the coarsest fan) containing the line `L`, the smallest face of them containing it, whether it
lies in the interior of a maximal cell, and the kind of its face.
"""
function line_cone(L)
    P = QQ.(L)
    cells = [i for (i, c) in enumerate(SOURCE.cells) if P in cone_on(collect(c))]
    isempty(cells) && error("line not in the fan")
    c = collect(SOURCE.cells[first(cells)])
    faces = vcat([Int[]], [[r] for r in c], [[a, b] for (a, b) in SOURCE.edges if a in c && b in c], [c])
    face = first(f for f in faces if P in cone_on(f))
    maximal = length(cells) == 1 && length(face) == length(c)
    kind = maximal ? :maximal : length(face) == 0 ? :lineality : length(face) == 1 ? :ray : :wall
    return (cells=cells, face=face, maximal=maximal, kind=kind)
end

# ---------------------------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------------------------

as_rational(v) = Rational{BigInt}[BigInt(numerator(c)) // BigInt(denominator(c)) for c in v]

"All data of one hexagon, from its points and lines (`raw_*` are optional unshifted values for the detail pages)."
function make_example(points, lines; title, inputs, raw_points=points, raw_lines=lines, details=false)
    trees = [tree_from_min_pluecker(QQ.(L)) for L in lines]
    positions = [[point_position(trees[k], points[k]), point_position(trees[k], points[mod1(k + 1, 6)])] for k in 1:6]
    t = @elapsed cones = [line_cone(L) for L in lines]
    println("    cones of the lines: $(round(t; digits=1)) s"); flush(stdout)
    t = @elapsed check = is_tropical_hexagon([(points[k], lines[k]) for k in 1:6])
    println("    is_tropical_hexagon: $check ($(round(t; digits=1)) s)"); flush(stdout)
    return (points=points, lines=lines, raw_points=raw_points, raw_lines=raw_lines, trees=trees, positions=positions,
            cones=cones, check=check, title=title, inputs=inputs, details=details)
end

println("computing the examples ...")
const EXAMPLES = Any[]
for (n, seed) in enumerate(NEW_SEEDS)
    Random.seed!(seed)
    f1, f2 = random_tropical_flag(coeff=10), random_tropical_flag(coeff=10)
    t = @elapsed H = tropical_complete_hexagon(f1, f2)       # vectors shifted to minimum 0
    println("  example $n: completed in $(round(t; digits=1)) s"); flush(stdout)
    push!(EXAMPLES, make_example(as_rational.(H.points), as_rational.(H.lines);
        title="Example $n",
        inputs="Input flags \\((p_1,L_1)\\), \\((p_4,L_4)\\): two random points of \\(\\trop(\\mathcal F)\\), " *
               "\\texttt{random\\_tropical\\_flag(coeff=10)} after \\texttt{Random.seed!($seed)}."))
    println("  example $n (seed $seed): is_tropical_hexagon = ", EXAMPLES[end].check)
end
let
    Random.seed!(1)
    K, _ = rational_function_field(QQ, "t")
    sample = () -> random_lift_element(K)
    L = complete_hexagon(random_flag(K; sample), random_flag(K; sample); normalize=false)
    f1 = (tropicalize(L.points[1]), tropicalize(L.lines[1]))
    f4 = (tropicalize(L.points[4]), tropicalize(L.lines[4]))
    raw = tropical_complete_hexagon(f1, f4; normalize=false)
    raw.points == tropicalize.(L.points) && shift_to_zero.(raw.lines) == shift_to_zero.(tropicalize.(L.lines)) ||
        error("example 6 does not agree with the valuations of its lift")
    H = tropical_complete_hexagon(f1, f4)
    push!(EXAMPLES, make_example(as_rational.(H.points), as_rational.(H.lines);
        title="Example 6",
        inputs="\\textbf{The input flags of this example were computed differently:} they are the \\(t\\)-adic valuations " *
               "of two random flags over \\(\\mathbb Q(t)\\), \\texttt{random\\_flag} with \\texttt{random\\_lift\\_element} " *
               "after \\texttt{Random.seed!(1)}, not random points of the fan.",
        raw_points=[Int.(p) for p in raw.points], raw_lines=[Int.(l) for l in raw.lines], details=true))
    println("  example 6: is_tropical_hexagon = ", EXAMPLES[end].check)
end

# ---------------------------------------------------------------------------------------------
# TeX helpers
# ---------------------------------------------------------------------------------------------

function tex_number(q)
    q = QQ(q)
    denominator(q) == 1 && return string(numerator(q))
    s = q < 0 ? "-" : ""
    return "$(s)\\tfrac{$(abs(numerator(q)))}{$(denominator(q))}"
end
"A rational as `a/b` (for the small tables)."
function slash_number(q)
    q = QQ(q)
    s = denominator(q) == 1 ? string(abs(numerator(q))) : "$(abs(numerator(q)))/$(denominator(q))"
    return q < 0 ? "\\(-\\)" * s : s
end
plain_number(q) = (q = QQ(q); denominator(q) == 1 ? string(numerator(q)) : "$(numerator(q))/$(denominator(q))")
taxa(set) = join(set)
split_tex(A, B) = "$(taxa(A))\\,|\\,$(taxa(B))"
const POINT_COLORS = ("2E7D32", "7B1FA2", "C2185B", "00838F", "6D4C41", "455A64")   # p1..p6
point_color(i) = "hexP$i"

"Leaf sets of the components of the tree around the inner vertex `v`."
vertex_components(T, v) = sort!([side_leaves(T, w, v) for (w, _) in T.adjacency[v]]; by=c -> (length(c), c))

function position_tex(T, pos)
    if pos.kind == :leaf
        return "leaf edge of \\(\\mathbf{$(pos.leaf)}\\), at distance \\($(tex_number(pos.distance))\\) from its inner vertex"
    elseif pos.kind == :edge
        a, b = pos.from_u, pos.length - pos.from_u
        return "edge \\($(split_tex(pos.side_u, pos.side_v))\\) of length \\($(tex_number(pos.length))\\), at distance " *
               "\\($(tex_number(a))\\) from the end at \\($(taxa(pos.side_u))\\) and \\($(tex_number(b))\\) from the end at \\($(taxa(pos.side_v))\\)"
    else
        return "inner vertex \\(" * join(taxa.(vertex_components(T, pos.vertex)), "\\,|\\,") * "\\)"
    end
end

function position_text(T, pos)
    if pos.kind == :leaf
        return "leaf edge of taxon $(pos.leaf), at distance $(plain_number(pos.distance)) from its inner vertex"
    elseif pos.kind == :edge
        return "edge $(taxa(pos.side_u))|$(taxa(pos.side_v)) of length $(plain_number(pos.length)), at distance " *
               "$(plain_number(pos.from_u)) from the end at $(taxa(pos.side_u))"
    else
        return "inner vertex " * join(taxa.(vertex_components(T, pos.vertex)), "|")
    end
end

vector_tex(v) = join(string.(v), ",\\,")

"Short description of the cone of a line: the maximal cell and the paper's tree, or what it is instead."
function cone_short(c; tex=true)
    trees = join(["#$i (tree $(paper_tree(i)))" for i in c.cells], ", ")
    tex && (trees = replace(trees, "#" => "\\#"))
    c.maximal && return "maximal cell $(trees)"
    what = c.kind == :wall ? "on a wall" : c.kind == :ray ? "on a ray" : "in the lineality space"
    return "not maximal: $(what) of the maximal cells $(trees)"
end

# ---------------------------------------------------------------------------------------------
# One page: six decorated G2-trees around the hexagon, the information and the table of the flags
# ---------------------------------------------------------------------------------------------

const TREE_SCALE = 1.0
const TREE_BOUNDED = 0.95
const CENTER_X = -3.85            # the ring of trees fills the left part of the page
const RING_X, RING_Y = 7.05, 6.0
const HEX_RADIUS = 1.8
# Hexagon vertex p_k at angle 90 - 60(k-1) (clockwise from the top); line L_k at the middle of p_k p_{k+1}.
vertex_angle(k) = deg2rad(90 - 60 * (k - 1))
side_angle(k) = deg2rad(60 - 60 * (k - 1))

lerp(a, b, f) = (a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2]))
unit_normal(a, b) = (d = (b[1] - a[1], b[2] - a[2]); n = hypot(d...); (-d[2] / n, d[1] / n))

"Drawn position of a point on the tree with drawn coordinates `C`, and the direction of its edge."
function marker_geometry(C, pos, other_on_same_leaf)
    if pos.kind == :leaf
        a, b = C[pos.inner], C[pos.leaf]
        # Several points on one leaf edge keep their order; the drawn position is only qualitative there.
        distances = sort!(unique(vcat([pos.distance], other_on_same_leaf)))
        f = length(distances) == 1 ? 0.48 : 0.30 + 0.38 * (findfirst(==(pos.distance), distances) - 1) / (length(distances) - 1)
        return lerp(a, b, f), (a, b)
    elseif pos.kind == :edge
        a, b = C[pos.u], C[pos.v]
        return lerp(a, b, 0.12 + 0.76 * Float64(pos.from_u / pos.length)), (a, b)
    else
        p = C[pos.vertex]
        e = first(pos.edges)
        q = e.kind == :leaf ? C[e.leaf] : C[e.u == pos.vertex ? e.v : e.u]
        return p, (p, q)
    end
end

function tree_body(ex, k, tag)
    T = ex.trees[k]
    tikz = g2_tree_tikz(T; palette=:color, label=true, scale=TREE_SCALE, taxon_scale=0.7,
                        bounded_length=TREE_BOUNDED, edge_labels=:none)
    body = split(split(tikz, r"\\begin\{tikzpicture\}\[[^\n]*\]\n")[2], "\\end{tikzpicture}")[1]
    body = replace(body, "(v" => "($(tag)v")
    C = g2_minimal_tree_coordinates(T; scale=TREE_SCALE, bounded_length=TREE_BOUNDED)
    io = IOBuffer()
    println(io, body)
    for u in sort!(collect(keys(T.adjacency))), (v, l) in T.adjacency[u]
        (u < v && u > 7 && v > 7) || continue
        a, b = C[u], C[v]
        m = lerp(a, b, 0.5)
        n = unit_normal(a, b)
        sgn = m[1] * n[1] + m[2] * n[2] >= 0 ? 1 : -1
        @printf(io, "\\node[hex length] at (%.4f,%.4f) {\\(%s\\)};\n", m[1] - sgn * 0.25n[1], m[2] - sgn * 0.25n[2], tex_number(l))
    end
    for (j, idx) in enumerate((k, mod1(k + 1, 6)))
        pos, other = ex.positions[k][j], ex.positions[k][3 - j]
        same_leaf = pos.kind == :leaf && other.kind == :leaf && other.leaf == pos.leaf ? [other.distance] : QQFieldElem[]
        p, (a, b) = marker_geometry(C, pos, same_leaf)
        n = unit_normal(a, b)
        sgn = p[1] * n[1] + p[2] * n[2] >= 0 ? 1 : -1
        style = pos.kind == :vertex ? "hex vertex point" : "hex point"
        @printf(io, "\\node[%s,fill=%s] at (%.4f,%.4f) {};\n", style, point_color(idx), p[1], p[2])
        @printf(io, "\\node[hex point label,text=%s] at (%.4f,%.4f) {\\(p_{%d}\\)};\n",
                point_color(idx), p[1] + sgn * 0.38n[1], p[2] + sgn * 0.38n[2], idx)
    end
    return String(take!(io))
end

"The table of the six flags: one column per flag (p_k, L_k), rows for the 7 point and 21 Plücker coordinates."
function flag_columns_table(ex)
    io = IOBuffer()
    println(io, raw"{\scriptsize\setlength{\tabcolsep}{2.6pt}\renewcommand{\arraystretch}{0.9}%")
    println(io, raw"\begin{tabular}{@{}r*{6}{r}@{}}\toprule")
    println(io, " & " * join(["\\((p_{$k},L_{$k})\\)" for k in 1:6], " & ") * " \\\\")
    println(io, " & " * join([k in (1, 4) ? "\\textit{input}" : "" for k in 1:6], " & ") * " \\\\\\midrule")
    println(io, "\\multicolumn{7}{@{}l}{\\textit{point} \\(p_k\\)} \\\\")
    for i in 1:7
        println(io, "\\(x_{$i}\\) & " * join([slash_number(ex.points[k][i]) for k in 1:6], " & ") * " \\\\")
    end
    println(io, "\\midrule\\multicolumn{7}{@{}l}{\\textit{line} \\(L_k\\), Pl\\\"ucker coordinates} \\\\")
    for (c, (i, j)) in enumerate(PAIRS)
        println(io, "\\($i$j\\) & " * join([slash_number(ex.lines[k][c]) for k in 1:6], " & ") * " \\\\")
    end
    println(io, raw"\bottomrule\end{tabular}}")
    return String(take!(io))
end

function info_block(ex)
    nonmax = [k for k in 1:6 if !ex.cones[k].maximal]
    io = IOBuffer()
    print(io, "{\\LARGE\\bfseries $(ex.title)}\\\\[5pt]\n{\\small ")
    print(io, ex.inputs, " The other four flags are computed by \\texttt{tropical\\_complete\\_hexagon}.\\\\[3pt]\n")
    status = ex.check ? "\\textcolor{hexGood}{\\textbf{true}}" : "\\textcolor{hexBad}{\\textbf{false}}"
    print(io, "\\texttt{is\\_tropical\\_hexagon}: $(status)\\\\[3pt]\n")
    if isempty(nonmax)
        print(io, "All six lines lie in the interior of maximal cells.")
    else
        print(io, "\\textcolor{hexBad}{\\textbf{Not in a maximal cell:} " * join(["\\(L_{$k}\\)" for k in nonmax], ", ") *
                  " (red captions).}")
    end
    print(io, "\\\\[3pt]\nVectors are valuations (min convention), each shifted to minimum \\(0\\). " *
              "Trees: \\(d_{ij} = -\\mathrm{val}(L_k)_{ij}\\), bold numbers are bounded edge lengths, the points are " *
              "decorations; on leaf edges only their order is drawn faithfully.}")
    return String(take!(io))
end

function example_page(ex)
    io = IOBuffer()
    println(io, raw"\begin{tikzpicture}[x=1cm,y=1cm]")
    println(io, raw"\useasboundingbox (-13.7,-9.2) rectangle (13.7,9.2);")
    for k in 1:6
        a, b = vertex_angle(k), vertex_angle(k + 1)
        @printf(io, "\\draw[g2EdgeBlue,line width=1.2pt] (%.4f,%.4f) -- node[font=\\normalsize,text=g2EdgeBlue,fill=white,inner sep=1pt] {\\(L_{%d}\\)} (%.4f,%.4f);\n",
                CENTER_X + HEX_RADIUS * cos(a), HEX_RADIUS * sin(a), k, CENTER_X + HEX_RADIUS * cos(b), HEX_RADIUS * sin(b))
    end
    for k in 1:6
        a = vertex_angle(k)
        @printf(io, "\\node[hex point,fill=%s,minimum size=9pt] at (%.4f,%.4f) {};\n", point_color(k),
                CENTER_X + HEX_RADIUS * cos(a), HEX_RADIUS * sin(a))
        @printf(io, "\\node[hex point label,text=%s] at (%.4f,%.4f) {\\(p_{%d}\\)};\n", point_color(k),
                CENTER_X + 1.36HEX_RADIUS * cos(a), 1.36HEX_RADIUS * sin(a), k)
    end
    @printf(io, "\\node[font=\\footnotesize,text=black!60,align=center] at (%.4f,0) {inputs\\\\[-1pt]\\((p_1,L_1)\\), \\((p_4,L_4)\\)};\n", CENTER_X)
    for k in 1:6
        φ = side_angle(k)
        x, y = CENTER_X + RING_X * cos(φ), RING_Y * sin(φ)
        c = ex.cones[k]
        @printf(io, "\\begin{scope}[shift={(%.4f,%.4f)}]\n", x, y)
        print(io, tree_body(ex, k, "E$(k)"))
        tag = k in (1, 4) ? " (input)" : ""
        title_style = c.maximal ? "hex title" : "hex title,text=hexBad"
        mark = c.maximal ? "" : " \\(\\boldsymbol{!}\\)"
        @printf(io, "\\node[%s] at (0,%.3f) {\\(L_{%d}\\)%s%s};\n", title_style, 2.12TREE_SCALE + 0.55, k, tag, mark)
        caption_style = c.maximal ? "hex caption" : "hex caption,text=hexBad"
        @printf(io, "\\node[%s] at (0,%.3f) {%s};\n", caption_style, -2.12TREE_SCALE - 0.55, cone_short(c))
        println(io, raw"\end{scope}")
    end
    println(io, "\\node[anchor=north west,text width=7.25cm,align=left,inner sep=0pt] at (6.35,9.1) {" * info_block(ex) * "};")
    println(io, "\\node[anchor=south east,inner sep=0pt] at (13.65,-9.15) {" * flag_columns_table(ex) * "};")
    println(io, raw"\end{tikzpicture}")
    return String(take!(io))
end

# ---------------------------------------------------------------------------------------------
# The detail pages of the earlier example (as before)
# ---------------------------------------------------------------------------------------------

function flags_table(ex)
    io = IOBuffer()
    println(io, raw"\begin{tabular}{@{}c l l@{}}\toprule")
    println(io, raw"& point \(p_k\in\RR^7\) & line \(L_k\in\RR^{21}\) (coordinates \(12,13,\dots,17,23,\dots,67\)) \\\midrule")
    for k in 1:6
        role = k in (1, 4) ? "input" : "output"
        println(io, "\\((p_{$k},L_{$k})\\) & \\(($(vector_tex(ex.raw_points[k])))\\) & \\(($(vector_tex(ex.raw_lines[k])))\\) \\\\[-1pt]")
        println(io, "\\multicolumn{1}{c}{\\scriptsize $role} &&\\\\[2pt]")
    end
    println(io, raw"\bottomrule\end{tabular}")
    return String(take!(io))
end

function cones_table(ex)
    io = IOBuffer()
    println(io, raw"\begin{tabular}{@{}c p{4.1cm} p{3.7cm} p{6.9cm} p{6.9cm}@{}}\toprule")
    println(io, raw"line & cone in \(\trop(\mathcal G)\) & bounded edges of the \(\Ggroup\)-tree & first point & second point \\\midrule")
    for k in 1:6
        c = ex.cones[k]
        rays = join(["\\($(ray_label(r))\\)\\,(\\#$r)" for r in c.face], ", ")
        cells = join(["\\#$i (tree \\(\\mathbf{$(paper_tree(i))}\\), orbit size $(ORBIT_SIZE[CELL_ORBIT[i]]))" for i in c.cells], "; ")
        kind = c.maximal ? "interior of the maximal cell $cells, spanned by $rays" :
                           "relative interior of the \\(2\\)-dimensional face spanned by $rays, a wall of the maximal cells $cells"
        splits = join(["\\($(split_tex(s.side, s.complement))\\): \\($(tex_number(s.length))\\)" for s in internal_splits(ex.trees[k])], "\\newline ")
        a, b = k, mod1(k + 1, 6)
        println(io, "\\(L_{$k}\\) & $kind & $splits & \\(p_{$a}\\): $(position_tex(ex.trees[k], ex.positions[k][1])) & " *
                    "\\(p_{$b}\\): $(position_tex(ex.trees[k], ex.positions[k][2])) \\\\\\addlinespace")
    end
    println(io, raw"\bottomrule\end{tabular}")
    return String(take!(io))
end

function detail_pages(ex)
    io = IOBuffer()
    print(io, raw"""
\newpage
\section*{Example 6 in detail}

The hexagon \(p_1 - L_1 - p_2 - L_2 - \cdots - p_6 - L_6 - p_1\) with flags \((p_k, L_k)\), \(p_k\in L_k\cap L_{k-1}\)
(\(L_0=L_6\)). \textbf{Unlike in Examples 1--5, the input flags \((p_1,L_1)\) and \((p_4,L_4)\) were computed differently:}
they are the \(t\)-adic valuations of two random flags over \(\mathbb Q(t)\),
\texttt{random\_flag} with \texttt{random\_lift\_element} after \texttt{Random.seed!(1)}, not random points of the fan.
The other four flags were computed purely tropically by \texttt{tropical\_complete\_hexagon}, and agree with the
valuations of the classical hexagon of the lift. Coordinates below are the unshifted valuations (min convention).

\medskip
{\small
""")
    print(io, flags_table(ex))
    print(io, raw"""
}

\medskip
\textbf{The pictures.} As in the paper, each line \(L_k\) is drawn as its \(\Ggroup\)-tree with tree distances
\(d_{ij} = -\mathrm{val}(L_k)_{ij}\) (max convention), and the taxa \(1,\dots,7\) are the shapes of the paper (blue
\(1,2,3\), yellow \(4,5,6\), heart \(7\)). Its two points \(p_k, p_{k+1}\) are decorations: \(d_{i8} = -\mathrm{val}(p)_i\), with the
pendant edge at the decoration normalized to \(0\). Bold numbers are the lengths of the bounded edges. Pendant
lengths are not shown: they change with the additive constant of \(L_k\), and the leaves are really rays. On a
bounded edge a point is drawn in proportion to its distances; on a leaf edge only the order of the points is
faithful. Exact positions are in the table below.

\newpage
\section*{Example 6: cones, trees and positions of the points}

The cone of each line is the smallest cone of the coarsest fan of \(\trop(\mathcal G)\) containing it (cell numbers
\#\,of the census, rays with the labels of the paper's Figure of ray trees, maximal \(\Ggroup\)-trees numbered as in the
paper's Figures of maximal trees). A bounded edge is written as the split \(A\,|\,B\) of the taxa.

\medskip
{\small
""")
    print(io, cones_table(ex))
    print(io, raw"""
}

\medskip
\textbf{Observations.}
""")
    walls = [k for k in 1:6 if !ex.cones[k].maximal]
    maxi = [k for k in 1:6 if ex.cones[k].maximal]
    println(io, "The lines \\(" * join(["L_{$k}" for k in maxi], ", ") * "\\) lie in the interior of maximal cells; the lines \\(" *
                join(["L_{$k}" for k in walls], ", ") * "\\) lie on walls, so their trees have one bounded edge less.")
    1 in walls && println(io, "In particular the input line \\(L_1\\) is not generic: the lift has put it on a wall.")
    for k in 1:6, j in 1:2
        ex.positions[k][j].kind == :vertex || continue
        idx = j == 1 ? k : mod1(k + 1, 6)
        println(io, "The point \\(p_{$idx}\\) sits exactly at an inner vertex of the tree of \\(L_{$k}\\).")
    end
    return String(take!(io))
end

# ---------------------------------------------------------------------------------------------
# The documents
# ---------------------------------------------------------------------------------------------

function styles()
    tikz = g2_tree_tikz(EXAMPLES[1].trees[1]; palette=:color, label=true, scale=TREE_SCALE, edge_labels=:none)
    block = split(tikz, "\\begin{tikzpicture}")[1]
    block = replace(block, "\\begingroup" => "")             # colours and g2tree styles, for the whole document
    io = IOBuffer()
    println(io, block)
    for i in 1:6
        println(io, "\\definecolor{$(point_color(i))}{HTML}{$(POINT_COLORS[i])}")
    end
    println(io, raw"""\definecolor{hexBad}{HTML}{C62828}
\definecolor{hexGood}{HTML}{2E7D32}
\tikzset{
  hex point/.style={circle,draw=white,line width=.9pt,inner sep=0pt,minimum size=7.4pt},
  hex vertex point/.style={circle,draw=white,line width=1.1pt,inner sep=0pt,minimum size=9.5pt},
  hex point label/.style={font=\small\bfseries,inner sep=1pt},
  hex length/.style={font=\sffamily\fontsize{8}{9}\selectfont\bfseries,text=g2EdgeBlue,inner sep=.5pt},
  hex title/.style={font=\large\bfseries,text=g2EdgeBlue},
  hex caption/.style={font=\footnotesize,text width=5.3cm,align=center,text=black!75},
}""")
    return String(take!(io))
end

function document()
    io = IOBuffer()
    print(io, raw"""\documentclass[10pt]{article}
\usepackage[a4paper,landscape,margin=1.1cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb,bm}
\usepackage{xcolor}
\usepackage{tikz}
\usetikzlibrary{calc,arrows.meta}
\usepackage{booktabs}
\usepackage{array}
\pagestyle{empty}
\newcommand{\RR}{\mathbb{R}}
\newcommand{\Ggroup}{\operatorname{G}_2}
\DeclareMathOperator{\trop}{trop}
\setlength{\parindent}{0pt}
""")
    print(io, styles())
    println(io, raw"\begin{document}")
    for (n, ex) in enumerate(EXAMPLES)
        n > 1 && println(io, raw"\newpage")
        println(io, "{\\centering")
        print(io, example_page(ex))
        println(io, "\\par}")
        ex.details && print(io, detail_pages(ex))
    end
    println(io, "\n\\end{document}")
    return String(take!(io))
end

"The readable data file: per example the six flags (input and computed), the cones of the lines and the positions."
function data_file()
    io = IOBuffer()
    println(io, "# Tropical hexagons of `hexagon_example.tex`\n")
    println(io, "Generated by [hexagon_example_figure.jl](hexagon_example_figure.jl). The hexagon is",
                " `p_1 - L_1 - p_2 - ... - p_6 - L_6 - p_1`; the flags are `(p_k, L_k)` with `p_k` on `L_k` and `L_{k-1}`.",
                " `(p_1, L_1)` and `(p_4, L_4)` are the inputs, the other four flags are computed by `tropical_complete_hexagon`.",
                " Vectors are valuations (min convention), shifted to minimum 0; points have coordinates `x_1..x_7`, lines the",
                " Plücker coordinates `12, 13, ..., 17, 23, ..., 67`.\n")
    for ex in EXAMPLES
        println(io, "## ", ex.title, "\n")
        inputs = replace(ex.inputs, r"\\textbf\{([^}]*)\}" => s"**\1**", r"\\texttt\{([^}]*)\}" => s"`\1`", "\\_" => "_",
                         "\\(" => "", "\\)" => "", "\\trop(\\mathcal F)" => "trop(F)", "\\mathbb Q" => "Q")
        println(io, inputs, "\n")
        println(io, "`is_tropical_hexagon`: **", ex.check, "**\n")
        println(io, "```")
        for k in 1:6
            role = k in (1, 4) ? "input   " : "computed"
            println(io, "($role) p_$k = [", join(plain_number.(ex.points[k]), ", "), "]")
            println(io, "           L_$k = [", join(plain_number.(ex.lines[k]), ", "), "]")
        end
        println(io, "```\n")
        println(io, "| line | cone | bounded edges | first point | second point |")
        println(io, "|---|---|---|---|---|")
        for k in 1:6
            splits = join(["$(taxa(s.side))|$(taxa(s.complement)): $(plain_number(s.length))" for s in internal_splits(ex.trees[k])], "; ")
            a, b = k, mod1(k + 1, 6)
            println(io, "| L_$k | ", cone_short(ex.cones[k]; tex=false), " | ", splits, " | p_$a: ",
                        position_text(ex.trees[k], ex.positions[k][1]), " | p_$b: ", position_text(ex.trees[k], ex.positions[k][2]), " |")
        end
        println(io)
    end
    return String(take!(io))
end

all(ex -> ex.check, EXAMPLES) || @warn "some example is not a tropical hexagon according to is_tropical_hexagon"
open(io -> print(io, document()), OUTPUT, "w")
open(io -> print(io, data_file()), DATA_OUTPUT, "w")
println("wrote ", OUTPUT, " and ", DATA_OUTPUT)
