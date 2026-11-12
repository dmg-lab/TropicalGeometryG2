###############################################################################
# Three decorated G2-trees side by side: the star tree of the lineality space
# with an apple on one of its pendant edges and, of a triangular cell, the tree
# of a ray with an apple at one of its inner vertices, both vertices of the
# tropical G2 flag variety, and the tree of the cell with an apple on one of
# its bounded edges, a maximal cell. The apple is the free taxon 8, whose
# pendant edge has length 0. The trees are those of triangle_cell_figure.jl,
# the rays at the size of the cell, without letters; the letters of that
# figure only name the vertex, and the bounded edge of the cell by the letters
# at its two ends. The star has its taxa where the other two have them.
#
#   julia --project=julia julia/to_tex/decorated_trees_figure.jl \
#     [--orbit=N] [--pendant=K] [--vertex=A] [--edge=AB] [--scale=FACTOR] \
#     [--taxon-scale=FACTOR] [--apple-style=flat|outlined|glossy] [--standalone] \
#     [--output=FILE]
#
# The orbit numbers are those of the tables of maximal trees, from the atlas.
###############################################################################

include(joinpath(@__DIR__, "triangle_cell_figure.jl"))

# The taxon whose pendant edge of the star carries the apple: the dark square.
const G2_DECORATED_PENDANT = 1
const G2_DECORATED_VERTEX = "A"
const G2_DECORATED_EDGE = "AB"
# The space between neighbouring trees, in centimetres.
const G2_DECORATED_GAP_CM = 0.6
# The apple in millimetres at taxon scale 0.75, hanging from its stem at the
# origin: the radius and centre of its body, and the clearance it keeps from
# the taxa and the edges of its tree, and from the edges through its stem.
const G2_APPLE_RADIUS_MM = 2.3
const G2_APPLE_CENTRE_MM = 3.5
const G2_APPLE_CLEARANCE_MM = 0.6
const G2_APPLE_STEM_CLEARANCE_MM = 0.2
# The colors of the apple, as the TikZ colors g2AppleBody and so on: a green
# body with a darker outline in the outlined styles, a dark green leaf, and a
# brown stem.
const G2_APPLE_COLORS = (Body="7CB342", Dark="33691E", Leaf="1B5E20", Stem="5D4037")

"The definitions of the colors of the apple."
g2_apple_colors_tex() =
    join("\\definecolor{g2Apple$name}{HTML}{$color}\n" for (name, color) in pairs(G2_APPLE_COLORS))

"""
How an apple of `scale` hangs from `point` of the drawn `tree`: its `angle`
in degrees, straight down or the nearest direction, in steps of 5 degrees, in
which its body keeps clear of the edges of the tree and of the `points` of its
taxa, and whether its leaf is `mirrored`, to the clockwise side, where the
edges through `point` leave more room.
"""
function g2_decorated_apple(tree, coordinates, points, point, scale)
    mm = scale / 10
    edges = [(coordinates[v], coordinates[w]) for v in keys(tree.adjacency)
             for (w, _) in tree.adjacency[v] if v < w]
    through(a, b) = g2_segment_distance(point, a, b) < 1e-9
    # The directions from `point` along the edges through it.
    directions = [atand(q[2] - point[2], q[1] - point[1]) for (a, b) in edges if through(a, b)
                  for q in (a, b) if g2_simple_distance(point, q) > 1e-9]
    for turn in vcat(0, [side * t for t in 5:5:90 for side in (1, -1)])
        angle = -90 + turn
        centre = point .+ (G2_APPLE_CENTRE_MM * mm) .* (cosd(angle), sind(angle))
        clear = all(edges) do (a, b)
            clearance = through(a, b) ? G2_APPLE_STEM_CLEARANCE_MM : G2_APPLE_CLEARANCE_MM
            g2_segment_distance(centre, a, b) >= (G2_APPLE_RADIUS_MM + clearance) * mm
        end
        clear && all(p -> g2_simple_distance(centre, p) >= (G2_APPLE_RADIUS_MM + G2_APPLE_CLEARANCE_MM) * mm,
                     points) || continue
        counterclockwise = minimum(mod(d - angle, 360) for d in directions)
        return (; angle, mirrored=counterclockwise < 360 - maximum(mod(d - angle, 360) for d in directions))
    end
    error("no room for an apple at $point")
end

# The styles of the apple: glossy, its body and leaf outlined and a white
# highlight on its body; outlined, without the highlight; flat, without the
# outlines either, as the taxa. The figures of main.tex use the flat one.
const G2_APPLE_STYLES = (:glossy, :outlined, :flat)
const G2_APPLE_STYLE = :flat

"""
    g2_apple_pic(style=G2_APPLE_STYLE)

The apple as a TikZ pic, hanging from its stem at the origin along -y, in
millimetres: a green body, a brown stem and a leaf, drawn in one of the
`G2_APPLE_STYLES`.
"""
function g2_apple_pic(style::Symbol=G2_APPLE_STYLE)
    style in G2_APPLE_STYLES ||
        throw(ArgumentError("the apple style is one of $(join(G2_APPLE_STYLES, ", "))"))
    outline(draw) = style == :flat ? "" : draw
    lines = [
        raw"pics/g2apple/.style={code={",
        raw"  \path[fill=g2AppleBody" * outline(",draw=g2AppleDark,line width=.35pt") * "]",
        raw"    (0mm,-1.95mm)..controls(.6mm,-1.25mm)and(2.3mm,-1.3mm)..(2.3mm,-3mm)",
        raw"    ..controls(2.3mm,-4.7mm)and(1.3mm,-5.8mm)..(.6mm,-5.7mm)",
        raw"    ..controls(.3mm,-5.65mm)and(-.3mm,-5.65mm)..(-.6mm,-5.7mm)",
        raw"    ..controls(-1.3mm,-5.8mm)and(-2.3mm,-4.7mm)..(-2.3mm,-3mm)",
        raw"    ..controls(-2.3mm,-1.3mm)and(-.6mm,-1.25mm)..(0mm,-1.95mm)--cycle;",
    ]
    style == :glossy && append!(lines, [
        raw"  \fill[white,opacity=.45,rotate around={25:(-.95mm,-2.9mm)}]",
        raw"    (-.95mm,-2.9mm) ellipse[x radius=.4mm,y radius=.7mm];",
    ])
    append!(lines, [
        raw"  \draw[g2AppleStem,line width=.6pt,line cap=round]",
        raw"    (0mm,0mm)..controls(.08mm,-.8mm)..(-.05mm,-2mm);",
        raw"  \path[fill=g2AppleLeaf" * outline(",draw=g2AppleLeaf!70!black,line width=.25pt") * "]",
        raw"    (.03mm,-1.15mm)..controls(.45mm,-.5mm)and(1.25mm,-.45mm)..(1.75mm,-.85mm)",
        raw"    ..controls(1.25mm,-1.35mm)and(.55mm,-1.5mm)..(.03mm,-1.15mm)--cycle;",
        raw"}},",
    ])
    return join(lines, '\n')
end

"The apple style named by the option `--apple-style` of `args`."
g2_apple_style_option(args) = Symbol(something(g2_tree_option(args, "apple-style"), string(G2_APPLE_STYLE)))

"""
    g2_decorated_trees_tex(number=G2_TRIANGLE_ORBIT; pendant=1, vertex="A", edge="AB",
                           scale=1.0, taxon_scale=0.75, apple_style=G2_APPLE_STYLE,
                           standalone=false)

The figure of the triangle of orbit `number` that the atlas draws: on the
left the star tree of the lineality space with an apple in the middle of its
pendant edge to the taxon `pendant`, named lineality; in the middle the tree
of the ray with the inner vertex lettered `vertex`, with an apple at it,
named by its orbit in the ray catalog; on the right the tree of the cell with
an apple on the bounded edge whose ends are lettered `edge`, named by its
orbit in the tables of maximal trees. The letters are those of
`g2_triangle_letters`, and are not drawn; the trees are those of
`g2_triangle_trees` with the rays at the `scale` and `taxon_scale` of the
cell, and the star is drawn as they are. The apples are drawn in `apple_style`.
"""
function g2_decorated_trees_tex(number::Int=G2_TRIANGLE_ORBIT; pendant::Integer=G2_DECORATED_PENDANT,
                                vertex::AbstractString=G2_DECORATED_VERTEX,
                                edge::AbstractString=G2_DECORATED_EDGE, scale::Real=1.0,
                                taxon_scale::Real=0.75, apple_style::Symbol=G2_APPLE_STYLE,
                                standalone::Bool=false)
    length(edge) == 2 || throw(ArgumentError("an edge is named by the letters of its two ends, as AB"))
    1 <= pendant <= 7 || throw(ArgumentError("the pendant edge of the star is named by its taxon in 1:7"))
    field, cell = g2_triangle_orbit_support(number)
    trees = g2_triangle_trees(cell; field, scale, taxon_scale, ray_scale=1.0)

    # The inner vertex `vertex` of a ray.
    found = [(index, v) for (index, letters) in enumerate(trees.ray_letters)
             for ((v, _), text) in letters if text == vertex]
    isempty(found) && throw(ArgumentError("no ray has an inner vertex lettered $vertex"))
    index, v = first(found)
    all(==((index, v)), found) || error("the letter $vertex names two vertices")
    # The trees again, without their letters.
    bare(tree, bounded_length) =
        g2_triangle_tree(tree; order=trees.order, scale, taxon_scale, bounded_length,
                         letters=Dict{Tuple{Int,Int},String}())
    ray = trees.rays[index]
    ray_drawing = bare(trees.ray_trees[index], trees.ray_bounded_length)
    ray_point = ray_drawing.coordinates[v]
    # The bounded edge of the cell whose ends carry the letters of `edge`.
    has(end_, letter) = letter in split(get(trees.letters, end_, ""), ',')
    ends = [(a, b) for (a, b) in keys(trees.letters)
            if has((a, b), string(edge[1])) && has((b, a), string(edge[2]))]
    length(ends) == 1 || throw(ArgumentError("the tree of the cell has no single bounded edge $edge"))
    a, b = only(ends)
    middle = bare(trees.selected.tree, trees.bounded_length)
    edge_point = (middle.coordinates[a] .+ middle.coordinates[b]) ./ 2
    # The star, its pendant edge to `pendant` from the centre.
    star = g2_star_tree(zeros(Int, 7))
    star_drawing = bare(star, trees.bounded_length)
    centre = first(only(star.adjacency[pendant]))
    star_point = (star_drawing.coordinates[pendant] .+ star_drawing.coordinates[centre]) ./ 2

    apple_scale = taxon_scale / 0.75
    star_apple = g2_decorated_apple(star, star_drawing.coordinates, star_drawing.points,
                                    star_point, apple_scale)
    ray_apple = g2_decorated_apple(trees.ray_trees[index], ray_drawing.coordinates,
                                   ray_drawing.points, ray_point, apple_scale)
    edge_apple = g2_decorated_apple(trees.selected.tree, middle.coordinates, middle.points,
                                    edge_point, apple_scale)

    # The trees side by side, their centres apart by their radii and the gap:
    # the ray and the cell at -shift and shift, the star at -3shift.
    shift = max(star_drawing.radius, ray_drawing.radius, middle.radius) + G2_DECORATED_GAP_CM / 2
    apple_bottom(point, apple) = point[2] + (G2_APPLE_CENTRE_MM * sind(apple.angle) - G2_APPLE_RADIUS_MM) *
                                            apple_scale / 10
    bottom = minimum(vcat([p[2] for p in star_drawing.points], [p[2] for p in ray_drawing.points],
                          [p[2] for p in middle.points],
                          [apple_bottom(star_point, star_apple), apple_bottom(ray_point, ray_apple),
                           apple_bottom(edge_point, edge_apple)]))
    point((x, y)) = "($(g2_simple_number(x)),$(g2_simple_number(y)))"
    # The mirror acts first, on the apple hanging along -y.
    apple(place, apple) = "\\pic[rotate=$(apple.angle + 90),scale=$(g2_simple_number(apple_scale))" *
                          (apple.mirrored ? ",xscale=-1" : "") * "] at $(point(place)) {g2apple};"

    tikz = sprint() do io
        println(io, "% The star with an apple on its pendant edge to taxon $pendant; ",
                "ray $ray with an apple at its inner vertex $vertex (vertex $v of its tree); ",
                "the cell with an apple on its edge $edge (vertices $a and $b of its tree).")
        println(io, "% Letters of the inner vertices of the rays, in the order of the rays $(join(trees.rays, ",")):")
        foreach(line -> println(io, "%   ", line), trees.record)
        println(io, raw"\begingroup")
        println(io, "\\definecolor{g2RayLabel}{HTML}{$G2_SQUARE_RAY_COLOR}")
        println(io, "\\definecolor{g2CellLabel}{HTML}{$G2_SQUARE_CELL_COLOR}")
        print(io, g2_apple_colors_tex())
        println(io, raw"""\tikzset{
  g2square name/.style={text=g2RayLabel,inner sep=0pt,font=\boldmath},
  g2square cell name/.style={g2square name,text=g2CellLabel},""")
        println(io, g2_apple_pic(apple_style))
        println(io, "}")
        println(io, raw"\begin{tikzpicture}[x=1cm,y=1cm,baseline=(current bounding box.center)]")
        for (place, drawing, apple_point, decoration, name) in (
                ((-3shift, 0.0), star_drawing, star_point, star_apple,
                 "g2square name] at $(point((-3shift, bottom - G2_SQUARE_NAME_GAP_CM))) {" *
                 g2_square_table_font("lineality") * "};"),
                ((-shift, 0.0), ray_drawing, ray_point, ray_apple,
                 "g2square name] at $(point((-shift, bottom - G2_SQUARE_NAME_GAP_CM))) {" *
                 g2_square_table_font("\\($(g2_square_ray_type(ray))\\)") * "};"),
                ((shift, 0.0), middle, edge_point, edge_apple,
                 "g2square cell name] at $(point((shift, bottom - G2_SQUARE_NAME_GAP_CM))) {" *
                 g2_square_table_font("\\oldstylenums{$(trees.number)}") * "};"))
            println(io, "\\begin{scope}[shift={$(point(place))}]")
            drawing.write(io)
            println(io, apple(apple_point, decoration))
            println(io, raw"\end{scope}")
            println(io, "\\node[anchor=north,", name)
        end
        println(io, raw"\end{tikzpicture}")
        println(io, raw"\endgroup")
    end
    description = "Decorated G2-trees of orbit $number ($(trees.selected.description)): the star " *
                  "with an apple on its pendant edge to taxon $pendant, the ray $ray with an apple at " *
                  "$vertex and the cell with an apple on $edge, generated by " *
                  "julia/to_tex/decorated_trees_figure.jl"
    return g2_wrap_tree_tikz(tikz, description; standalone=standalone)
end

function g2_decorated_trees_main(args=ARGS)
    script = "julia --project=julia julia/to_tex/decorated_trees_figure.jl"
    if "--help" in args
        println("Usage: $script [--orbit=N] [--pendant=K] [--vertex=A] [--edge=AB] [--scale=FACTOR] " *
                "[--taxon-scale=FACTOR] [--apple-style=flat|outlined|glossy] [--standalone] " *
                "[--output=FILE]")
        return nothing
    end
    valued = ("orbit", "pendant", "vertex", "edge", "scale", "taxon-scale", "apple-style", "output")
    for arg in args
        arg == "--standalone" || any(name -> startswith(arg, "--$name="), valued) ||
            throw(ArgumentError("unknown option $arg; use --help"))
    end
    g2_tree_option(args, "output") == "" && throw(ArgumentError("output path must not be empty"))
    number(name, default) = parse(Float64, something(g2_tree_option(args, name), default))
    tex = g2_decorated_trees_tex(
        parse(Int, something(g2_tree_option(args, "orbit"), string(G2_TRIANGLE_ORBIT)));
        pendant=parse(Int, something(g2_tree_option(args, "pendant"), string(G2_DECORATED_PENDANT))),
        vertex=something(g2_tree_option(args, "vertex"), G2_DECORATED_VERTEX),
        edge=something(g2_tree_option(args, "edge"), G2_DECORATED_EDGE),
        scale=number("scale", "1.0"),
        taxon_scale=number("taxon-scale", "0.75"),
        apple_style=g2_apple_style_option(args),
        standalone="--standalone" in args,
    )
    return g2_tree_emit(args, tex)
end

if abspath(PROGRAM_FILE) == @__FILE__
    g2_decorated_trees_main()
end
