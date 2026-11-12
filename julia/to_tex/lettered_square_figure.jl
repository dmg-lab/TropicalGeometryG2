###############################################################################
# A quadrilateral cell of the tropical G2-Grassmannian with lettered trees, as
# triangle_cell_figure.jl draws a triangle: a square whose corners are circles
# holding the G2-trees of its four rays, clockwise from the top left corner in
# their cyclic order from the smallest, and the G2-tree of the cell in its
# middle, without lengths. Every ray keeps its own tree, unlike in
# square_cell_figure.jl, and its inner vertices are lettered A, B, C, ... in
# the order of the rays (g2_triangle_letters); with --cell-letters the tree of
# the cell has them too, at the same ends of the same splits.
#
#   julia --project=julia julia/to_tex/lettered_square_figure.jl \
#     [--rays=I,J,K,L] [--scale=FACTOR] [--taxon-scale=FACTOR] \
#     [--bounded-length=LENGTH] [--ray-scale=FACTOR] [--ray-bounded-length=LENGTH] \
#     [--name-below] [--cell-letters] [--standalone] [--output=FILE]
#
# The cell is over QQ, by default that of square_cell_figure.jl.
###############################################################################

include(joinpath(@__DIR__, "triangle_cell_figure.jl"))

"""
    g2_lettered_square_tex(cell=G2_SQUARE_RAYS; scale=1.0, taxon_scale=0.75,
                           bounded_length=nothing, ray_scale=0.8, ray_bounded_length=nothing,
                           name_below=false, cell_letters=false, standalone=false)

The figure of the square `cell` of the fan over QQ, given by its index or its
four rays: the trees of `g2_triangle_trees`, with their options, at the
corners of a square as in `g2_square_cell_tex`, each ray labeled outside its
circle on the diagonal by its orbit in the ray catalog, and the cell by its
orbit in the tables of maximal trees, inside the bottom side, or below it,
between the circles, when `name_below`. The tree of the cell has the
letters of the rays only with `cell_letters`; without them the letters of a
ray stand at their vertices, by `g2_vertex_letter_places`. The square is the smallest whose
sides between the circles are `G2_SQUARE_SIDE_MIN_CM` long and that holds the
tree of the cell clear of the circles, the sides, and the label of the cell
inside. The picture ends with the hook `\\gtwocelloverlay` of
`g2_triangle_cell_tex`, the centre of the square at the origin.
"""
function g2_lettered_square_tex(cell=G2_SQUARE_RAYS; scale::Real=1.0, taxon_scale::Real=0.75,
                                bounded_length=nothing, ray_scale::Real=0.8, ray_bounded_length=nothing,
                                name_below::Bool=false, cell_letters::Bool=false, standalone::Bool=false,
                                source=nothing)
    (; rays, number, selected, order, middle, corners, record, bounded_length, ray_bounded_length) =
        g2_triangle_trees(cell; field=:qq, scale, taxon_scale, bounded_length, ray_scale,
                          ray_bounded_length, cell_letters, source)
    length(rays) == 4 || throw(ArgumentError("$(Tuple(rays)) is a triangle, not a square"))

    radius = maximum(corner -> corner.radius, corners) + G2_SQUARE_CIRCLE_PAD_CM
    name_width, name_height = G2_SQUARE_CELL_NAME_CM
    # Half the side of the square, by steps of a tenth of a millimetre.
    clear(half) = all(middle.points) do point
        max(abs(point[1]), abs(point[2])) <= half - G2_SQUARE_CENTRE_GAP_CM &&
            all(g2_simple_distance(point, (x * half, y * half)) >= radius + G2_SQUARE_CENTRE_GAP_CM
                for x in (-1, 1) for y in (-1, 1)) &&
            (name_below || abs(point[1]) > name_width / 2 + G2_SQUARE_CENTRE_GAP_CM ||
             point[2] >= -half + G2_SQUARE_NAME_GAP_CM + name_height + G2_SQUARE_CENTRE_GAP_CM)
    end
    half = ceil((radius + G2_SQUARE_SIDE_MIN_CM / 2) * 100) / 100
    while !clear(half)
        half += 0.01
    end
    places = [(-half, half), (half, half), (half, -half), (-half, -half)]
    point((x, y)) = "($(g2_simple_number(x)),$(g2_simple_number(y)))"

    tikz = sprint() do io
        println(io, "% Bounded edges of length $(g2_simple_number(bounded_length)) in the cell and " *
                    "$(g2_simple_number(ray_bounded_length)) in the rays, before scaling.")
        println(io, "% Letters of the inner vertices of the rays, in the order of the rays:")
        foreach(line -> println(io, "%   ", line), record)
        println(io, raw"\begingroup")
        println(io, "\\definecolor{g2RayLabel}{HTML}{$G2_SQUARE_RAY_COLOR}")
        println(io, "\\definecolor{g2CellLabel}{HTML}{$G2_SQUARE_CELL_COLOR}")
        println(io, G2_SQUARE_FILL)
        println(io, raw"""\tikzset{
  g2square side/.style={draw=g2RayLabel,line width=.8pt},
  g2square circle/.style={draw=g2RayLabel,fill=white,line width=.8pt},
  g2square cell/.style={fill=g2SquareFill},
  g2square name/.style={text=g2RayLabel,inner sep=0pt,font=\boldmath},
  g2square cell name/.style={g2square name,text=g2CellLabel},""")
        println(io, "  g2triangle letter/.style={text=g2RayLabel,inner sep=0pt," *
                    "font={\\sffamily\\bfseries\\fontsize{$G2_TRIANGLE_LETTER_FONT_PT}" *
                    "{$G2_TRIANGLE_LETTER_FONT_PT}\\selectfont}},")
        println(io, "}")
        println(io, raw"\begin{tikzpicture}[x=1cm,y=1cm,baseline=(current bounding box.center)]")
        square = join(point.(places), "--") * "--cycle;"
        println(io, "\\fill[g2square cell] ", square)
        println(io, "\\draw[g2square side] ", square)
        for (ray, place, corner) in zip(rays, places, corners)
            println(io, "% Ray $ray")
            println(io, "\\draw[g2square circle] $(point(place)) " *
                        "circle[radius=$(g2_simple_number(radius))];")
            # The label touches the circle on the diagonal, away from the square.
            outward = (radius + G2_SQUARE_NAME_GAP_CM) / sqrt(2) .* sign.(place)
            anchor = (place[2] > 0 ? "south" : "north") * " " * (place[1] > 0 ? "west" : "east")
            println(io, "\\node[g2square name,anchor=$anchor] at $(point(place .+ outward)) {" *
                        g2_square_table_font("\\($(g2_square_ray_type(ray))\\)") * "};")
            println(io, "\\begin{scope}[shift={$(point(place))}]")
            corner.write(io)
            println(io, raw"\end{scope}")
        end
        println(io, "% The cell")
        name_y = name_below ? -half - G2_SQUARE_NAME_GAP_CM : -half + G2_SQUARE_NAME_GAP_CM
        println(io, "\\node[g2square cell name,anchor=$(name_below ? "north" : "south")] at " *
                    "$(point((0.0, name_y))) {" *
                    g2_square_table_font("\\oldstylenums{$number}") * "};")
        middle.write(io)
        # A document draws over the figure, in its coordinates with the centre
        # of the square at the origin, by defining this macro.
        println(io, raw"\ifdefined\gtwocelloverlay\gtwocelloverlay\fi")
        println(io, raw"\end{tikzpicture}")
        println(io, raw"\endgroup")
    end
    description = selected.description * ", orbit $number, with its rays clockwise from the top " *
                  "left corner; G2-trees with lettered inner vertices, generated by " *
                  "julia/to_tex/lettered_square_figure.jl"
    return g2_wrap_tree_tikz(tikz, description; standalone=standalone)
end

function g2_lettered_square_main(args=ARGS)
    script = "julia --project=julia julia/to_tex/lettered_square_figure.jl"
    if "--help" in args
        println("Usage: $script [--rays=I,J,K,L] [--scale=FACTOR] [--taxon-scale=FACTOR] " *
                "[--bounded-length=LENGTH] [--ray-scale=FACTOR] [--ray-bounded-length=LENGTH] " *
                "[--name-below] [--cell-letters] [--standalone] [--output=FILE]")
        return nothing
    end
    valued = ("rays", "scale", "taxon-scale", "bounded-length", "ray-scale", "ray-bounded-length", "output")
    for arg in args
        arg in ("--standalone", "--name-below", "--cell-letters") || any(name -> startswith(arg, "--$name="), valued) ||
            throw(ArgumentError("unknown option $arg; use --help"))
    end
    g2_tree_option(args, "output") == "" && throw(ArgumentError("output path must not be empty"))
    rays = g2_tree_option(args, "rays")
    cell = rays === nothing ? G2_SQUARE_RAYS : Tuple(parse.(Int, split(rays, ',')))
    number(name, default) = parse(Float64, something(g2_tree_option(args, name), default))
    length_option(name) = (value = g2_tree_option(args, name); value === nothing ? nothing : parse(Float64, value))
    tex = g2_lettered_square_tex(
        cell;
        scale=number("scale", "1.0"),
        taxon_scale=number("taxon-scale", "0.75"),
        bounded_length=length_option("bounded-length"),
        ray_scale=number("ray-scale", "0.8"),
        ray_bounded_length=length_option("ray-bounded-length"),
        name_below="--name-below" in args,
        cell_letters="--cell-letters" in args,
        standalone="--standalone" in args,
    )
    return g2_tree_emit(args, tex)
end

if abspath(PROGRAM_FILE) == @__FILE__
    g2_lettered_square_main()
end
