###############################################################################
# A standalone page of the G2-trees of toric_degeneration_tree_table.jl, in
# color and without edge lengths, with an apple on every edge that has a prime
# cone in the tropical G2 flag variety: the flag cell over the support cell of
# the panel whose decorated trees carry the apple, the free taxon 8, on that edge.
#
# The flag cell of an edge is that of the G2-tree atlas, and it is prime when
# the stage-6 flag census of a field of the panel says so. The fields agree on
# every edge, and the script checks that they do. It also checks that the edge
# lengths of the atlas are those of the trees of the panels, location by
# location and up to one shift of all pendant edges, a lineality vector, so
# that the vertex numbers of the atlas are those of the trees.
#
#   julia --project=. to_tex/toric_degeneration_apple_trees.jl \
#     [--scale=FACTOR] [--taxon-scale=FACTOR] [--bounded-length=LENGTH] \
#     [--apple-scale=FACTOR] [--apple-style=flat|outlined|glossy] [--pdf] [--output=FILE]
#
# With --pdf the page is compiled to ../output/pdf/g2_toric_degeneration_apple_trees.pdf.
###############################################################################

include(joinpath(@__DIR__, "toric_degeneration_tree_table.jl"))
include(joinpath(@__DIR__, "decorated_trees_figure.jl"))

# The coarsest flag fans and their stage-6 censuses, by field.
const APPLE_FLAG_DATA = (qq="julia/data/g2/flag", gf2="julia/data/g2/flag_gf2")
# The page is larger than the figure of main.tex, to leave room for the apples.
const APPLE_TREE_SCALE = 1.6
const APPLE_TAXON_SCALE = 0.75
const APPLE_BOUNDED_LENGTH = 0.85
# The apples, relative to those of decorated_trees_figure.jl.
const APPLE_SCALE = 0.8
# The apple at the heading of a panel, before the number of apples.
const APPLE_HEADING_SCALE = 0.62
# The clearance in millimetres that an apple aims for; beyond it, an apple
# prefers the middle of its edge, and then to hang straight down.
const APPLE_TARGET_CLEARANCE_MM = 0.4
# The distance in millimetres from the stem of an apple to the inner vertices,
# beyond their dots.
const APPLE_STEM_VERTEX_GAP_MM = 1.2

"""
The coarsest flag cells of a field by their rays, and whether each is prime by
the census. The census lists the rays of the representative of an orbit; they
must be those of the cell of that number.
"""
function apple_flag_census(field::Symbol)
    directory = APPLE_FLAG_DATA[field]
    ints(text) = parse.(Int, split(text, ';'))
    cell_of = Dict(Set(ints(row["extremal_original_rays"])) => parse(Int, row["coarsest_cell"])
                   for row in toric_table_csv_rows("$directory/02_coarsest_cells.csv"))
    prime = Dict{Int,Bool}()
    for row in toric_table_csv_rows("$directory/06_canonical_interior_degenerations.csv")
        get(cell_of, Set(ints(row["rays"])), 0) == parse(Int, row["representative_cell"]) ||
            error("$field flag census orbit $(row["orbit"]) uses another cell numbering")
        foreach(cell -> prime[cell] = row["prime"] == "true", ints(row["cells"]))
    end
    sort!(collect(keys(prime))) == sort!(collect(values(cell_of))) ||
        error("the $field flag census orbits do not partition the coarsest cells")
    return (cell_of=cell_of, prime=prime)
end

"""The edge of `tree` at an atlas `location`, `(:ray, k)` or `(:edge, i, j)`."""
function apple_atlas_edge(tree, location::AbstractString)
    pendant = match(r"^\(:ray, (\d+)\)$", location)
    if pendant !== nothing
        taxon = parse(Int, pendant[1])
        return (only(tree.adjacency[taxon])[1], taxon)
    end
    bounded = match(r"^\(:edge, (\d+), (\d+)\)$", location)
    bounded === nothing && error("unknown atlas location $location")
    return (parse(Int, bounded[1]), parse(Int, bounded[2]))
end

"""
The edges of the tree of `row` whose flag cells are prime, as pairs (inner
vertex, other end). The atlas lists the flag cells of the edges of its tree
over every field of the support orbit; the fields must agree.
"""
function apple_prime_edges(row, censuses)
    lengths = Dict(minmax(v, w) => QQ(length) for v in keys(row.tree.adjacency)
                   for (w, length) in row.tree.adjacency[v])
    seen = Set{Tuple{Int,Int}}()
    # The atlas tree may differ from the drawn one by the same length on every
    # pendant edge, the all-ones vector of the lineality space.
    shifts = Set{QQFieldElem}()
    prime = Tuple{Int,Int}[]
    for entry in toric_table_csv_rows(TORIC_TABLE_ATLAS)
        entry["label"] == row.label && entry["location_kind"] == "edge" || continue
        edge = apple_atlas_edge(row.tree, entry["location"])
        haskey(lengths, minmax(edge...)) ||
            error("atlas label $(row.label): $(entry["location"]) is not an edge of the drawn tree")
        shift = QQ(parse(Int, entry["length"])) - lengths[minmax(edge...)]
        edge[2] <= row.tree.ntaxa ? push!(shifts, shift) : iszero(shift) ||
            error("atlas label $(row.label): the bounded edge $(entry["location"]) has length " *
                  "$(entry["length"]) in the atlas, and not in the drawn tree")
        push!(seen, minmax(edge...))
        verdicts = Bool[]
        for field in Symbol.(split(entry["fields"], '+'))
            census = censuses[field]
            rays = Set(parse.(Int, split(entry["$(field)_flag_rays"], ';')))
            push!(verdicts, census.prime[census.cell_of[rays]])
        end
        length(unique(verdicts)) == 1 ||
            error("atlas label $(row.label): the fields disagree on $(entry["location"])")
        only(unique(verdicts)) && push!(prime, edge)
    end
    seen == Set(keys(lengths)) ||
        error("the atlas does not list every edge of the tree of label $(row.label)")
    length(shifts) == 1 ||
        error("atlas label $(row.label): the pendant lengths differ from the drawn ones unevenly")
    return prime
end

"""
Where the apples on the `edges` of `tree` hang: the point of its edge from
which an apple hangs, its direction in degrees, and whether its leaf is
mirrored. The tree is drawn at `coordinates`, in centimetres, without its
edge lengths. An apple keeps clear of the edges, the taxa, the inner
vertices, and the other apples, and stays inside the circle of the taxa. Its
stem starts on its edge, as near the middle of the visible part of the edge,
between the dot of its inner vertex and the dot or the taxon at its other
end, as keeps it `APPLE_TARGET_CLEARANCE_MM` clear; there it hangs as close
to straight down as it can.
"""
function apple_placements(tree, coordinates, edges, edges_on, options, apple_scale)
    pt = 2.54 / 72.27
    mm = 0.1
    radius = G2_APPLE_RADIUS_MM * apple_scale * mm
    centre_offset = G2_APPLE_CENTRE_MM * apple_scale * mm
    # The leaf, as a small disk, in the frame of the apple hanging along -y.
    leaf_centre = (1.1, -0.95) .* (apple_scale * mm)
    leaf_radius = 0.45 * apple_scale * mm
    symbol_scale = G2_TREE_SYMBOL_SCALE * options.taxon_scale
    bound = G2_TREE_LEAF_RADIUS * PAPER_CATALOG_TREE_SCALE * options.scale
    segments = [(coordinates[v], coordinates[w], (v <= 7 || w <= 7 ? .3 : 1.) * pt)
                for (v, w) in edges]
    dots = [coordinates[v] for v in keys(tree.adjacency) if v > tree.ntaxa]
    dot_radius = G2_TREE_DOT_RADIUS_PT * pt
    taxon_clearance(p) = minimum(
        g2_minimal_taxon_distance(taxon, (p .- coordinates[taxon]) ./ pt, symbol_scale)
        for taxon in 1:tree.ntaxa) * pt
    # The clearance of a disk of radius `r` at `p` from the drawing, without
    # the edge `own` from which its apple hangs.
    function clearance(p, r, own, others)
        clear = taxon_clearance(p) - r
        for (index, (a, b, width)) in enumerate(segments)
            index == own && continue
            clear = min(clear, g2_segment_distance(p, a, b) - width - r)
        end
        for dot in dots
            clear = min(clear, hypot((p .- dot)...) - dot_radius - r)
        end
        for other in others
            clear = min(clear, hypot((p .- other.centre)...) - radius - r,
                        hypot((p .- other.stem)...) - r)
        end
        return min(clear, bound - hypot(p...) - r)
    end
    rotate((x, y), angle) = (x * cosd(angle + 90) - y * sind(angle + 90),
                             x * sind(angle + 90) + y * cosd(angle + 90))
    steps = 0:.0125:1
    visible(stem, gap) = taxon_clearance(stem) >= .5pt &&
                         all(dot -> hypot((stem .- dot)...) >= dot_radius + gap, dots)
    function apple_candidates(index, others, vertex_gap)
        a, b, _ = segments[index]
        heading = atand(b[2] - a[2], b[1] - a[1])
        span = [t for t in steps if visible(a .+ t .* (b .- a), .5pt)]
        middle = (first(span) + last(span)) / 2
        candidates = []
        for t in steps
            stem = a .+ t .* (b .- a)
            visible(stem, vertex_gap) || continue
            # The distance of the stem from the middle, in half millimetres.
            off = round(abs(t - middle) * hypot((b .- a)...) / (.5mm))
            for angle in 0:5:355
                # An apple leaves its edge at an angle of at least 45 degrees.
                abs(sind(angle - heading)) >= sind(45) || continue
                body = stem .+ centre_offset .* (cosd(angle), sind(angle))
                clear = min(clearance(body, radius, 0, others),
                            g2_segment_distance(body, a, b) - radius - .2mm,
                            minimum((hypot((stem .- other.centre)...) - radius for other in others);
                                    init=Inf))
                for mirrored in (false, true)
                    leaf = stem .+ rotate(leaf_centre .* (mirrored ? -1 : 1, 1), angle)
                    total = min(clear, clearance(leaf, leaf_radius, index, others))
                    turn = abs(mod(angle + 90 + 180, 360) - 180)
                    push!(candidates, ((round(min(total, APPLE_TARGET_CLEARANCE_MM * mm) / mm; digits=1),
                                        -off, -round(turn / 15), round(total / mm; digits=1)),
                                       (; stem, centre=body, angle, mirrored, clear=total)))
                end
            end
        end
        return candidates
    end
    function best(index, others)
        # The stem keeps away from the inner vertices, so that it starts on
        # its edge and not at an end; on a short edge it may come closer.
        for gap in (APPLE_STEM_VERTEX_GAP_MM, APPLE_STEM_VERTEX_GAP_MM / 2, 0.0)
            candidates = apple_candidates(index, others, gap * mm)
            isempty(candidates) || return last(argmax(first, candidates))
        end
        error("no visible point on the edge $(edges[index])")
    end
    # The bounded edges first: the middle of the tree has the least room.
    order = sort(edges_on; by=index -> (edges[index][1] <= 7 || edges[index][2] <= 7, index))
    placed = Dict{Int,Any}()
    for _ in 1:3, index in order
        placed[index] = best(index, [placed[other] for other in keys(placed) if other != index])
    end
    return [placed[index] for index in order]
end

"""The TikZ of the apples on the `prime` edges of `tree`, and their smallest clearance in mm."""
function apple_tree_tikz(tree, prime, options, apple_scale)
    coordinates = g2_minimal_tree_coordinates(tree; scale=PAPER_CATALOG_TREE_SCALE * options.scale,
                                               bounded_length=options.bounded_length)
    edges = [(left, right) for left in sort!(collect(keys(tree.adjacency)))
             for (right, _) in tree.adjacency[left] if left < right]
    on = [findfirst(==(minmax(v, w)), edges) for (v, w) in prime]
    apples = isempty(on) ? [] : apple_placements(tree, coordinates, edges, on, options, apple_scale)
    point((x, y)) = @sprintf("(%.5f,%.5f)", x, y)
    tikz = join(("\\pic[rotate=$(apple.angle + 90),scale=$(g2_simple_number(apple_scale))" *
                 (apple.mirrored ? ",xscale=-1" : "") * "] at $(point(apple.stem)) {g2apple};\n"
                 for apple in apples))
    return (tikz=tikz, clearance=isempty(apples) ? Inf : minimum(apple.clear for apple in apples) * 10)
end

"""The panel of a row, as in toric_tree_table_panel, with its apples and their number."""
function apple_tree_panel(io::IO, row, options, prime, apple_scale)
    fields = first.(row.orbits)
    apples = apple_tree_tikz(row.tree, prime, options, apple_scale)
    tikz = styled_tree_tikz(row.tree, options)
    count(raw"\end{tikzpicture}", tikz) == 1 || error("unexpected tree TikZ")
    tikz = replace(tikz, raw"\end{tikzpicture}" => apples.tikz * raw"\end{tikzpicture}")
    height = 5.8 * APPLE_HEADING_SCALE
    corner = isempty(prime) ? "" :
             "\\textcolor{g2AppleDark}{\\large$(length(prime))}\\," *
             "\\smash[b]{\\tikz[baseline=0pt]\\pic[scale=$APPLE_HEADING_SCALE] at " *
             "(0,$(g2_simple_number(height))mm) {g2apple};}"
    styled_tree_panel(
        io, styled_tree_label("\\oldstylenums{$(toric_table_label_number(row.label))}",
                              row.orbit_size, options),
        tikz;
        note=length(fields) == length(TORIC_TABLE_FIELDS) ? "" : TORIC_TREE_TABLE_NOTE[only(fields)],
        corner=corner,
    )
    return apples.clearance
end

"""
    toric_degeneration_apple_trees_tex(; scale=1.6, taxon_scale=0.75,
                                       bounded_length=0.85, apple_scale=0.8,
                                       apple_style=G2_APPLE_STYLE)

The standalone page: the table of toric_degeneration_tree_table.jl in color,
without edge lengths, with apples, and a caption that counts the prime cones
of the flag variety. The apples are drawn in `apple_style`, one of
`G2_APPLE_STYLES`.
"""
function toric_degeneration_apple_trees_tex(; scale::Real=APPLE_TREE_SCALE,
                                            taxon_scale::Real=APPLE_TAXON_SCALE,
                                            bounded_length::Real=APPLE_BOUNDED_LENGTH,
                                            apple_scale::Real=APPLE_SCALE,
                                            apple_style::Symbol=G2_APPLE_STYLE)
    # The trees are drawn without their edge lengths.
    options = merge(styled_tree_table_options(; palette=:color, scale, taxon_scale, bounded_length,
                                              color=TORIC_TABLE_LABEL_COLOR),
                    (edge_labels=:none,))
    censuses = Dict(field => apple_flag_census(field) for field in TORIC_TABLE_FIELDS)
    rows = toric_tree_table_rows()
    primes = [apple_prime_edges(row, censuses) for row in rows]
    # The prime cones of the flag variety, by field: an apple counts once per
    # cell of its orbit.
    totals = Dict(field => sum(length(prime) * row.orbit_size
                               for (row, prime) in zip(rows, primes) if field in first.(row.orbits))
                  for field in TORIC_TABLE_FIELDS)
    layout = TORIC_TREE_TABLE_LAYOUT.toric
    width = layout.columns * max(styled_tree_width(options), layout.minimum_width) +
            (layout.columns - 1) * .25 * options.scale
    clearances = Float64[]
    table = sprint() do io
        styled_tree_grid(io, length(rows), options; tight=true, layout...,
                         groups=[first.(row.orbits) for row in rows]) do io, position
            push!(clearances, apple_tree_panel(io, rows[position], options, primes[position],
                                               apple_scale))
        end
    end
    @info "apples per panel" Dict(row.label => length(prime) for (row, prime) in zip(rows, primes))
    @info "prime cones of the flag variety" totals
    @info "smallest clearance of an apple in mm" minimum(clearances)
    colors = g2_apple_colors_tex()
    apple = "\\tikz[baseline=0pt]\\pic[scale=$APPLE_HEADING_SCALE] at (0,$(g2_simple_number(5.8 * APPLE_HEADING_SCALE))mm) {g2apple};"
    return """
% Generated by julia/to_tex/toric_degeneration_apple_trees.jl: the trees of
% figures/toric_degeneration_trees.tex, with an apple on every edge whose flag
% cell is a prime cone of the tropical G2 flag variety.
\\documentclass[12pt,border=5mm]{standalone}
\\usepackage{amsmath,amssymb}
\\usepackage{xcolor}
\\usepackage{tikz}
\\usepackage{pdfrender}
$(colors)\\tikzset{
$(g2_apple_pic(apple_style))
}
\\begin{document}
\\begin{tabular}{@{}p{$(g2_simple_number(width))cm}@{}}
$(table)\\\\[4mm]
\\small The \$20\$ types of maximal \$\\mathrm{G}_2\$-trees giving toric degenerations of the
\$\\mathrm{G}_2\$-Grassmannian. An apple $(apple) hangs from an edge when placing the apple
(taxon \$8\$) on that edge gives a prime cone of \$\\operatorname{trop}(\\mathcal{F})\$, that is,
a toric degeneration of the \$\\mathrm{G}_2\$ flag variety. A heading counts the apples of its
tree; weighted by the orbit sizes, they add up to \$$(totals[:qq])\$ in characteristic
\$\\neq 2\$ and to \$$(totals[:gf2])\$ in characteristic \$2\$.
\\end{tabular}
\\end{document}
"""
end

function toric_degeneration_apple_trees_main(args=ARGS)
    script = "julia --project=. to_tex/toric_degeneration_apple_trees.jl"
    if "--help" in args
        println("Usage: $script [--scale=FACTOR] [--taxon-scale=FACTOR] [--bounded-length=LENGTH] " *
                "[--apple-scale=FACTOR] [--apple-style=flat|outlined|glossy] [--pdf] [--output=FILE]")
        return nothing
    end
    valued = ("scale", "taxon-scale", "bounded-length", "apple-scale", "apple-style", "output")
    for arg in args
        arg == "--pdf" || any(name -> startswith(arg, "--$name="), valued) ||
            throw(ArgumentError("unknown option $arg; use --help"))
    end
    number(name, default) = parse(Float64, something(g2_tree_option(args, name), string(default)))
    tex = toric_degeneration_apple_trees_tex(;
        scale=number("scale", APPLE_TREE_SCALE),
        taxon_scale=number("taxon-scale", APPLE_TAXON_SCALE),
        bounded_length=number("bounded-length", APPLE_BOUNDED_LENGTH),
        apple_scale=number("apple-scale", APPLE_SCALE),
        apple_style=g2_apple_style_option(args),
    )
    "--pdf" in args || return g2_tree_emit(args, tex)
    g2_tree_option(args, "output") === nothing || g2_tree_emit(args, tex)
    println("PDF: ", styled_tree_compile_preview(tex; standalone=true,
                                                stem="g2_toric_degeneration_apple_trees", caption=""))
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    toric_degeneration_apple_trees_main()
end
