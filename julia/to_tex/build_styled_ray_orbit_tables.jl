###############################################################################
# The paper ray-orbit catalog, with the shared G2 tree renderer's options:
# the twelve ray orbits of the tropical G2 Grassmannian, the last of which
# occurs only in characteristic 2.
#
# Run from julia/ with `julia --project=.`. The script prints the tabular, or
# writes it to --output=FILE; main.tex includes it as the PDF of --standalone
# and supplies the float and its caption.
###############################################################################

include(joinpath(@__DIR__, "render_maximal_cell_tree.jl"))

using Printf

const PAPER_RAY_PDF_DIR = normpath(joinpath(@__DIR__, "..", "..", "output", "pdf"))
const PAPER_MAIN_TEX = normpath(joinpath(@__DIR__, "..", "..", "main.tex"))

# Data-orbit indices 1:11 stay fixed in ray_label_orbits.csv; orbit 12 holds
# the rays 52:54, which occur only over GF(2). Article labels follow the
# left-to-right, top-to-bottom order of the catalog. Both are read from
# ray_orbit_labels.csv: the label of each data orbit and its catalog position.
const PAPER_RAY_LABELS_CSV = joinpath(dirname(RAY_ORBITS_CSV), "ray_orbit_labels.csv")
const PAPER_RAY_LABEL_ROWS = sort!(TropicalGeometryG2.csv_rows(PAPER_RAY_LABELS_CSV); by=row -> parse(Int, row["orbit"]))
[parse(Int, row["orbit"]) for row in PAPER_RAY_LABEL_ROWS] == 1:12 &&
    sort([parse(Int, row["catalog_position"]) for row in PAPER_RAY_LABEL_ROWS]) == 1:12 ||
    error("$PAPER_RAY_LABELS_CSV must list the ray orbits 1 to 12 and the catalog positions 1 to 12 once each")
const PAPER_RAY_CATALOG_ORDER =
    Tuple(sortperm([parse(Int, row["catalog_position"]) for row in PAPER_RAY_LABEL_ROWS]))
const PAPER_RAY_ARTICLE_LABELS = Tuple("\\" * row["label"] for row in PAPER_RAY_LABEL_ROWS)
const PAPER_RAY_CHARACTERISTIC_2_ORBIT = 12

function paper_ray_article_label(index::Int)
    1 <= index <= length(PAPER_RAY_ARTICLE_LABELS) ||
        throw(ArgumentError("unknown ray orbit $index"))
    return PAPER_RAY_ARTICLE_LABELS[index]
end

function paper_ray_catalog_entries(orbits)
    length(PAPER_RAY_ARTICLE_LABELS) == length(orbits) &&
        sort(collect(PAPER_RAY_CATALOG_ORDER)) == collect(1:length(orbits)) ||
        error("article labels must cover every ray orbit exactly once")
    return PAPER_RAY_CATALOG_ORDER
end

# Clockwise after the heart at the top. Every bounded split stays contiguous,
# so each tree has a planar straight-line drawing on a regular heptagon.
const PAPER_CATALOG_LEAF_ORDERS = (
    (7, 6, 5, 4, 1, 2, 3),
    (7, 6, 5, 4, 1, 2, 3),
    (7, 5, 4, 3, 2, 1, 6),
    (7, 4, 6, 5, 3, 2, 1),
    (7, 6, 5, 4, 2, 1, 3),
    (7, 5, 6, 4, 1, 3, 2),
    (7, 4, 5, 3, 6, 2, 1),
    (7, 6, 4, 5, 2, 1, 3),
    (7, 4, 5, 6, 3, 2, 1),
    (7, 4, 5, 3, 6, 2, 1),
    (7, 4, 1, 3, 6, 2, 5),
    (7, 4, 1, 3, 5, 2, 6),
)
const PAPER_CATALOG_BOUNDED_EDGE_LENGTH = 0.90
const PAPER_CATALOG_INNER_DISK_RADIUS = 0.95
# The four bounded edges of orbit 12 do not fit the disk of the other panels.
const PAPER_CATALOG_WIDE_DISK_RADIUS = 1.20
const PAPER_CATALOG_LEAF_RADIUS = 2.12

function paper_catalog_tree_coordinates(tree::TreeMetric, orbit_index::Int;
                                        radius::Float64=PAPER_CATALOG_LEAF_RADIUS)
    order = PAPER_CATALOG_LEAF_ORDERS[orbit_index]
    sort(collect(order)) == collect(1:7) || error("invalid leaf order for orbit $orbit_index")
    positions = Dict(taxon => index for (index, taxon) in enumerate(order))
    for edge in rooted_internal_edges(tree)
        slots = sort([positions[taxon] for taxon in edge.cluster])
        slots == collect(first(slots):last(slots)) ||
            error("nonplanar leaf order for orbit $orbit_index and split $(edge.cluster)")
    end

    coordinates = Dict{Int,Tuple{Float64,Float64}}()
    for (position, taxon) in enumerate(order)
        angle = pi / 2 - 2pi * (position - 1) / 7
        coordinates[taxon] = (radius * cos(angle), radius * sin(angle))
    end
    internal = sort([vertex for vertex in keys(tree.adjacency) if vertex > tree.ntaxa])
    lookup = Dict(vertex => index for (index, vertex) in enumerate(internal))
    A = zeros(Float64, length(internal), length(internal))
    bx = zeros(Float64, length(internal))
    by = zeros(Float64, length(internal))
    for vertex in internal
        index = lookup[vertex]
        A[index, index] = length(tree.adjacency[vertex])
        for (neighbor, _) in tree.adjacency[vertex]
            if neighbor <= tree.ntaxa
                x, y = coordinates[neighbor]
                bx[index] += x
                by[index] += y
            else
                A[index, lookup[neighbor]] -= 1
            end
        end
    end
    xs, ys = A \ bx, A \ by
    for vertex in internal
        coordinates[vertex] = (xs[lookup[vertex]], ys[lookup[vertex]])
    end

    if orbit_index == 4
        (positions[1], positions[2]) == (7, 6) ||
            error("taxa 1 and 2 are no longer next to the heart on the left")
    elseif orbit_index == 9
        # First arrange taxa 1,...,6 in mirror-image fans, then attach 7.
        # The unique bounded edge becomes horizontal.
        component_taxa(tree, 8, 9) == [1, 2, 3] ||
            error("orbit 9's split changed")
        coordinates[8] = (-0.75, -0.25)
        coordinates[9] = (0.45, -0.25)
    elseif orbit_index == 10
        # Equal shapes (1/4, 2/5, 6/3) have equal heights on opposite sides.
        component_taxa(tree, 8, 9) == [1, 2, 6] ||
            error("orbit 10's left split changed")
        component_taxa(tree, 10, 9) == [3, 4, 5] ||
            error("orbit 10's right split changed")
        coordinates[8] = (-0.72, -0.23)
        coordinates[9] = (0.0, 0.47)
        coordinates[10] = (0.72, -0.23)
        all(abs(coordinates[left][2] - coordinates[right][2]) < 1e-10
            for (left, right) in ((1, 4), (2, 5), (6, 3))) ||
            error("matching shapes are no longer level in orbit 10")
    elseif orbit_index == 11
        # Leave room for the disks on the two upper pendant edges.
        coordinates[11] = (0.0, 0.10)
    elseif orbit_index == 12
        # A snowflake with hub 12 and the heart on the arm of taxa 1 and 4.
        component_taxa(tree, 11, 12) == [1, 4, 7] &&
            component_taxa(tree, 8, 11) == [1, 4] &&
            component_taxa(tree, 9, 12) == [2, 6] &&
            component_taxa(tree, 10, 12) == [3, 5] ||
            error("orbit 12's splits changed")
        coordinates[12] = (0.0, 0.0)
        coordinates[11] = (0.64, 0.77)
        coordinates[8] = (1.55, 0.35)
        coordinates[9] = (-0.97, 0.26)
        coordinates[10] = (0.26, -0.97)
    end

    # The edge numbers carry the metric data. Give every bounded edge the
    # same drawn length so the diagrams do not suggest a second scale.
    original = copy(coordinates)
    root = only([neighbor for (neighbor, _) in tree.adjacency[7]
                 if neighbor > tree.ntaxa])
    pending = [(root, 0)]
    while !isempty(pending)
        vertex, parent = pop!(pending)
        for (neighbor, _) in tree.adjacency[vertex]
            (neighbor <= tree.ntaxa || neighbor == parent) && continue
            vx, vy = original[vertex]
            nx, ny = original[neighbor]
            dx, dy = nx - vx, ny - vy
            distance = hypot(dx, dy)
            distance > 0 || error("collapsed bounded edge in orbit $orbit_index")
            x, y = coordinates[vertex]
            coordinates[neighbor] =
                (x + PAPER_CATALOG_BOUNDED_EDGE_LENGTH * dx / distance,
                 y + PAPER_CATALOG_BOUNDED_EDGE_LENGTH * dy / distance)
            push!(pending, (neighbor, vertex))
        end
    end
    # Center the whole internal tree in the same small disk in every panel.
    disk_radius = orbit_index == 12 ? PAPER_CATALOG_WIDE_DISK_RADIUS :
                  PAPER_CATALOG_INNER_DISK_RADIUS
    center_x = sum(coordinates[vertex][1] for vertex in internal) / length(internal)
    center_y = sum(coordinates[vertex][2] for vertex in internal) / length(internal)
    for vertex in internal
        x, y = coordinates[vertex]
        coordinates[vertex] = (x - center_x, y - center_y)
        hypot(coordinates[vertex]...) <= disk_radius + 1e-10 ||
            error("internal node left the central disk in orbit $orbit_index")
    end
    if orbit_index == 9
        # Removing the heart leaves a reflection-symmetric six-taxon tree.
        abs(coordinates[8][1] + coordinates[9][1]) < 1e-10 ||
            error("orbit 9's bounded edge is not centered")
        abs(coordinates[8][2] - coordinates[9][2]) < 1e-10 ||
            error("orbit 9's bounded edge is not horizontal")
        all(abs(coordinates[left][1] + coordinates[right][1]) < 1e-10 &&
            abs(coordinates[left][2] - coordinates[right][2]) < 1e-10
            for (left, right) in ((1, 4), (2, 5), (3, 6))) ||
            error("orbit 9's six non-heart taxa lost their symmetry")
    end
    for vertex in internal, (neighbor, _) in tree.adjacency[vertex]
        neighbor > vertex || continue
        distance = hypot(coordinates[vertex][1] - coordinates[neighbor][1],
                         coordinates[vertex][2] - coordinates[neighbor][2])
        abs(distance - PAPER_CATALOG_BOUNDED_EDGE_LENGTH) < 1e-10 ||
            error("bounded edge length changed in orbit $orbit_index")
    end
    for taxon in 1:tree.ntaxa
        x, y = coordinates[taxon]
        abs(hypot(x, y) - radius) < 1e-10 ||
            error("taxon $taxon left the regular heptagon in orbit $orbit_index")
    end
    return coordinates
end

"""Taxa in the component containing `start` after deleting `start--blocked`."""
function component_taxa(tree::TreeMetric, start::Int, blocked::Int)
    result = Int[]
    stack = [(start, blocked)]
    visited = Set{Int}()
    while !isempty(stack)
        vertex, parent = pop!(stack)
        vertex in visited && continue
        push!(visited, vertex)
        vertex <= tree.ntaxa && push!(result, vertex)
        for (neighbor, _) in tree.adjacency[vertex]
            neighbor == parent && continue
            push!(stack, (neighbor, vertex))
        end
    end
    sort!(result)
    return result
end

"""
Load and certify the twelve ray orbits and their exact metrics: the eleven
declared orbits of the 51 rays over QQ, and the orbit of the three rays that
occur only over GF(2).
"""
function styled_ray_orbit_data()
    orbits = load_index_orbits(RAY_ORBITS_CSV, :ray)
    certify_declared_representatives(RAY_ORBITS_CSV, orbits)
    length(orbits) == 11 || error("expected 11 ray orbits")
    sum(length, orbits) == 51 || error("expected 51 rays")
    sort!(vcat(orbits...)) == collect(1:51) ||
        error("ray orbits must partition the 51 ray indices")
    source = load_g2_ray_source(; field=:gf2)
    records = [balance_outer_edges(row_as_vector(source.balanced_rays, index))
               for index in 1:54]
    signature(ray) = rooted_metric_signature(records[ray].tree)
    for (orbit_number, line) in enumerate(readlines(RAY_ORBITS_CSV)[2:end])
        declared = split(chomp(line), ','; keepempty=true)[5]
        representative = first(orbits[orbit_number])
        signature(representative) == declared ||
            error("ray orbit $orbit_number has an incorrect representative")
    end
    # The twelve taxon permutations of W(G2) carry ray 52 to rays 52:54.
    permutations = [parse.(Int, split(split(line, ',')[2], ';'))
                    for line in readlines(LABEL_GROUP_CSV)[2:end]]
    length(permutations) == 12 || error("expected the 12 elements of W(G2)")
    images = Set(rooted_metric_signature(records[52].tree, permutation)
                 for permutation in permutations)
    length(images) == 3 && images == Set(signature.(52:54)) ||
        error("the rays that occur only over GF(2) are not one orbit")
    push!(orbits, collect(52:54))
    return (orbits=orbits, records=records)
end

###############################################################################
# Tables of trees, shared with toric_degeneration_tree_table.jl
###############################################################################

# The color of a table, with the color palette: of the labels of its orbits,
# and of its rules and lines. It is black, as the labels of the maximal trees
# of toric_degeneration_tree_table.jl and \raytype and \celltype in main.tex.
const STYLED_TREE_COLOR = "000000"

"""
Check the drawing options of a table of trees and collect them. `color` is
the HTML color of the table.
"""
function styled_tree_table_options(;
        palette::Symbol=:grayscale,
        label::Bool=false,
        scale::Real=1.0,
        taxon_scale::Real=1.0,
        bounded_length::Real,
        edge_labels::Symbol=:disk,
        color::AbstractString=STYLED_TREE_COLOR)
    palette in (:color, :grayscale) ||
        throw(ArgumentError("palette must be :color or :grayscale"))
    edge_labels in (:disk, :aligned) ||
        throw(ArgumentError("edge_labels must be :disk or :aligned"))
    scale_value = Float64(scale)
    isfinite(5.25 * scale_value) && scale_value > 0 ||
        throw(ArgumentError("scale must be a positive finite number"))
    taxon_scale_value = Float64(taxon_scale)
    isfinite(taxon_scale_value) && taxon_scale_value > 0 ||
        throw(ArgumentError("taxon_scale must be a positive finite number"))
    return (palette=palette, label=label, scale=scale_value,
            taxon_scale=taxon_scale_value,
            bounded_length=Float64(bounded_length),
            edge_labels=edge_labels, color=color)
end

# A tree is drawn at this multiple of the scale of its table.
const PAPER_CATALOG_TREE_SCALE = .78

"""The tree of a panel, drawn by `g2_tree_tikz` with the `options` of its table."""
function styled_tree_tikz(tree::TreeMetric, options; layout_coordinates=nothing)
    return g2_tree_tikz(tree;
        palette=options.palette,
        label=options.label,
        scale=PAPER_CATALOG_TREE_SCALE * options.scale,
        taxon_scale=options.taxon_scale,
        bounded_length=options.bounded_length,
        edge_labels=options.edge_labels,
        layout_coordinates=layout_coordinates,
    )
end

# The labels of the orbits in a table: bold, thickened by an outline of this
# width (pdfrender, loaded by main.tex).
const STYLED_TREE_LABEL_OUTLINE = ".5pt"
# The sizes of the label of an orbit and of its size: at the size of the text,
# so that a heading ends before the heart at the top of a narrow panel.
const STYLED_TREE_LABEL_FONTS = (label=raw"\normalsize", size=raw"\small")

"""
The heading TeX of an orbit, in the color of its table: its label `tex`,
followed by its `size` in parentheses, not bold. Its depth is dropped: a
descender, as of a parenthesis, would lower the tree of its panel.
"""
function styled_tree_label(tex::AbstractString, size::Int, options)
    fonts = STYLED_TREE_LABEL_FONTS
    label = "\\smash[b]{\\textpdfrender{TextRenderingMode=FillStroke," *
            "LineWidth=$STYLED_TREE_LABEL_OUTLINE}{$(fonts.label)\\bfseries\\boldmath $tex}" *
            "{$(fonts.size)($size)}}"
    options.palette == :color || return label
    return "\\textcolor[HTML]{$(options.color)}{$label}"
end

# The thickness of the rule over a panel and of the line between two panels.
const STYLED_TREE_RULE = "1pt"
# The tree of a panel rises by this length in millimetres beside its heading.
const STYLED_TREE_RAISE_MM = 2.5
# The height in centimetres of a panel of the ray catalog at scale 1.
const STYLED_TREE_PANEL_HEIGHT = 4.35

"""
One panel of a table of trees: a rule, the heading TeX `label` of an orbit,
and its tree `tikz`, which rises by `raise` millimetres beside the heading;
a negative `raise` lowers it below the heading line. A `note` stands at the
bottom left of the panel, below the tree: at the bottom of a panel of fixed
height, and on a line of its own in a panel fitted to its tree, whose leaves
would meet it beside the tree. With `note_beside`, the note takes no line: it
stands beside the bottom of the tree, which must leave that corner free.
A `corner` mark stands at the right end of the heading line.
The rule has the color `g2TreeRule` of `styled_tree_grid`.
"""
function styled_tree_panel(io::IO, label::AbstractString, tikz::AbstractString;
                           note::AbstractString="", note_beside::Bool=false,
                           corner::AbstractString="", raise::Real=STYLED_TREE_RAISE_MM)
    println(io, "\\noindent\\textcolor{g2TreeRule}{\\rule{\\linewidth}{$STYLED_TREE_RULE}}\\par")
    println(io, "\\noindent$label", isempty(corner) ? "" : "\\hfill$corner", "\\par")
    @printf(io, "\\vskip %gmm\n", -raise)
    println(io, raw"{\centering")
    print(io, tikz)
    println(io, raw"\par}")
    isempty(note) && return
    # This glue takes the room left in a panel of fixed height.
    println(io, raw"\vskip 0pt plus 1filll")
    note_beside ? println(io, "\\nointerlineskip\\noindent\\smash{\\rlap{$note}}\\par") :
                  println(io, "\\noindent\\rlap{$note}\\par")
end

"""
The width in centimetres that no tree of a table exceeds: the taxa lie on the
heptagon of `g2_tree_tikz`, and half a symbol reaches up to 3 mm beyond it at
symbol scale 1. A millimetre is left for the outlines.
"""
function styled_tree_width(options)
    radius = G2_TREE_UNIT_CM * G2_TREE_LEAF_RADIUS * PAPER_CATALOG_TREE_SCALE * options.scale
    symbol = .3 * G2_TREE_SYMBOL_SCALE * options.taxon_scale
    return 2 * (radius * cos(pi / 14) + symbol) + .1
end

"""
The tabular of a table of trees: `count` panels in rows of `columns`, where
`panel(io, position)` writes the panel at `position`. The panels have the
width and the height of the ray catalog, which leave room around a tree. A
`tight` table fits them to their trees instead: a panel is as wide as the
trees, and at least `minimum_width` centimetres for its heading, and it is as
high as its content; its rows are `row_gap` centimetres apart, by default a
quarter of the scale. `groups` has the group of every panel: a vertical line
separates two panels of a row that belong to different groups. The rules of
the panels and the lines have the color of the table, as `g2TreeRule`.
"""
function styled_tree_grid(panel, io::IO, count::Int, options; columns::Int=3,
                          tight::Bool=false, minimum_width::Real=0, groups=nothing,
                          row_gap::Union{Real,Nothing}=nothing)
    # Keep the table geometry proportional to the tree renderer's scale.
    width = tight ? max(styled_tree_width(options), minimum_width) : 5.25 * options.scale
    gap = (tight ? .25 : .35) * options.scale
    column_spec = @sprintf("p{%gcm}", width)
    gap_spec = @sprintf("@{\\hspace{%gcm}}", gap)
    # The line stands in the middle of the gap, from the rules of the panels
    # to the bottom of their row.
    half_gap = @sprintf("\\hspace{\\dimexpr%gcm-%s/2\\relax}", gap / 2, STYLED_TREE_RULE)
    line_spec = "@{$(half_gap){\\color{g2TreeRule}\\vrule width $STYLED_TREE_RULE " *
                "height $STYLED_TREE_RULE}$(half_gap)}"
    lined(position) = groups !== nothing && position < count &&
                      groups[position] != groups[position + 1]
    minipage = tight ? raw"\begin{minipage}[t]{\linewidth}" :
               @sprintf("\\begin{minipage}[t][%gcm][t]{\\linewidth}",
                        STYLED_TREE_PANEL_HEIGHT * options.scale)
    print(io, raw"""\begingroup
\setlength{\tabcolsep}{0pt}
""")
    println(io, options.palette == :color ? "\\definecolor{g2TreeRule}{HTML}{$(options.color)}" :
                                            raw"\definecolor{g2TreeRule}{gray}{0}")
    # The strut of a tabular row would stand above every panel.
    tight && println(io, raw"\renewcommand{\arraystretch}{0}")
    println(io, "\\begin{tabular}{@{}", join(fill(column_spec, columns), gap_spec), "@{}}")
    for first_position in 1:columns:count
        for column in 1:columns
            position = first_position + column - 1
            if position <= count
                # The cell brings its own column, whose gap has the line.
                line = column < columns && lined(position)
                line && println(io, "\\multicolumn{1}{", column == 1 ? "@{}" : "",
                                column_spec, line_spec, "}{%")
                println(io, minipage)
                panel(io, position)
                println(io, raw"\end{minipage}", line ? "}" : "")
            else
                println(io, raw"\mbox{}")
            end
            column < columns && println(io, " & ")
        end
        if !tight
            @printf(io, "\\\\[%gmm]\n", .2 * options.scale)
        elseif first_position + columns <= count
            @printf(io, "\\\\\\noalign{\\vskip %gcm}\n", something(row_gap, .25 * options.scale))
        else
            println(io, "\\\\")
        end
    end
    println(io, raw"\end{tabular}")
    println(io, raw"\endgroup")
end

"""A complete document: the preamble of main.tex and the `fragment`."""
function styled_tree_paper_document(fragment::AbstractString)
    paper = read(PAPER_MAIN_TEX, String)
    pieces = split(paper, raw"\begin{document}"; limit=2)
    length(pieces) == 2 || error("cannot locate the main.tex preamble")
    return pieces[1] * "\\begin{document}\n" * fragment *
           "\\end{document}\n"
end

"""
Compile the preview `stem`.pdf of a table of trees: a `standalone` document,
or a fragment, which is set in a float with this `caption`.
"""
function styled_tree_compile_preview(tex::AbstractString; standalone::Bool=false,
                                    stem::AbstractString, caption::AbstractString)
    latexmk = Sys.which("latexmk")
    latexmk === nothing && error("latexmk is required for PDF previews")
    mkpath(PAPER_RAY_PDF_DIR)
    pdf_path = joinpath(PAPER_RAY_PDF_DIR, "$stem.pdf")
    mktempdir() do build_dir
        preview_tex = joinpath(build_dir, "$stem.tex")
        if standalone
            write(preview_tex, tex)
        else
            write(joinpath(build_dir, "catalog_fragment.tex"), tex)
            float = raw"""\begin{table}[p]
\centering
\captionsetup{skip=3pt}
\input{catalog_fragment.tex}
""" * "\\caption{$caption}\n\\end{table}\n"
            write(preview_tex, styled_tree_paper_document(float))
        end
        command = `$(latexmk) -silent -pdf -interaction=nonstopmode -halt-on-error -file-line-error $stem.tex`
        run(Cmd(command; dir=build_dir))
        cp(joinpath(build_dir, "$stem.pdf"), pdf_path; force=true)
    end
    return pdf_path
end

"""
The command line of a table of trees. `table_tex(; options...)` returns the
table; `script` is its file, `bounded_length` its default bounded length, and
`stem` and `caption` are those of its preview. The script reads its own
`flags`, and its own `options` with a value, named without their dashes, from
`args`.
"""
function styled_tree_table_main(table_tex, args; script::AbstractString,
                                bounded_length::Real, stem::AbstractString,
                                caption::AbstractString, flags=(), options=())
    if "--help" in args
        println("Usage: julia --project=. to_tex/$script " *
                join("[$flag] " for flag in flags) *
                join("[--$name=$(uppercase(name))] " for name in options) *
                "[--color|--grayscale] " *
                "[--label=true|false] [--scale=FACTOR] [--taxon-scale=FACTOR] " *
                "[--bounded-length=LENGTH] [--edge-labels=disk|aligned] " *
                "[--standalone] [--pdf] [--output=FILE]")
        println("--scale scales tree distances, column widths, row heights, and gaps.")
        println("--pdf compiles a preview to ../output/pdf/ instead of printing the TeX.")
        return nothing
    end
    allowed = ("--color", "--grayscale", "--label=true", "--label=false",
               "--standalone", "--pdf", flags...)
    valued = ("scale", "taxon-scale", "bounded-length", "edge-labels", "output", options...)
    for arg in args
        arg in allowed || any(name -> startswith(arg, "--$name="), valued) ||
            throw(ArgumentError("unknown option $arg; use --help"))
    end
    "--color" in args && "--grayscale" in args &&
        throw(ArgumentError("choose either color or grayscale"))
    "--label=true" in args && "--label=false" in args &&
        throw(ArgumentError("choose one label setting"))
    number(name, default) = parse(Float64, something(g2_tree_option(args, name), default))
    standalone = "--standalone" in args
    tex = table_tex(;
        palette="--color" in args ? :color : :grayscale,
        label="--label=true" in args,
        scale=number("scale", "1.0"),
        taxon_scale=number("taxon-scale", "1.0"),
        bounded_length=number("bounded-length", string(bounded_length)),
        edge_labels=Symbol(something(g2_tree_option(args, "edge-labels"), "disk")),
        standalone=standalone,
    )
    "--pdf" in args || return g2_tree_emit(args, tex)
    g2_tree_option(args, "output") === nothing || g2_tree_emit(args, tex)
    println("PDF: ", styled_tree_compile_preview(tex; standalone=standalone,
                                                stem=stem, caption=caption))
    return nothing
end

###############################################################################
# The ray catalog
###############################################################################

# The panels of the ray catalog: as high as their trees, and at least
# `minimum_width` centimetres wide, which spreads four of them over the text
# width of main.tex. A tree, in the middle of its panel, rises by `raise`
# millimetres, up to the rule above it, beside its heading at the left; the
# rows are `row_gap` centimetres apart.
const PAPER_RAY_CATALOG_PANEL = (minimum_width=3.9, raise=4.0, row_gap=0.2)

function styled_ray_catalog_entry(io::IO, data, index::Int, options)
    orbit = data.orbits[index]
    tree = data.records[first(orbit)].tree
    layout = paper_catalog_tree_coordinates(tree, index)
    styled_tree_panel(
        io, styled_tree_label("\\($(paper_ray_article_label(index))\\)", length(orbit), options),
        styled_tree_tikz(tree, options; layout_coordinates=layout);
        note=index == PAPER_RAY_CHARACTERISTIC_2_ORBIT ? raw"{\small char\ 2}" : "",
        note_beside=true, raise=PAPER_RAY_CATALOG_PANEL.raise,
    )
end

"""
    styled_ray_orbit_tables_tex(; palette=:grayscale, label=false, scale=1.0,
                                taxon_scale=1.0, bounded_length=0.90,
                                edge_labels=:disk, columns=3, standalone=false)

Return the input-ready ray catalog, in rows of `columns` panels. The drawing
options are those of `g2_tree_tikz`; `scale` also scales the table columns,
row heights, and gaps. Pass `standalone=true` for the document of
`g2_standalone_figure`.
"""
function styled_ray_orbit_tables_tex(; columns::Int=3, standalone::Bool=false, drawing...)
    columns >= 1 || throw(ArgumentError("the catalog has at least one column"))
    options = styled_tree_table_options(;
        bounded_length=PAPER_CATALOG_BOUNDED_EDGE_LENGTH, drawing...,
    )
    data = styled_ray_orbit_data()
    entries = paper_ray_catalog_entries(data.orbits)
    fragment = sprint() do io
        println(io, "% Generated by build_styled_ray_orbit_tables.jl from exact G2 ray metrics.")
        styled_tree_grid(io, length(entries), options; columns, tight=true,
                         minimum_width=PAPER_RAY_CATALOG_PANEL.minimum_width,
                         row_gap=PAPER_RAY_CATALOG_PANEL.row_gap,
                         groups=entries .== PAPER_RAY_CHARACTERISTIC_2_ORBIT) do io, position
            styled_ray_catalog_entry(io, data, entries[position], options)
        end
    end
    return standalone ? g2_standalone_figure(fragment) : fragment
end

styled_ray_tables_main(args=ARGS) = styled_tree_table_main(
    (; drawing...) -> styled_ray_orbit_tables_tex(;
        columns=parse(Int, something(g2_tree_option(args, "columns"), "3")), drawing...),
    args;
    script="build_styled_ray_orbit_tables.jl",
    options=("columns",),
    bounded_length=PAPER_CATALOG_BOUNDED_EDGE_LENGTH,
    stem="g2_ray_orbits_styled_catalog",
    caption="Representatives of the 12 Grassmannian ray orbits; " *
            "the last one occurs only in characteristic 2.",
)

if abspath(PROGRAM_FILE) == @__FILE__
    styled_ray_tables_main()
end
