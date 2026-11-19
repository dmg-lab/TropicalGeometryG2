###############################################################################
# A standalone copy of ../figures/toric_degeneration_trees.tex, the figure
# fig:atlas_toric_cells_prime of main.tex, in which an edge and its length are
# green when the edge has a prime cone over it in the tropical G2 flag variety,
# the edges with an apple in toric_degeneration_apple_trees.jl, and red when it
# bears no fruit. At the top right of every panel stands its number of prime
# cones, followed by a green apple.
#
# The panels are drawn with the options of the figure in main.tex, and the
# script checks that, without these colors, they are
# ../figures/toric_degeneration_trees.tex.
#
#   julia --project=. to_tex/toric_degeneration_prime_edge_trees.jl \
#     [--apple-style=flat|outlined|glossy] [--fragment|--standalone] [--pdf] [--output=FILE]
#
# With --fragment the script writes the figure of main.tex, which supplies the
# float and the caption; with --standalone, that figure as the document of
# standalone_figure.jl, whose PDF main.tex includes. With --pdf the page is compiled to
# ../output/pdf/g2_toric_degeneration_prime_edge_trees.pdf, and the fragment, in
# the preamble of main.tex, to g2_toric_degeneration_prime_edge_trees_fragment.pdf.
###############################################################################

include(joinpath(@__DIR__, "toric_degeneration_apple_trees.jl"))

# The figure of main.tex and its drawing options, from ../figures/generate_figures.sh.
const PRIME_EDGE_PAPER_FIGURE = "figures/toric_degeneration_trees.tex"
const PRIME_EDGE_OPTIONS = (scale=0.75, taxon_scale=0.6, bounded_length=1.0)
# The edges with a prime cone over them, and the others.
const PRIME_EDGE_COLOR = "2E7D32"
const BARREN_EDGE_COLOR = "C62828"
# The scale of the apple of decorated_trees_figure.jl at the top right of a panel.
const PRIME_EDGE_APPLE_SCALE = 0.62

"""
The TikZ of the panel tree. With `recolor`, its `prime` edges and their
lengths have the color `g2Prime`, the others `g2Barren`. The lengths stand in
the order of the edges of `g2_tree_commands`.
"""
function prime_edge_tree_tikz(tree, prime, options; recolor::Bool=true)
    tikz = styled_tree_tikz(tree, options)
    recolor || return tikz
    lines = String.(split(tikz, '\n'))
    edges = [(left, right) for left in sort!(collect(keys(tree.adjacency)))
             for (right, _) in tree.adjacency[left] if left < right]
    sort!(edges; by=edge -> (edge[1] <= tree.ntaxa || edge[2] <= tree.ntaxa, edge...))
    weights = findall(line -> startswith(line, "\\node[g2tree weight,"), lines)
    length(weights) == length(edges) || error("the tree does not have one length per edge")
    marked = Set(minmax(v, w) for (v, w) in prime)
    for (index, (left, right)) in enumerate(edges)
        color = (left, right) in marked ? "g2Prime" : "g2Barren"
        kind = left <= tree.ntaxa ? "pendant" : "internal"
        drawn = findall(==("\\draw[g2tree $kind] (v$left)--(v$right);"), lines)
        length(drawn) == 1 || error("the edge ($left, $right) is not drawn once")
        lines[only(drawn)] = "\\draw[g2tree $kind,draw=$color] (v$left)--(v$right);"
        lines[weights[index]] = replace(lines[weights[index]],
                                        "\\node[g2tree weight," => "\\node[g2tree weight,text=$color,")
    end
    return join(lines, '\n')
end

"""
The panel of a row, as toric_tree_table_panel draws it, with its edges
colored by `prime`, and, with `recolor`, the number of prime edges, in the
green of these edges, and a green apple at its top right.
"""
function prime_edge_tree_panel(io::IO, row, options, prime; recolor::Bool=true)
    fields = first.(row.orbits)
    square = length(row.support) == 4 ? TORIC_TREE_TABLE_SQUARE : ""
    apples = "\\textcolor{g2Prime}{$(STYLED_TREE_LABEL_FONTS.label)$(length(prime))}\\,\\smash[b]{\\tikz[baseline=0pt]\\pic[scale=$PRIME_EDGE_APPLE_SCALE] at " *
             "(0,$(g2_simple_number(5.8 * PRIME_EDGE_APPLE_SCALE))mm) {g2apple};}"
    styled_tree_panel(
        io, styled_tree_label("\\oldstylenums{$(toric_table_label_number(row.label))}",
                              row.orbit_size, options),
        prime_edge_tree_tikz(row.tree, prime, options; recolor);
        note=length(fields) == length(TORIC_TABLE_FIELDS) ? "" : TORIC_TREE_TABLE_NOTE[only(fields)],
        corner=recolor ? (isempty(square) ? "" : square * "\\quad") * apples : square,
        raise=TORIC_TREE_TABLE_RAISE_MM,
    )
end

"""
The caption of the figure in main.tex, whose macros it uses, with the
numbers of prime cones in `totals`.
"""
prime_edge_caption(totals) =
    "The \$20\$ types of maximal \$\\Ggroup\$-trees giving toric degenerations, see " *
    "Table~\\ref{tab:g2-grassmannian-newton-okounkov}. An edge and its length are green when " *
    "it has a prime cone of \$\\trop(\\flagspace)\$ over it, and red when it bears no fruit. " *
    "The number at the top right of a tree counts its prime cones. Weighted by the orbit " *
    "sizes, they number \$$(totals[:qq])\$ in characteristic \$\\neq 2\$ and " *
    "\$$(totals[:gf2])\$ in characteristic \$2\$."

"""
    toric_degeneration_prime_edge_trees_tex(; fragment=false)

The standalone page: the figure of main.tex in black, its prime edges in
green and the others in red, the number of prime edges and a green apple at
the top right of every panel, and its caption with a sentence on the colors.
With `fragment`, the figure of main.tex instead, with its colors
in a group and without its caption, and the caption for main.tex. The apples
are drawn in `apple_style`, one of `G2_APPLE_STYLES`.
"""
function toric_degeneration_prime_edge_trees_tex(; fragment::Bool=false,
                                                 apple_style::Symbol=G2_APPLE_STYLE)
    options = styled_tree_table_options(; palette=:color, PRIME_EDGE_OPTIONS...,
                                        color=TORIC_TABLE_LABEL_COLOR)
    censuses = Dict(field => apple_flag_census(field) for field in TORIC_TABLE_FIELDS)
    rows = toric_tree_table_rows()
    primes = [apple_prime_edges(row, censuses) for row in rows]
    totals = Dict(field => sum(length(prime) * row.orbit_size
                               for (row, prime) in zip(rows, primes) if field in first.(row.orbits))
                  for field in TORIC_TABLE_FIELDS)
    layout = TORIC_TREE_TABLE_LAYOUT.toric
    table(options, recolor) = sprint() do io
        styled_tree_grid(io, length(rows), options; tight=true, layout...,
                         row_gap=TORIC_TREE_TABLE_ROW_GAP_CM,
                         groups=[first.(row.orbits) for row in rows]) do io, position
            prime_edge_tree_panel(io, rows[position], options, primes[position]; recolor)
        end
    end
    # Without these colors, the table is the figure of main.tex after its
    # comment line.
    paper = read(joinpath(TORIC_TABLE_ROOT, PRIME_EDGE_PAPER_FIGURE), String)
    last(split(paper, '\n'; limit=2)) == table(options, false) ||
        error("$PRIME_EDGE_PAPER_FIGURE was drawn with other options or data; regenerate it")
    @info "prime edges per panel" Dict(row.label => length(prime) for (row, prime) in zip(rows, primes))
    @info "prime cones of the flag variety" totals
    colors = """
\\definecolor{g2Prime}{HTML}{$PRIME_EDGE_COLOR}
\\definecolor{g2Barren}{HTML}{$BARREN_EDGE_COLOR}
$(g2_apple_colors_tex())\\tikzset{
$(g2_apple_pic(apple_style))
}
"""
    fragment && return (tex="""
% Generated by julia/to_tex/toric_degeneration_prime_edge_trees.jl --fragment: the
% figure $PRIME_EDGE_PAPER_FIGURE in black, with the edges that have a prime cone
% over them in the tropical G2 flag variety, and their lengths, in green, the other
% edges in red, and the number of prime cones at the top right of a panel.
\\begingroup
$(colors)$(table(options, true))\\endgroup
""", caption=prime_edge_caption(totals))
    width = layout.columns * max(styled_tree_width(options), layout.minimum_width) +
            (layout.columns - 1) * .25 * options.scale
    return (tex="""
% Generated by julia/to_tex/toric_degeneration_prime_edge_trees.jl: the figure
% $PRIME_EDGE_PAPER_FIGURE of main.tex in black, with the edges that have a prime
% cone over them in the tropical G2 flag variety, and their lengths, in green, the
% other edges in red, and the number of prime cones at the top right of a panel.
\\documentclass[12pt,border=5mm]{standalone}
\\usepackage{amsmath,amssymb}
\\usepackage{xcolor}
\\usepackage{tikz}
\\usepackage{pdfrender}
$(colors)\\begin{document}
\\begin{tabular}{@{}p{$(g2_simple_number(width))cm}@{}}
$(table(options, true))\\\\[4mm]
The \$20\$ types of maximal \$\\mathrm{G}_2\$-trees giving toric degenerations, see Table~2.
An edge and its length are \\textcolor{g2Prime}{green} when it has a prime cone of
\$\\operatorname{trop}(\\mathcal{F})\$ over it, that is, a toric degeneration of the
\$\\mathrm{G}_2\$ flag variety, and \\textcolor{g2Barren}{red} when it bears no fruit.
The number at the top right of a tree counts its prime cones, its green edges.
Weighted by the orbit sizes, the green edges number
\$$(totals[:qq])\$ in characteristic \$\\neq 2\$ and \$$(totals[:gf2])\$ in characteristic \$2\$.
\\end{tabular}
\\end{document}
""", caption="")
end

function toric_degeneration_prime_edge_trees_main(args=ARGS)
    if "--help" in args
        println("Usage: julia --project=. to_tex/toric_degeneration_prime_edge_trees.jl " *
                "[--apple-style=flat|outlined|glossy] [--fragment|--standalone] [--pdf] [--output=FILE]")
        println("--fragment writes the figure for \\input in main.tex; with --pdf its preview has ",
                "the preamble of main.tex and the caption of the figure.")
        println("--standalone writes that figure as the document of g2_standalone_figure.")
        return nothing
    end
    for arg in args
        arg in ("--pdf", "--fragment", "--standalone") || startswith(arg, "--output=") ||
            startswith(arg, "--apple-style=") ||
            throw(ArgumentError("unknown option $arg; use --help"))
    end
    standalone = "--standalone" in args
    fragment = "--fragment" in args
    fragment && standalone && throw(ArgumentError("choose either --fragment or --standalone"))
    figure = toric_degeneration_prime_edge_trees_tex(; fragment=fragment || standalone,
                                                     apple_style=g2_apple_style_option(args))
    standalone && (figure = (tex=g2_standalone_figure(figure.tex), caption=figure.caption))
    "--pdf" in args || return g2_tree_emit(args, figure.tex)
    g2_tree_option(args, "output") === nothing || g2_tree_emit(args, figure.tex)
    println("PDF: ", styled_tree_compile_preview(figure.tex; standalone=!fragment,
                                                stem="g2_toric_degeneration_prime_edge_trees" *
                                                     (fragment ? "_fragment" : standalone ? "_standalone" : ""),
                                                caption=figure.caption))
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    toric_degeneration_prime_edge_trees_main()
end
