###############################################################################
# Input-ready trees for the cones of the tropical G2 Grassmannian over QQ or
# over GF(2): a generator of the lineality space, a ray, an edge, or a
# maximal cell.
#
# This file may also be included by an atlas. The public entry points are
# g2_tree_tikz and the fragments g2_..._tree_tex of the trees (plus their file
# writers). The trees g2_lineality_tree, g2_ray_tree, g2_edge_tree and
# g2_maximal_cell_tree, and the sources load_g2_ray_source and
# load_g2_maximal_cell_source, are those of
# studies/g2/atlas/flag/cell_orbit_trees.jl, which this file includes:
# both fields share one ray numbering, rays 1:51 over QQ and 52:54 over GF(2)
# only, and a cell index is the index of the cell in the census of its field.
###############################################################################

include(joinpath(
    @__DIR__, "..", "studies", "g2", "atlas", "flag", "cell_orbit_trees.jl",
))
include(joinpath(@__DIR__, "standalone_figure.jl"))

using Printf

const G2_TREE_LEAF_RADIUS = 2.12
const G2_TREE_INNER_RADIUS = 0.95
const G2_TREE_BOUNDED_DRAW_LENGTH = 0.70
const G2_TREE_UNIT_CM = 1.0
const G2_TREE_SYMBOL_SCALE = 1.35
const G2_TREE_WEIGHT_FONT_PT = 7.5
const G2_TREE_TAXON_FONT_PT = 12.0

"""The field `:qq` or `:gf2` that `name` denotes on the command line."""
function g2_tree_field(name::AbstractString)
    lowercase(name) == "qq" && return :qq
    lowercase(name) in ("gf2", "gf(2)", "f2") && return :gf2
    throw(ArgumentError("unknown field $name; use QQ or GF2"))
end

"""
    g2_heptagon_layout(tree, order)

A layout of `tree` for the `layout_coordinates` of `g2_tree_tikz`: the taxa
on the regular heptagon in the cyclic `order`, clockwise from the top, and
every inner vertex at the mean of its neighbours, as in `symmetric_layout`.
The trees of the faces of a cell can thus be drawn with the taxa where the
tree of the cell has them, `topology_compatible_leaf_order` of that tree.
Every split of `tree` must be contiguous in `order`, so that the drawing is
planar.
"""
function g2_heptagon_layout(tree::TreeMetric, order)
    sort(collect(order)) == collect(1:tree.ntaxa) ||
        throw(ArgumentError("the order must list every taxon once"))
    for edge in rooted_internal_edges(tree)
        # A contiguous cluster has two ends around the heptagon.
        ends = count(index -> (order[index] in edge.cluster) !=
                              (order[mod1(index + 1, length(order))] in edge.cluster),
                     eachindex(order))
        ends == 2 || throw(ArgumentError(
            "the split $(edge.cluster) is not contiguous in the order $(Tuple(order))",
        ))
    end
    coordinates = Dict{Int,Tuple{Float64,Float64}}()
    for (position, taxon) in enumerate(order)
        angle = pi / 2 - 2pi * (position - 1) / tree.ntaxa
        coordinates[taxon] = (G2_TREE_LEAF_RADIUS * cos(angle), G2_TREE_LEAF_RADIUS * sin(angle))
    end
    internal = sort!([vertex for vertex in keys(tree.adjacency) if vertex > tree.ntaxa])
    index = Dict(vertex => i for (i, vertex) in enumerate(internal))
    A = zeros(Float64, length(internal), length(internal))
    bx = zeros(Float64, length(internal))
    by = zeros(Float64, length(internal))
    for vertex in internal
        i = index[vertex]
        A[i, i] = length(tree.adjacency[vertex])
        for (neighbor, _) in tree.adjacency[vertex]
            if neighbor <= tree.ntaxa
                bx[i] += coordinates[neighbor][1]
                by[i] += coordinates[neighbor][2]
            else
                A[i, index[neighbor]] -= 1
            end
        end
    end
    xs, ys = A \ bx, A \ by
    for vertex in internal
        coordinates[vertex] = (xs[index[vertex]], ys[index[vertex]])
    end
    return coordinates
end

"""
Place taxa on a regular heptagon and all inner vertices in a central disk.
Every bounded edge has drawn length `bounded_length` before the overall
`scale`; printed metric lengths remain exact.
"""
function g2_minimal_tree_coordinates(tree::TreeMetric;
                                     scale::Real=1.0,
                                     bounded_length::Real=G2_TREE_BOUNDED_DRAW_LENGTH,
                                     layout_coordinates=nothing)
    scale_value = Float64(scale)
    isfinite(scale_value * G2_TREE_LEAF_RADIUS) && scale_value > 0 ||
        throw(ArgumentError("scale must be a positive finite number"))
    bounded_value = Float64(bounded_length)
    isfinite(bounded_value * scale_value) && bounded_value > 0 ||
        throw(ArgumentError("bounded_length must be a positive finite number"))
    tree.ntaxa == 7 || error("the G2 renderer expects seven taxa")
    original = layout_coordinates === nothing ?
               symmetric_layout(tree; radius=G2_TREE_LEAF_RADIUS) :
               layout_coordinates
    coordinates = copy(original)
    root = only(neighbor for (neighbor, _) in tree.adjacency[7]
                if neighbor > tree.ntaxa)
    internal = Dict{Int,Tuple{Float64,Float64}}(root => (0.0, 0.0))
    pending = [(root, 0)]
    while !isempty(pending)
        vertex, parent = pop!(pending)
        for (neighbor, _) in tree.adjacency[vertex]
            (neighbor <= tree.ntaxa || neighbor == parent) && continue
            dx = original[neighbor][1] - original[vertex][1]
            dy = original[neighbor][2] - original[vertex][2]
            norm = hypot(dx, dy)
            norm > 0 || error("collapsed bounded edge")
            x, y = internal[vertex]
            internal[neighbor] =
                (x + bounded_value * dx / norm,
                 y + bounded_value * dy / norm)
            push!(pending, (neighbor, vertex))
        end
    end
    n = length(internal)
    n == count(vertex -> vertex > tree.ntaxa, keys(tree.adjacency)) ||
        error("the internal tree is disconnected")
    cx = sum(point[1] for point in values(internal)) / n
    cy = sum(point[2] for point in values(internal)) / n
    inner_radius = G2_TREE_INNER_RADIUS *
                   bounded_value / G2_TREE_BOUNDED_DRAW_LENGTH
    for (vertex, (x, y)) in internal
        coordinates[vertex] = (x - cx, y - cy)
        radius = hypot(coordinates[vertex]...)
        radius <= inner_radius + 1e-10 ||
            error("inner vertex $vertex left the central disk")
        radius < G2_TREE_LEAF_RADIUS ||
            throw(ArgumentError("bounded_length places inner vertex $vertex outside the heptagon"))
    end
    for taxon in 1:7
        abs(hypot(coordinates[taxon]...) - G2_TREE_LEAF_RADIUS) < 1e-10 ||
            error("taxon $taxon left the heptagon")
    end
    return Dict(vertex => (scale_value * x, scale_value * y)
                for (vertex, (x, y)) in coordinates)
end

function g2_minimal_tree_styles(io::IO, palette::Symbol)
    palette in (:grayscale, :color) ||
        throw(ArgumentError("palette must be :grayscale or :color"))
    # The edges, the inner vertices, and the edge lengths are black.
    edge = "black"
    weight_text = edge
    dark_fill = palette == :color ? "g2TaxonBlue" : "black!42"
    light_fill = palette == :color ? "g2TaxonYellow" : "white"
    heart_fill = palette == :color ? "g2HeartRed" : "black"
    dark_text = palette == :color ? "white" : "black"
    light_text = palette == :color ? "white" : "black"
    outline = palette == :color ? "white" : "black"
    outline_width = palette == :color ? ".55pt" : "1.2pt"
    aligned_pendant_text = edge
    aligned_internal_text = edge
    println(io, raw"\begingroup")
    if palette == :color
        print(io, raw"""\definecolor{g2TaxonBlue}{HTML}{1676A3}
\definecolor{g2TaxonYellow}{HTML}{E39A24}
\definecolor{g2HeartRed}{HTML}{D7443E}
""")
    end
    print(io, "\\tikzset{\n" *
              "  g2tree pendant/.style={draw=$edge,line width=.6pt,line cap=round},\n" *
              "  g2tree internal/.style={draw=$edge,line width=2pt,line cap=round},\n" *
              "  g2tree vertex/.style={circle,fill=$edge,draw=white,line width=.4pt,inner sep=0pt,minimum size=$(G2_TREE_DOT_SIZE_PT)pt},\n" *
              "  g2tree dark/.style={fill=$dark_fill,draw=$outline,line width=$outline_width,text=$dark_text,inner sep=0pt},\n" *
              "  g2tree light/.style={fill=$light_fill,draw=$outline,line width=$outline_width,text=$light_text,inner sep=0pt},\n" *
              "  g2tree weight/.style={circle,fill=white,draw=none,text=$weight_text,inner sep=.15pt,minimum size=8.5pt},\n" *
              "  g2tree aligned pendant/.style={fill=white,text=$aligned_pendant_text,inner sep=.7pt,sloped,above=1pt,font={\\sffamily\\fontsize{7.5}{8}\\selectfont}},\n" *
              "  g2tree aligned internal/.style={fill=white,text=$aligned_internal_text,inner sep=1pt,sloped,above=1pt,font={\\sffamily\\fontsize{7.5}{8}\\selectfont\\bfseries}},\n" *
              "}\n")
    return (dark_fill=dark_fill, light_fill=light_fill,
            heart_fill=heart_fill, dark_text=dark_text,
            light_text=light_text, outline=outline,
            outline_width=outline_width)
end

function g2_minimal_tree_taxon(io::IO, taxon::Int, colors;
                               label::Bool, taxon_scale::Real=1.0)
    text = label ? string(taxon) : ""
    symbol_scale = G2_TREE_SYMBOL_SCALE * taxon_scale
    font_size = G2_TREE_TAXON_FONT_PT * taxon_scale
    font = @sprintf("\\sffamily\\fontsize{%g}{%g}\\selectfont\\bfseries",
                    font_size, font_size)
    style = taxon <= 3 ? "g2tree dark" : "g2tree light"
    if taxon in (1, 4)
        size = 4.6symbol_scale
        println(io, "\\node[$style,rectangle,minimum size=$(size)mm,font={$font}] " *
                    "at (v$taxon) {$text};")
    elseif taxon in (3, 6)
        size = 4.8symbol_scale
        println(io, "\\node[$style,circle,minimum size=$(size)mm,font={$font}] " *
                    "at (v$taxon) {$text};")
    elseif taxon in (2, 5)
        fill = taxon == 2 ? colors.dark_fill : colors.light_fill
        ink = taxon == 2 ? colors.dark_text : colors.light_text
        width = colors.outline == "white" ? colors.outline_width :
                taxon == 2 ? "1.2pt" : "1pt"
        scoped_font = G2_TREE_TAXON_FONT_PT / G2_TREE_SYMBOL_SCALE
        println(io, "\\begin{scope}[shift={(v$taxon)}," *
                    "scale=$symbol_scale,transform shape]")
        println(io, "\\path[fill=$fill,draw=$(colors.outline),line width=$width,line join=round] " *
                    "(0mm,3.0mm)--(-3.0mm,-2.2mm)--(3.0mm,-2.2mm)--cycle;")
        if label
            println(io, "\\node[font={\\sffamily\\fontsize{$scoped_font}{$scoped_font}" *
                        "\\selectfont\\bfseries},text=$ink,inner sep=0pt] " *
                        "at (0mm,-.35mm) {$taxon};")
        end
        println(io, raw"\end{scope}")
    elseif taxon == 7
        scoped_font = G2_TREE_TAXON_FONT_PT / G2_TREE_SYMBOL_SCALE
        println(io, "\\begin{scope}[shift={(v7)}," *
                    "scale=$symbol_scale,transform shape]")
        heart_outline = colors.outline == "white" ? "white" : colors.heart_fill
        heart_width = colors.outline == "white" ? ".55pt" : ".5pt"
        print(io, "\\path[fill=$(colors.heart_fill),draw=$heart_outline,line width=$heart_width]\n")
        print(io, raw"""  (0mm,-2.7mm)..controls(-.45mm,-2.1mm)and(-3.0mm,.35mm)..(-3.0mm,1.45mm)
  ..controls(-3.0mm,3.2mm)and(-.9mm,3.4mm)..(0mm,1.75mm)
  ..controls(.9mm,3.4mm)and(3.0mm,3.2mm)..(3.0mm,1.45mm)
  ..controls(3.0mm,.35mm)and(.45mm,-2.1mm)..(0mm,-2.7mm)--cycle;
""")
        if label
            println(io, "\\node[font={\\sffamily\\fontsize{$scoped_font}{$scoped_font}" *
                        "\\selectfont\\bfseries},text=white,inner sep=0pt] " *
                        "at (0mm,.35mm) {7};")
        end
        println(io, raw"\end{scope}")
    else
        error("unexpected taxon $taxon")
    end
end

function g2_minimal_edge_label_position(left::Int, right::Int, coordinates;
                                       taxon_scale::Real=1.0,
                                       weight_font_pt::Real=G2_TREE_WEIGHT_FONT_PT)
    if left > 7 && right > 7
        return ((coordinates[left][1] + coordinates[right][1]) / 2,
                (coordinates[left][2] + coordinates[right][2]) / 2)
    end
    leaf = left <= 7 ? left : right
    inner = left <= 7 ? right : left
    lx, ly = coordinates[leaf]
    ix, iy = coordinates[inner]
    distance_pt = hypot(ix - lx, iy - ly) * G2_TREE_UNIT_CM * (72.27 / 2.54)
    shape_radius_mm = leaf == 7 ? 3.4G2_TREE_SYMBOL_SCALE :
                      leaf in (2, 5) ? 3.0G2_TREE_SYMBOL_SCALE :
                      leaf in (3, 6) ? 2.4G2_TREE_SYMBOL_SCALE :
                      2.3G2_TREE_SYMBOL_SCALE
    shape_radius_mm *= taxon_scale
    shape_radius_pt = shape_radius_mm * (72.27 / 25.4)
    disk_radius_pt = (weight_font_pt + 1.0) / 2
    fraction = clamp((shape_radius_pt + disk_radius_pt + .3) /
                     distance_pt, .5, .75)
    return (lx + fraction * (ix - lx), ly + fraction * (iy - ly))
end

# The diameter of an inner vertex dot, in points.
const G2_TREE_DOT_SIZE_PT = 3.8
# The geometry of the edge lengths on disks, in points: the room that they
# keep around a vertex dot, and the size of a digit of the lengths. The room
# is the radius with its outline of the dots of 4.8pt that the figures had
# first, so that the lengths stay where they were.
const G2_TREE_DOT_RADIUS_PT = 2.6
# A lettered inner vertex of `g2_tree_commands`: its diameter, the room that
# the lengths keep around it, as above that of a dot of 9.5pt, and the size of
# its letter, in points.
const G2_TREE_LETTERED_DOT_SIZE_PT = 8.5
const G2_TREE_LETTERED_DOT_RADIUS_PT = 4.95
const G2_TREE_LETTER_FONT_PT = 7
const G2_TREE_DIGIT_WIDTH_PT = .5G2_TREE_WEIGHT_FONT_PT
const G2_TREE_DIGIT_HEIGHT_PT = .7G2_TREE_WEIGHT_FONT_PT
const G2_TREE_PT_PER_MM = 72.27 / 25.4

# The outline of the heart of taxon 7 in mm at symbol scale 1, as a polygon
# through its curves.
const G2_TREE_HEART_OUTLINE = let
    curves = (((0., -2.7), (-.45, -2.1), (-3., .35), (-3., 1.45)),
              ((-3., 1.45), (-3., 3.2), (-.9, 3.4), (0., 1.75)),
              ((0., 1.75), (.9, 3.4), (3., 3.2), (3., 1.45)),
              ((3., 1.45), (3., .35), (.45, -2.1), (0., -2.7)))
    [(1 - t)^3 .* a .+ 3(1 - t)^2 * t .* b .+ 3(1 - t) * t^2 .* c .+ t^3 .* d
     for (a, b, c, d) in curves for t in 0:1/16:15/16]
end

"""The distance from the point `p` to the segment from `a` to `b`."""
function g2_segment_distance(p, a, b)
    d = b .- a
    t = clamp(sum((p .- a) .* d) / max(sum(d .* d), eps()), 0, 1)
    return hypot((p .- a .- t .* d)...)
end

"""The distance from `p` to a `polygon`, negative inside it."""
function g2_polygon_distance(p, polygon)
    n = length(polygon)
    distance = minimum(g2_segment_distance(p, polygon[i], polygon[mod1(i + 1, n)])
                       for i in 1:n)
    inside = false
    for i in 1:n
        (x1, y1), (x2, y2) = polygon[i], polygon[mod1(i + 1, n)]
        (y1 > p[2]) != (y2 > p[2]) &&
            p[1] < x1 + (p[2] - y1) * (x2 - x1) / (y2 - y1) && (inside = !inside)
    end
    return inside ? -distance : distance
end

"""
The distance in points from `p`, relative to the centre of the symbol of
`taxon`, to the outline of the symbol that `g2_minimal_tree_taxon` draws,
negative inside it.
"""
function g2_minimal_taxon_distance(taxon::Int, p, symbol_scale::Real)
    mm = G2_TREE_PT_PER_MM * symbol_scale
    if taxon in (1, 4)
        q = abs.(p) .- 2.3mm
        return hypot(max.(q, 0)...) + min(maximum(q), 0)
    elseif taxon in (3, 6)
        return hypot(p...) - 2.4mm
    end
    polygon = taxon == 7 ? G2_TREE_HEART_OUTLINE : ((0., 3.), (-3., -2.2), (3., -2.2))
    return g2_polygon_distance(p, [mm .* v for v in polygon])
end

"""
The width and the height in points of the TeX of the edge length `q`, in a
font of `font_pt` points.
"""
function g2_minimal_metric_size(q; font_pt::Real=G2_TREE_WEIGHT_FONT_PT)
    q = QQ(q)
    digits(n) = length(string(abs(n)))
    minus = q < 0 ? 1.5 : 0.0
    factor = font_pt / G2_TREE_WEIGHT_FONT_PT
    denominator(q) == 1 &&
        return (factor * G2_TREE_DIGIT_WIDTH_PT * (digits(numerator(q)) + minus),
                factor * G2_TREE_DIGIT_HEIGHT_PT)
    # A fraction in text style stacks two rows of smaller digits.
    rows = max(digits(numerator(q)), digits(denominator(q)))
    return (factor * G2_TREE_DIGIT_WIDTH_PT * (.7rows + minus),
            2factor * G2_TREE_DIGIT_HEIGHT_PT)
end

"""
The positions of the lengths of the `edges` on disks, as `g2_tree_tikz` draws
them. A length stands where `g2_minimal_edge_label_position` puts it, unless
its digits come within a point of a taxon symbol, of a vertex dot, or of
another disk, which covers them when it is drawn later: it then moves along
its edge to the nearest position that keeps them clear, and beside its edge
only when every position on the edge brings them within a third of a point.
Among the positions that keep them clear, its disk avoids covering the other
edges. The vertex dots have the radius `dot_radius_pt`, and the lengths the
font size `weight_font_pt`. With `beside=false` every length stays on its
edge, at the position that keeps it clearest.
"""
function g2_minimal_edge_label_positions(edges, coordinates, ntaxa::Int;
                                         taxon_scale::Real=1.0,
                                         dot_radius_pt::Real=G2_TREE_DOT_RADIUS_PT,
                                         weight_font_pt::Real=G2_TREE_WEIGHT_FONT_PT,
                                         beside::Bool=true)
    unit = G2_TREE_UNIT_CM * 72.27 / 2.54
    point(vertex) = unit .* coordinates[vertex]
    symbol_scale = G2_TREE_SYMBOL_SCALE * taxon_scale
    dots = [point(vertex) for vertex in keys(coordinates) if vertex > ntaxa]
    segments = [(point(left), point(right), left <= ntaxa || right <= ntaxa ? .6 : 2.)
                for (left, right, _) in edges]
    sizes = [g2_minimal_metric_size(length; font_pt=weight_font_pt) for (_, _, length) in edges]
    glyph = [hypot(size...) / 2 for size in sizes]
    disk = [max((weight_font_pt + 1) / 2, hypot(size...) / 2 + .15) for size in sizes]
    natural = [unit .* g2_minimal_edge_label_position(left, right, coordinates;
                                                      taxon_scale=taxon_scale,
                                                      weight_font_pt=weight_font_pt)
               for (left, right, _) in edges]
    positions = copy(natural)
    # The clearance of the digits, capped at a point and lowered beside the
    # edge, then whether the disk stays on its edge, then the clearance of the
    # disk from the other edges, then the distance moved.
    function score(i, p, beside)
        clear = minimum(g2_minimal_taxon_distance(taxon, p .- point(taxon), symbol_scale)
                        for taxon in 1:ntaxa) - .3 - glyph[i]
        for dot in dots
            clear = min(clear, hypot((p .- dot)...) - dot_radius_pt - glyph[i])
        end
        cover = 0.0
        for j in eachindex(edges)
            j == i && continue
            clear = min(clear, hypot((p .- positions[j])...) -
                               max(disk[j] + glyph[i], disk[i] + glyph[j]))
            a, b, width = segments[j]
            cover = min(cover, g2_segment_distance(p, a, b) - disk[i] - width / 2)
        end
        return (round(min(clear, 1.0) - (beside ? .7 : 0.0); digits=1), !beside,
                round(cover; digits=1),
                -hypot((p .- natural[i])...))
    end
    for _ in 1:3, i in eachindex(edges)
        a, b, _ = segments[i]
        direction = b .- a
        normal = (-direction[2], direction[1]) ./ hypot(direction...)
        candidates = [(natural[i], false)]
        append!(candidates, ((a .+ t .* direction, false) for t in .05:.025:.95))
        beside && append!(candidates,
                          ((a .+ t .* direction .+ side * (disk[i] + .5) .* normal, true)
                           for t in .3:.1:.7 for side in (-1, 1)))
        positions[i] = first(argmax(candidate -> score(i, candidate...), candidates))
    end
    return [position ./ unit for position in positions]
end

function g2_minimal_metric_label(length)
    q = QQ(length)
    # The pendant lengths of a lineality star may be negative; in text mode
    # the sign would be a hyphen.
    denominator(q) == 1 && q < 0 && return "\\(-\\)" * string(-numerator(q))
    denominator(q) == 1 && return string(numerator(q))
    return "\\($(tex_rational(q))\\)"
end

"""
    g2_tree_tikz(tree::TreeMetric; palette=:grayscale, label=false,
                 scale=1.0, taxon_scale=1.0, bounded_length=0.70,
                 edge_labels=:disk)

Return a self-contained TikZ fragment for a labelled G2 phylogenetic tree.
The caller's LaTeX document needs TikZ and xcolor. The `label` option shows
taxon numbers; exact edge lengths are printed unless `edge_labels=:none`.
`scale` changes only
node distances, not strokes, fonts, or shapes. `taxon_scale` multiplies the
sizes of all seven taxon symbols and their numbers relative to the default,
without changing tree coordinates or metric-label fonts. It must be positive
and finite. `bounded_length` changes the
drawn length of each bounded edge before `scale`, without changing the
printed metric. `edge_labels` may be `:disk`, `:aligned` (sloped above
each edge), or `:none` for a purely
combinatorial picture. The edges, the inner vertices, and the edge lengths
are black.
"""
function g2_tree_tikz(tree::TreeMetric;
                      palette::Symbol=:grayscale, label::Bool=false,
                      scale::Real=1.0,
                      taxon_scale::Real=1.0,
                      bounded_length::Real=G2_TREE_BOUNDED_DRAW_LENGTH,
                      edge_labels::Symbol=:disk,
                      layout_coordinates=nothing)
    return sprint() do io
        colors = g2_minimal_tree_styles(io, palette)
        println(io, "\\begin{tikzpicture}[x=$(G2_TREE_UNIT_CM)cm,y=$(G2_TREE_UNIT_CM)cm," *
                    "baseline=(current bounding box.center)]")
        g2_tree_commands(io, tree, colors; label=label, scale=scale,
                         taxon_scale=taxon_scale, bounded_length=bounded_length,
                         edge_labels=edge_labels, layout_coordinates=layout_coordinates)
        println(io, raw"\end{tikzpicture}")
        println(io, raw"\endgroup")
    end
end

"""
    g2_tree_commands(io, tree, colors; label=false, scale=1.0, taxon_scale=1.0,
                     bounded_length=0.70, edge_labels=:disk, layout_coordinates=nothing,
                     vertex_labels=nothing, weight_font_pt=7.5, labels_beside=true)

Write the commands of `g2_tree_tikz` that draw `tree` to `io`, without the
picture and its styles: `colors` is what `g2_minimal_tree_styles` returns,
after it has opened a group and set the styles. The heptagon of the taxa is
centred at the origin. The coordinates are named `v1`, `v2`, ..., and the
commands use them at once, so several trees can share a picture, each in a
shifted scope.

`vertex_labels` maps inner vertices to labels, such as letters, that stand
in white in their dots; the dots are then larger, and the lengths on disks
keep clear of them. `weight_font_pt` is the font size of the lengths on
disks, whose disks the caller enlarges to its size plus a point. With
`labels_beside=false` no length on a disk is moved beside its edge.
"""
function g2_tree_commands(io::IO, tree::TreeMetric, colors;
                          label::Bool=false,
                          scale::Real=1.0,
                          taxon_scale::Real=1.0,
                          bounded_length::Real=G2_TREE_BOUNDED_DRAW_LENGTH,
                          edge_labels::Symbol=:disk,
                          layout_coordinates=nothing,
                          vertex_labels=nothing,
                          weight_font_pt::Real=G2_TREE_WEIGHT_FONT_PT,
                          labels_beside::Bool=true)
    edge_labels in (:disk, :aligned, :none) ||
        throw(ArgumentError("edge_labels must be :disk, :aligned, or :none"))
    taxon_scale_value = Float64(taxon_scale)
    isfinite(taxon_scale_value * G2_TREE_TAXON_FONT_PT) && taxon_scale_value > 0 ||
        throw(ArgumentError("taxon_scale must be a positive finite number"))
    coordinates = g2_minimal_tree_coordinates(tree; scale=scale,
                                               bounded_length=bounded_length,
                                               layout_coordinates=layout_coordinates)
    for vertex in sort!(collect(keys(coordinates)))
        x, y = coordinates[vertex]
        @printf(io, "\\coordinate (v%d) at (%.5f,%.5f);\n", vertex, x, y)
    end
    edges = Tuple{Int,Int,QQFieldElem}[]
    for left in sort!(collect(keys(tree.adjacency)))
        for (right, length) in tree.adjacency[left]
            left < right && push!(edges, (left, right, QQ(length)))
        end
    end
    sort!(edges; by=edge -> (edge[1] <= tree.ntaxa || edge[2] <= tree.ntaxa,
                             edge[1], edge[2]))
    for (left, right, length) in edges
        pendant = left <= tree.ntaxa || right <= tree.ntaxa
        style = pendant ? "g2tree pendant" : "g2tree internal"
        if edge_labels == :aligned
            aligned_style = pendant ? "g2tree aligned pendant" :
                                      "g2tree aligned internal"
            position = pendant ? "0.58" : "0.50"
            if pendant && taxon_scale_value > 1
                x, y = g2_minimal_edge_label_position(left, right, coordinates;
                                                      taxon_scale=taxon_scale_value)
                lx, ly = coordinates[left]
                rx, ry = coordinates[right]
                fraction = hypot(x - lx, y - ly) / hypot(rx - lx, ry - ly)
                position = string(max(0.58, fraction))
            end
            metric = g2_minimal_metric_label(length)
            println(io, "\\draw[$style] (v$left) -- " *
                        "node[$aligned_style,pos=$position] {$metric} (v$right);")
        else
            println(io, "\\draw[$style] (v$left)--(v$right);")
        end
    end
    if edge_labels == :disk
        positions = g2_minimal_edge_label_positions(
            edges, coordinates, tree.ntaxa; taxon_scale=taxon_scale_value,
            dot_radius_pt=vertex_labels === nothing ? G2_TREE_DOT_RADIUS_PT :
                          G2_TREE_LETTERED_DOT_RADIUS_PT,
            weight_font_pt=weight_font_pt, beside=labels_beside,
        )
        font = @sprintf("\\fontsize{%g}{%g}", weight_font_pt, weight_font_pt + .5)
        for ((left, right, length), (x, y)) in zip(edges, positions)
            metric = g2_minimal_metric_label(length)
            @printf(io, "\\node[g2tree weight,font={\\sffamily%s\\selectfont}] at (%.5f,%.5f) {%s};\n",
                    font, x, y, metric)
        end
    end
    vertex_labels === nothing || println(io,
        "\\tikzset{g2tree lettered/.style={minimum size=$(G2_TREE_LETTERED_DOT_SIZE_PT)pt," *
        "text=white,font={\\sffamily\\bfseries\\fontsize{$G2_TREE_LETTER_FONT_PT}" *
        "{$G2_TREE_LETTER_FONT_PT}\\selectfont}}}")
    for vertex in sort!(collect(keys(tree.adjacency)))
        vertex > tree.ntaxa || continue
        text = vertex_labels === nothing ? nothing : get(vertex_labels, vertex, nothing)
        if text === nothing
            println(io, "\\node[g2tree vertex] at (v$vertex) {};")
        else
            println(io, "\\node[g2tree vertex,g2tree lettered] at (v$vertex) {$text};")
        end
    end
    for taxon in 1:tree.ntaxa
        g2_minimal_tree_taxon(io, taxon, colors; label=label,
                              taxon_scale=taxon_scale_value)
    end
    return nothing
end

"""
    g2_wrap_tree_tikz(tikz, description; standalone=false)

The fragment of a picture: its `description` as a comment, then `tikz`. With
`standalone=true`, the document of `g2_standalone_figure`.
"""
function g2_wrap_tree_tikz(tikz::AbstractString, description::AbstractString;
                           standalone::Bool=false)
    fragment = "% " * description * "\n" * tikz
    return standalone ? g2_standalone_figure(fragment) : fragment
end

###############################################################################
# Fragments
###############################################################################

"""The fragment of a selected tree; `options` are the keywords of `g2_tree_tikz`."""
function g2_described_tree_tex(selected; standalone::Bool=false, options...)
    return g2_wrap_tree_tikz(g2_tree_tikz(selected.tree; options...),
                             selected.description; standalone=standalone)
end

"""
    g2_lineality_tree_tex(generator; palette=:grayscale, label=false,
                          scale=1.0, taxon_scale=1.0, bounded_length=0.70,
                          edge_labels=:disk, standalone=false)

Return the input-ready TikZ fragment for the star tree of
`g2_lineality_tree(generator)`, labelled by its pendant lengths. Pass
`standalone=true` for a complete document.
"""
g2_lineality_tree_tex(generator::Integer; options...) =
    g2_described_tree_tex(g2_lineality_tree(generator); options...)

"""
    g2_ray_tree_tex(ray; field=:qq, source=nothing, palette=:grayscale,
                    label=false, scale=1.0, taxon_scale=1.0,
                    bounded_length=0.70, edge_labels=:disk, standalone=false)

Return the input-ready TikZ fragment for the tree of `g2_ray_tree`. Pass
`standalone=true` for a complete document, and a preloaded `source` to render
many rays.
"""
g2_ray_tree_tex(ray::Integer; field::Symbol=:qq, source=nothing, options...) =
    g2_described_tree_tex(g2_ray_tree(ray; field=field, source=source); options...)

"""
    g2_edge_tree_tex(edge; field=:qq, source=nothing, palette=:grayscale,
                     label=false, scale=1.0, taxon_scale=1.0,
                     bounded_length=0.70, edge_labels=:disk, standalone=false)

Return the input-ready TikZ fragment for the tree of `g2_edge_tree`.
"""
g2_edge_tree_tex(edge; field::Symbol=:qq, source=nothing, options...) =
    g2_described_tree_tex(g2_edge_tree(edge; field=field, source=source); options...)

"""
    g2_maximal_cell_tree_tex(cell; field=:qq, source=nothing,
                             palette=:grayscale, label=false, scale=1.0,
                             taxon_scale=1.0, bounded_length=0.70,
                             edge_labels=:disk, standalone=false)

Return a TeX fragment ready for LaTeX's input command, for the tree of
`g2_maximal_cell_tree`. Pass `standalone=true` to wrap it in a complete
compilable document.
"""
g2_maximal_cell_tree_tex(cell; field::Symbol=:qq, source=nothing, options...) =
    g2_described_tree_tex(g2_maximal_cell_tree(cell; field=field, source=source);
                          options...)

"""Write the TeX returned by `tex()` to `path`."""
function g2_write_tree_tex(tex, path::AbstractString)
    isempty(path) && throw(ArgumentError("output path must not be empty"))
    fragment = tex()
    mkpath(dirname(abspath(path)))
    write(path, fragment)
    return path
end

"""Write the input-ready star tree of lineality generator `generator` to `path`."""
write_g2_lineality_tree_tex(path::AbstractString, generator::Integer; kwargs...) =
    g2_write_tree_tex(() -> g2_lineality_tree_tex(generator; kwargs...), path)

"""Write the input-ready tree for ray `ray` to `path`."""
write_g2_ray_tree_tex(path::AbstractString, ray::Integer; kwargs...) =
    g2_write_tree_tex(() -> g2_ray_tree_tex(ray; kwargs...), path)

"""Write the input-ready tree for the edge `edge` to `path`."""
write_g2_edge_tree_tex(path::AbstractString, edge; kwargs...) =
    g2_write_tree_tex(() -> g2_edge_tree_tex(edge; kwargs...), path)

"""
    write_g2_maximal_cell_tree_tex(path, cell; kwargs...)

Write the input-ready TeX fragment to `path`. The same options as
`g2_maximal_cell_tree_tex` are accepted, including `standalone=true`.
"""
write_g2_maximal_cell_tree_tex(path::AbstractString, cell; kwargs...) =
    g2_write_tree_tex(() -> g2_maximal_cell_tree_tex(cell; kwargs...), path)

###############################################################################
# Command line, shared with render_simple_tree.jl
###############################################################################

const G2_TREE_SELECTIONS = ("lineality", "ray", "edge", "cell", "rays")
const G2_TREE_LISTINGS = ("--list-rays", "--list-edges", "--list-cells")
const G2_TREE_SELECTION_USAGE =
    "(--lineality=K | --ray=K | --edge=I,J | --cell=INDEX | --rays=I,J,K[,L]) " *
    "[--field=QQ|GF2]"

"""The value of the option `--name=VALUE` among `args`, or `nothing`."""
function g2_tree_option(args, name)
    values = [split(arg, '='; limit=2)[2] for arg in args
              if startswith(arg, "--$name=")]
    length(values) <= 1 || throw(ArgumentError("pass only one --$name"))
    return isempty(values) ? nothing : String(only(values))
end

"""
Reject every argument that is neither a listing, nor a selection with its
field, nor `--output=FILE`, nor one of the `flags`, nor a `--name=VALUE` with
`name` in `valued`.
"""
function g2_tree_check_arguments(args, flags, valued)
    "--lineality" in args && throw(ArgumentError(
        "choose a generator of the lineality space: " *
        "--lineality=K with K in 1:$(length(G2_TREE_LINEALITY_STARS))",
    ))
    names = ("field", "output", G2_TREE_SELECTIONS..., valued...)
    for arg in args
        arg in G2_TREE_LISTINGS || arg in flags ||
            any(name -> startswith(arg, "--$name="), names) ||
            throw(ArgumentError("unknown option $arg; use --help"))
    end
    g2_tree_option(args, "output") == "" &&
        throw(ArgumentError("output path must not be empty"))
    return nothing
end

"""Print the listing that `args` ask for; returns whether they asked for one."""
function g2_tree_listing(args)
    listed = [arg for arg in args if arg in G2_TREE_LISTINGS]
    isempty(listed) && return false
    length(listed) == 1 &&
        all(arg -> arg in G2_TREE_LISTINGS || startswith(arg, "--field="), args) ||
        throw(ArgumentError("$(first(listed)) takes no other option than --field"))
    field = g2_tree_field(something(g2_tree_option(args, "field"), "QQ"))
    if only(listed) == "--list-rays"
        source = load_g2_ray_source(; field=field)
        println("ray\tray of the $(g2_tree_field_name(field)) census")
        for ray in 1:nrows(source.balanced_rays)
            println(ray, '\t', findfirst(==(ray), source.ray_index))
        end
        return true
    end
    source = load_g2_maximal_cell_source(; field=field)
    if only(listed) == "--list-edges"
        foreach(edge -> println(join(edge, ',')), source.edges)
    else
        for (index, rays) in enumerate(source.cells)
            println(index, '\t', length(rays) == 3 ? "triangle" : "square",
                    '\t', join(rays, ','))
        end
    end
    return true
end

"""The tree that `args` select, with its description."""
function g2_tree_selected(args)
    chosen = [name for name in G2_TREE_SELECTIONS
              if g2_tree_option(args, name) !== nothing]
    length(chosen) == 1 || throw(ArgumentError(
        "specify exactly one of --lineality=K, --ray=K, --edge=I,J, --cell=INDEX, " *
        "or --rays=I,J,K[,L]",
    ))
    field = g2_tree_field(something(g2_tree_option(args, "field"), "QQ"))
    value = g2_tree_option(args, only(chosen))
    indices() = parse.(Int, split(value, ','))
    only(chosen) == "lineality" && return g2_lineality_tree(parse(Int, value))
    only(chosen) == "ray" && return g2_ray_tree(parse(Int, value); field=field)
    only(chosen) == "edge" && return g2_edge_tree(indices(); field=field)
    only(chosen) == "cell" && return g2_maximal_cell_tree(parse(Int, value); field=field)
    return g2_maximal_cell_tree(indices(); field=field)
end

"""Print the fragment `tex`, or write it to the `--output=FILE` of `args`."""
function g2_tree_emit(args, tex::AbstractString)
    path = g2_tree_option(args, "output")
    if path === nothing
        print(tex)
    else
        g2_write_tree_tex(() -> tex, path)
        println("Wrote ", path)
    end
    return nothing
end

function g2_maximal_cell_tree_main(args=ARGS)
    script = "julia --project=. to_tex/render_maximal_cell_tree.jl"
    if "--help" in args
        println("Usage: $script $G2_TREE_SELECTION_USAGE " *
                "[--color | --grayscale] [--label=true|false] " *
                "[--scale=FACTOR] [--taxon-scale=FACTOR] [--bounded-length=LENGTH] " *
                "[--edge-labels=disk|aligned|none] " *
                "[--standalone] [--output=FILE]")
        println("       $script (--list-rays | --list-edges | --list-cells) " *
                "[--field=QQ|GF2]")
        return nothing
    end
    g2_tree_check_arguments(
        args,
        ("--color", "--grayscale", "--label=true", "--label=false", "--standalone"),
        ("scale", "taxon-scale", "bounded-length", "edge-labels"),
    )
    g2_tree_listing(args) && return nothing
    !("--color" in args && "--grayscale" in args) ||
        throw(ArgumentError("choose either --color or --grayscale"))
    !("--label=true" in args && "--label=false" in args) ||
        throw(ArgumentError("choose one label setting"))
    number(name, default) = parse(Float64, something(g2_tree_option(args, name), default))
    tex = g2_described_tree_tex(
        g2_tree_selected(args);
        palette="--color" in args ? :color : :grayscale,
        label="--label=true" in args,
        scale=number("scale", "1.0"),
        taxon_scale=number("taxon-scale", "1.0"),
        bounded_length=number("bounded-length", string(G2_TREE_BOUNDED_DRAW_LENGTH)),
        edge_labels=Symbol(something(g2_tree_option(args, "edge-labels"), "disk")),
        standalone="--standalone" in args,
    )
    return g2_tree_emit(args, tex)
end

if abspath(PROGRAM_FILE) == @__FILE__
    g2_maximal_cell_tree_main()
end
