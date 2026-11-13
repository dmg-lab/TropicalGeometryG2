include(joinpath(@__DIR__, "find_better_gens.jl"))

using LinearAlgebra
using Printf

const DEFAULT_BALANCED_FAN_PATH =
    TropicalGeometryG2.data_path("g2", "reference", "fine_groebner", "normalized_fine_fan.mrdi")
const DEFAULT_BALANCED_RAYS_PATH =
    TropicalGeometryG2.data_path("g2", "reference", "fine_groebner", "normalized_fine_rays.mrdi")
const DEFAULT_ATLAS_TEX_PATH = joinpath(@__DIR__, "balanced_phylogenetic_trees.tex")

"""A full-rank, human-readable basis for the common lineality space."""
function balanced_lineality_basis()
    star_rows = star_lineality_matrix()
    M = Matrix{QQFieldElem}(undef, 3, length(PLUECKER_PAIRS))
    M[1, :] = star_rows[1, :] ./ 2  # the primitive all-ones Pluecker direction
    M[2, :] = star_rows[2, :]
    M[3, :] = star_rows[3, :]
    @assert star_rows[4, :] == -star_rows[2, :] - star_rows[3, :]
    @assert rank(matrix(QQ, M)) == 3
    return M
end

function vectors_as_matrix(rows::AbstractVector{<:AbstractVector})
    M = Matrix{QQFieldElem}(undef, length(rows), length(first(rows)))
    for (i, row) in enumerate(rows)
        M[i, :] = QQ.(row)
    end
    return M
end

"""
    rebuild_balanced_fan(fan, records)

Re-encode `fan` with the balanced ray representatives while retaining exactly
the original maximal-cone/ray incidence.  Replacing a ray `r_i` by `r_i+l_i`
with `l_i` in the common lineality space does not change any cone.
"""
function rebuild_balanced_fan(fan, records)
    incidence = maximal_cones(IncidenceMatrix, fan)
    balanced_rays = vectors_as_matrix([record.balanced_generator for record in records])
    lineality = balanced_lineality_basis()
    rebuilt = polyhedral_fan(
        incidence,
        matrix(QQ, balanced_rays),
        matrix(QQ, lineality);
        non_redundant=true,
    )

    @assert ambient_dim(rebuilt) == ambient_dim(fan)
    @assert dim(rebuilt) == dim(fan)
    @assert lineality_dim(rebuilt) == lineality_dim(fan) == 3
    @assert length(rays_modulo_lineality(rebuilt).rays_modulo_lineality) == length(records)
    @assert Matrix{Bool}(maximal_cones(IncidenceMatrix, rebuilt)) == Matrix{Bool}(incidence)
    return rebuilt
end

"""Ray-index orbits of the combinatorial automorphism group of the fan."""
function ray_orbits(fan)
    complex = simplicial_complex(maximal_cones(IncidenceMatrix, fan))
    group = automorphism_group(complex)
    result = [sort!(collect(orbit)) for orbit in orbits(group)]
    sort!(result; by=first)
    return result
end

"""Whether the internal graph of `tree` is a path (resolved or not)."""
function is_path_shaped_tree(tree::TreeMetric)
    internal = collect((tree.ntaxa + 1):maximum(keys(tree.adjacency)))
    return all(count(neighbor > tree.ntaxa for (neighbor, _) in tree.adjacency[vertex]) <= 2
               for vertex in internal)
end

"""
Return the internal vertices of a path-shaped tree from one endpoint to the
other. Internal node labels are used only to choose an initial orientation;
the final leaf order is made orientation-independent below.
"""
function internal_spine(tree::TreeMetric)
    internal = collect((tree.ntaxa + 1):maximum(keys(tree.adjacency)))
    length(internal) == 1 && return internal
    endpoints = sort!([
        vertex
        for vertex in internal
        if count(neighbor > tree.ntaxa for (neighbor, _) in tree.adjacency[vertex]) == 1
    ])
    @assert length(endpoints) == 2

    result = [first(endpoints)]
    previous = 0
    while last(result) != last(endpoints)
        next_vertex = only(
            neighbor
            for (neighbor, _) in tree.adjacency[last(result)]
            if neighbor > tree.ntaxa && neighbor != previous
        )
        previous = last(result)
        push!(result, next_vertex)
    end
    @assert sort(result) == sort(internal)
    return result
end

"""
Canonical convex-arc leaf order for a path-shaped tree. Starting at taxon 7,
visit the spine groups near-to-far on one side of its attachment and
far-to-near on the other. The reverse cyclic order is the only mirror choice;
choose the lexicographically smaller one. Harmonic internal vertices then form
a monotone convex arc, including when taxon 7 attaches in the spine interior.
"""
function path_arc_leaf_order_from_seven(tree::TreeMetric)
    special_taxon = 7
    spine = internal_spine(tree)
    leaf_groups = [
        sort!([neighbor for (neighbor, _) in tree.adjacency[vertex]
               if neighbor <= tree.ntaxa])
        for vertex in spine
    ]
    special_position = only(i for i in eachindex(leaf_groups)
                            if special_taxon in leaf_groups[i])

    left_near_to_far = reduce(
        vcat,
        (leaf_groups[i] for i in (special_position - 1):-1:1);
        init=Int[],
    )
    right_far_to_near = reduce(
        vcat,
        (leaf_groups[i] for i in length(spine):-1:(special_position + 1));
        init=Int[],
    )
    special_companions = filter(!=(special_taxon), leaf_groups[special_position])
    order = vcat([special_taxon], left_near_to_far, right_far_to_near, special_companions)
    reflected_order = vcat([special_taxon], reverse(order[2:end]))
    Tuple(reflected_order) < Tuple(order) && (order = reflected_order)

    @assert sort(order) == collect(1:tree.ntaxa)
    return order
end

"""Canonical taxon-7-rooted planar order for a branched (snowflake) tree."""
function snowflake_leaf_order_from_seven(tree::TreeMetric)
    special_taxon = 7
    root = only(neighbor for (neighbor, _) in tree.adjacency[special_taxon]
                if neighbor > tree.ntaxa)
    function visit(vertex, parent)
        vertex <= tree.ntaxa && return [vertex]
        branches = [
            (neighbor=neighbor, leaves=visit(neighbor, vertex))
            for (neighbor, _) in tree.adjacency[vertex]
            if neighbor != parent
        ]
        # Keep each arm contiguous. Continue through internal nodes before
        # pendant leaves, then use taxa only to resolve mirror symmetries.
        sort!(branches; by=branch -> (
            branch.neighbor <= tree.ntaxa,
            minimum(branch.leaves),
            join(branch.leaves, ","),
        ))
        return reduce(vcat, (branch.leaves for branch in branches); init=Int[])
    end

    order = vcat([special_taxon], visit(root, special_taxon))
    reflected_order = vcat([special_taxon], reverse(order[2:end]))
    Tuple(reflected_order) < Tuple(order) && (order = reflected_order)
    @assert sort(order) == collect(1:tree.ntaxa)
    return order
end

function topology_compatible_leaf_order(tree::TreeMetric)
    return is_path_shaped_tree(tree) ?
           path_arc_leaf_order_from_seven(tree) :
           snowflake_leaf_order_from_seven(tree)
end

"""
    symmetric_layout(tree; radius=2.55)

Place leaves at equal angles on a regular polygon in a topology-compatible
cyclic order.  Internal vertices receive their harmonic (barycentric)
positions.  This gives a symmetric, straight-line unrooted tree drawing.
"""
function symmetric_layout(tree::TreeMetric; radius=2.55)
    order = topology_compatible_leaf_order(tree)
    coordinates = Dict{Int, Tuple{Float64, Float64}}()
    for (position, taxon) in enumerate(order)
        angle = pi / 2 - 2pi * (position - 1) / tree.ntaxa
        coordinates[taxon] = (radius * cos(angle), radius * sin(angle))
    end

    internal = collect((tree.ntaxa + 1):maximum(keys(tree.adjacency)))
    index = Dict(vertex => i for (i, vertex) in enumerate(internal))
    A = zeros(Float64, length(internal), length(internal))
    bx = zeros(Float64, length(internal))
    by = zeros(Float64, length(internal))
    for vertex in internal
        i = index[vertex]
        A[i, i] = length(tree.adjacency[vertex])
        for (neighbor, _) in tree.adjacency[vertex]
            if neighbor <= tree.ntaxa
                x, y = coordinates[neighbor]
                bx[i] += x
                by[i] += y
            else
                A[i, index[neighbor]] -= 1
            end
        end
    end
    x_internal, y_internal = A \ bx, A \ by
    for vertex in internal
        coordinates[vertex] = (x_internal[index[vertex]], y_internal[index[vertex]])
    end
    return coordinates
end

function tex_rational(x)
    q = QQ(x)
    p, d = numerator(q), denominator(q)
    d == 1 && return string(p)
    return p < 0 ? "-\\frac{$(abs(p))}{$d}" : "\\frac{$p}{$d}"
end

tex_tuple(values) = "(" * join(tex_rational.(values), ",\\,") * ")"

function write_tree_tikz(io::IO, record; radius=2.55)
    tree = record.tree
    coordinates = symmetric_layout(tree; radius=radius)
    for vertex in sort!(collect(keys(coordinates)))
        x, y = coordinates[vertex]
        println(io, @sprintf("\\coordinate (v%d) at (%.5f,%.5f);", vertex, x, y))
    end

    edges = Tuple{Int, Int, QQFieldElem}[]
    for u in sort!(collect(keys(tree.adjacency)))
        for (v, edge_length) in tree.adjacency[u]
            u < v && push!(edges, (u, v, edge_length))
        end
    end
    # Draw internal edges first, then the radial pendant edges.
    sort!(edges; by=edge -> (edge[1] <= tree.ntaxa || edge[2] <= tree.ntaxa, edge[1], edge[2]))
    for (u, v, edge_length) in edges
        pendant = u <= tree.ntaxa || v <= tree.ntaxa
        style = pendant ? "pendant edge" : "internal edge"
        # Pendant edges are stored leaf-first (taxon labels are smaller than
        # internal-node labels). Put their lengths just beyond the midpoint,
        # toward the much smaller internal vertex: the old outward position
        # could place a short-edge label underneath the taxon circle.
        label_position = pendant ? "0.55" : "0.50"
        println(
            io,
            "\\draw[$style] (v$u) -- node[edge label,pos=$label_position] {\$$(tex_rational(edge_length))\$} (v$v);",
        )
    end

    for vertex in (tree.ntaxa + 1):maximum(keys(tree.adjacency))
        println(io, "\\node[internal vertex] at (v$vertex) {};")
    end
    for taxon in 1:tree.ntaxa
        style = taxon == 7 ? "special leaf" : "leaf"
        println(io, "\\node[$style] at (v$taxon) {$taxon};")
    end
end

function write_tree_card(io::IO, ray_index::Int, record)
    println(io, raw"\begin{minipage}[t][0.285\textheight][t]{0.485\textwidth}")
    println(io, raw"\centering")
    println(io, "{\\large\\bfseries Ray $ray_index}\\par")
    println(io, "{\\scriptsize\\color{muted} Coxeter shift \$c=$(tex_tuple(record.root_coefficients))\$}\\par")
    println(
        io,
        "{\\scriptsize\\color{muted} integral scale \$m=$(record.integer_scale)\$, " *
        "diagonal shift \$a=$(tex_rational(record.diagonal_leaf_shift))\$, " *
        "gcd \$g=$(record.primitive_divisor)\$}\\par\\vspace{1mm}",
    )
    println(io, raw"\begin{tikzpicture}[x=1cm,y=1cm]")
    write_tree_tikz(io, record)
    println(io, raw"\end{tikzpicture}")
    println(io, raw"\end{minipage}")
end


function write_atlas_preamble(io::IO)
    print(
        io,
        raw"""\documentclass[10pt]{article}
\usepackage[a4paper,left=11mm,right=11mm,top=18mm,bottom=11mm,
            includehead,headheight=14pt,headsep=4mm]{geometry}
\usepackage{amsmath}
\usepackage{microtype}
\usepackage{xcolor}
\usepackage{tikz}
\usepackage{fancyhdr}
\definecolor{ink}{HTML}{17212B}
\definecolor{muted}{HTML}{637083}
\definecolor{internalblue}{HTML}{245B78}
\definecolor{pendantgray}{HTML}{8795A5}
\definecolor{paperblue}{HTML}{EDF4F7}
\definecolor{specialseven}{HTML}{C4583C}
\tikzset{
  pendant edge/.style={draw=pendantgray,line width=0.8pt},
  internal edge/.style={draw=internalblue,line width=2.0pt},
  edge label/.style={font=\scriptsize,fill=white,fill opacity=.88,text opacity=1,
                     inner sep=1.1pt,sloped,above},
  internal vertex/.style={circle,fill=internalblue,draw=white,line width=.4pt,
                          minimum size=4.2pt,inner sep=0pt},
  leaf/.style={circle,fill=white,draw=ink,line width=.75pt,minimum size=5.2mm,
               inner sep=0pt,font=\small\bfseries,text=ink},
  special leaf/.style={leaf,fill=specialseven,draw=specialseven,text=white},
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small\color{muted} Primitive integral Coxeter-balanced $G_2$ trees}
\fancyhead[R]{\small\color{muted}\thepage}
\renewcommand{\headrulewidth}{0pt}
\setlength{\parindent}{0pt}
\color{ink}
\begin{document}
\begin{center}
  {\Huge\bfseries Primitive integral $G_2$ trees}\par\vspace{3mm}
  {\large A symmetric atlas of the 51 rays of $\operatorname{Trop}(G_2/ P)$}\par\vspace{7mm}
  \begin{tikzpicture}
    \node[rounded corners=3pt,fill=paperblue,inner xsep=10mm,inner ysep=5mm,
          text width=.78\textwidth,align=center] {
      For pendant lengths $\ell$ and Coxeter-root stars $r_1,r_2,r_3$, each ray is
      represented by
      \[
        c_i=\frac{\langle\ell,r_i\rangle}{6},\qquad
        \ell_{\mathrm{bal}}=\ell-\sum_{i=1}^3c_i r_i,
        \qquad
        \ell_{0}=\frac{m\ell_{\mathrm{bal}}-a\mathbf 1}{g},
\]
      where $m$ is the smallest common-denominator scale and
      $a=\min_i(m\ell_{\mathrm{bal},i})$ is the diagonal-torus shift; $g$ is
      the gcd of all shifted tree-edge lengths. The two star-tree operations
      are lineality translations, performed before the positive scaling by
      $1/g$. Thus $\ell_0$ is primitive integral, its shortest pendant edge is
      zero, and Coxeter balance is preserved. Thin grey edges are pendant;
      thick blue edges are the lineality-invariant internal splits.
    };
  \end{tikzpicture}
\end{center}
\vfill
\begin{center}\color{muted}
Taxon $7$ is the rust-colored top vertex of every regular heptagon. Internal
vertices use a harmonic layout, zero-length internal refinements are contracted,
and every label is integral.
\end{center}
""",
    )
end

function write_orbit_header(io::IO, orbit_number, orbit, continuation)
    println(io, raw"\clearpage")
    suffix = continuation ? " (continued)" : ""
    println(io, raw"\begin{tikzpicture}")
    println(
        io,
        "\\node[rounded corners=2pt,fill=paperblue,inner xsep=5mm,inner ysep=2.2mm,text width=.94\\textwidth] " *
        "{\\Large\\bfseries Orbit $orbit_number$suffix\\hfill\\normalsize $(length(orbit)) rays};",
    )
    println(io, raw"\end{tikzpicture}\par\vspace{1mm}")
    println(io, "{\\small\\color{muted} Ray indices: $(join(orbit, ", "))}\\par\\vspace{2mm}")
end

function write_phylogenetic_atlas(io::IO, records, orbit_list)
    write_atlas_preamble(io)
    per_page = 6
    for (orbit_number, orbit) in enumerate(orbit_list)
        for page_start in 1:per_page:length(orbit)
            write_orbit_header(io, orbit_number, orbit, page_start > 1)
            page_rays = orbit[page_start:min(page_start + per_page - 1, length(orbit))]
            for (position, ray_index) in enumerate(page_rays)
                lone_last_card = isodd(length(page_rays)) && position == length(page_rays)
                lone_last_card && println(io, raw"\hfill")
                write_tree_card(io, ray_index, records[ray_index])
                if lone_last_card
                    println(io, raw"\hfill\par\vspace{1mm}")
                elseif isodd(position)
                    println(io, raw"\hfill")
                else
                    println(io, raw"\par\vspace{1mm}")
                end
            end
        end
    end
    println(io, raw"\end{document}")
end

function compile_atlas(tex_path::AbstractString)
    latexmk = Sys.which("latexmk")
    latexmk === nothing && error("latexmk is required to compile the atlas PDF")
    pdf_name = replace(basename(tex_path), r"\.tex$" => ".pdf")
    pdf_path = joinpath(dirname(tex_path), pdf_name)
    mktempdir() do build_dir
        temporary_tex = joinpath(build_dir, basename(tex_path))
        cp(tex_path, temporary_tex; force=true)
        command = `$(latexmk) -silent -pdf -interaction=nonstopmode -halt-on-error -file-line-error $(basename(tex_path))`
        run(Cmd(command; dir=build_dir))
        cp(joinpath(build_dir, pdf_name), pdf_path; force=true)
    end
    return pdf_path
end

function build_balanced_fan_and_atlas(
    ;
    fan_path::AbstractString=DEFAULT_FAN_PATH,
    balanced_fan_path::AbstractString=DEFAULT_BALANCED_FAN_PATH,
    balanced_rays_path::AbstractString=DEFAULT_BALANCED_RAYS_PATH,
    tex_path::AbstractString=DEFAULT_ATLAS_TEX_PATH,
    compile_pdf::Bool=true,
)
    fan = load(fan_path)
    verify_star_lineality(fan) || error("the stated star rows do not span the fan lineality")
    records = balanced_ray_tree_data(fan)
    rebuilt = rebuild_balanced_fan(fan, records)
    mkpath(dirname(balanced_fan_path))
    save(balanced_fan_path, rebuilt)

    # Polymake regards ray generators up to positive scaling and normalizes
    # RAYS when a fan is constructed. Preserve the chosen primitive-integral
    # tree normalization in a companion exact matrix.
    balanced_ray_matrix = matrix(
        QQ,
        vectors_as_matrix([record.balanced_generator for record in records]),
    )
    mkpath(dirname(balanced_rays_path))
    save(balanced_rays_path, balanced_ray_matrix)
    @assert load(balanced_rays_path) == balanced_ray_matrix

    orbit_list = ray_orbits(fan)
    @assert sort!(reduce(vcat, orbit_list)) == collect(1:length(records))
    open(tex_path, "w") do io
        write_phylogenetic_atlas(io, records, orbit_list)
    end

    pdf_path = compile_pdf ? compile_atlas(tex_path) : nothing
    println("balanced fan: ", balanced_fan_path)
    println("primitive integral ray matrix: ", balanced_rays_path)
    println("fan data: ambient=$(ambient_dim(rebuilt)), dim=$(dim(rebuilt)), lineality=$(lineality_dim(rebuilt)), " *
            "rays=$(length(records)), maximal cones=$(length(maximal_cones(rebuilt)))")
    println("ray orbits: ", length.(orbit_list))
    println("LaTeX atlas: ", tex_path)
    pdf_path === nothing || println("PDF atlas: ", pdf_path)
    return (
        fan=rebuilt,
        ray_matrix=balanced_ray_matrix,
        records=records,
        orbits=orbit_list,
        tex_path=tex_path,
        pdf_path=pdf_path,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    build_balanced_fan_and_atlas(; compile_pdf=!("--no-pdf" in ARGS))
end
