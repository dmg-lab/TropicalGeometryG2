include(joinpath(@__DIR__, "build_balanced_fan_and_atlas.jl"))

const DEFAULT_BALANCED_TRIANGLES_PATH =
    TropicalGeometryG2.data_path("g2", "reference", "fine_groebner", "canonical_fine_weights.mrdi")
const DEFAULT_TRIANGLE_INDICES_PATH =
    TropicalGeometryG2.data_path("g2", "reference", "fine_groebner", "fine_cell_rays.mrdi")
const DEFAULT_TRIANGLE_ATLAS_TEX_PATH =
    joinpath(@__DIR__, "balanced_phylogenetic_triangles.tex")
const DEFAULT_TRIANGLE_ANALYSIS_CSV_PATH =
    TropicalGeometryG2.data_path(
        "g2", "reference", "fine_groebner", "noncanonical_fine_initial_ideals.csv",
    )

function parse_csv_bool(value::AbstractString)
    value == "true" && return true
    value == "false" && return false
    throw(ArgumentError("expected CSV boolean, received '$value'"))
end

function parse_csv_ints(value::AbstractString)
    isempty(value) && return Int[]
    return parse.(Int, split(value, ';'))
end

"""Read the exact per-triangle analysis, `noncanonical_fine_initial_ideals.csv`."""
function load_triangle_analysis_csv(path::AbstractString)
    isfile(path) || error("triangle analysis CSV not found at $path")
    lines = readlines(path)
    isempty(lines) && error("triangle analysis CSV is empty: $path")
    expected_header = [
        "triangle", "orbit", "rays", "m2_weight", "prime_over_QQ", "binomial",
        "monomial_free_tropical", "torus_saturated", "toric", "dimension", "degree",
        "groebner_generators", "maximum_terms", "relation_lattice_certified",
        "body_f_vector", "body_vertices", "body_facets", "body_degree_one_values",
        "body_lattice_points", "body_normalized_volume", "chevalley_simplex",
    ]
    header = split(chomp(first(lines)), ',')
    header == expected_header || error("unexpected triangle analysis CSV header")

    records = Any[]
    for (offset, line) in enumerate(lines[2:end])
        line_number = offset + 1
        fields = split(chomp(line), ','; keepempty=true)
        length(fields) == length(expected_header) ||
            error("CSV line $line_number has $(length(fields)) fields, expected $(length(expected_header))")
        optional_int(field) = isempty(field) ? nothing : parse(Int, field)
        triangle = parse(Int, fields[1])
        push!(records, (
            triangle=triangle,
            orbit=parse(Int, fields[2]),
            rays=parse_csv_ints(fields[3]),
            m2_weight=parse_csv_ints(fields[4]),
            prime=parse_csv_bool(fields[5]),
            binomial=parse_csv_bool(fields[6]),
            monomial_free=parse_csv_bool(fields[7]),
            torus_saturated=parse_csv_bool(fields[8]),
            toric=parse_csv_bool(fields[9]),
            dimension=parse(Int, fields[10]),
            degree=parse(Int, fields[11]),
            groebner_generators=parse(Int, fields[12]),
            maximum_terms=parse(Int, fields[13]),
            relation_lattice_certified=parse_csv_bool(fields[14]),
            body_f_vector=parse_csv_ints(fields[15]),
            body_vertices=optional_int(fields[16]),
            body_facets=optional_int(fields[17]),
            body_degree_one_values=optional_int(fields[18]),
            body_lattice_points=optional_int(fields[19]),
            body_normalized_volume=optional_int(fields[20]),
            chevalley_simplex=isempty(fields[21]) ? nothing : parse_csv_bool(fields[21]),
        ))
    end
    @assert length(records) == 251
    @assert [record.triangle for record in records] == collect(1:251)
    return records
end

"""The maximal triangles, in Oscar's facet order, as sorted ray-index triples."""
function maximal_triangle_ray_indices(fan)
    incidence = maximal_cones(IncidenceMatrix, fan)
    complex = simplicial_complex(incidence)
    triangles = [sort!(collect(facet)) for facet in facets(complex)]

    @assert length(triangles) == size(incidence, 1)
    @assert all(length(triangle) == 3 for triangle in triangles)

    incidence_triangles = Set(
        Tuple(findall(identity, collect(row)))
        for row in eachrow(Matrix{Bool}(incidence))
    )
    @assert Set(Tuple(triangle) for triangle in triangles) == incidence_triangles
    return triangles
end

"""Facet-index orbits, using the same facet order as `maximal_triangle_ray_indices`."""
function triangle_orbits(fan)
    complex = simplicial_complex(maximal_cones(IncidenceMatrix, fan))
    group = automorphism_group(complex, action=:on_facets)
    result = [sort!(collect(orbit)) for orbit in orbits(group)]
    sort!(result; by=first)
    return result
end

function row_vector(M, row::Int)
    return [QQ(M[row, column]) for column in 1:ncols(M)]
end

"""
    balanced_triangle_tree_data(fan, balanced_ray_matrix)

For every maximal triangle `tau={i,j,k}`, form the relative-interior point
`w_tau=w_i+w_j+w_k` from the chosen primitive-integral ray representatives.
Its tree is already integral and Coxeter-balanced. Apply the diagonal star
direction so that its shortest pendant edge has length zero, then divide all
tree lengths by their gcd.
"""
function balanced_triangle_tree_data(fan, balanced_ray_matrix)
    triangles = maximal_triangle_ray_indices(fan)
    @assert nrows(balanced_ray_matrix) == 51
    @assert ncols(balanced_ray_matrix) == length(PLUECKER_PAIRS) == 21

    roots = [rational_row(COXETER_STAR_ROOTS, i) for i in axes(COXETER_STAR_ROOTS, 1)]
    records = Any[]
    for (triangle_index, triangle) in enumerate(triangles)
        ray_sum = zeros(QQFieldElem, ncols(balanced_ray_matrix))
        for ray_index in triangle
            ray_sum .+= row_vector(balanced_ray_matrix, ray_index)
        end

        record = balance_outer_edges(ray_sum)
        @assert all(iszero, record.root_coefficients)
        @assert all(iszero, record.lineality_leaf_shift)
        @assert record.integer_scale == 1
        @assert minimum(record.pendant_lengths) == 0
        @assert all(iszero(inner_product(record.pendant_lengths, root)) for root in roots)
        @assert all(edge_length >= 0 && denominator(edge_length) == 1
                    for edge_length in tree_edge_lengths(record.tree))

        push!(records, merge(
            (triangle_index=triangle_index, rays=triangle, ray_sum=ray_sum),
            record,
        ))
    end
    return records
end

function triangle_status_tex(analysis)
    if analysis.toric
        return raw"\textcolor{toricgreen}{\bfseries PRIME TORIC}"
    elseif analysis.prime
        return raw"\textcolor{warningamber}{\bfseries PRIME NONBINOMIAL}"
    end
    return raw"\textcolor{nonprimered}{\bfseries NON-PRIME}"
end

function body_summary_tex(analysis)
    analysis.toric || return nothing
    f_vector = tex_tuple(analysis.body_f_vector)
    return "\$f=$f_vector\$\\quad " *
           "degree-one/lattice \$=$(analysis.body_degree_one_values)/$(analysis.body_lattice_points)\$"
end

function write_triangle_tree_card(io::IO, triangle_index::Int, record, analysis)
    @assert analysis.triangle == triangle_index
    @assert analysis.rays == record.rays
    @assert analysis.m2_weight == [-Int(QQ(x)) for x in record.balanced_generator]
    rays = join(record.rays, ",\\,")
    println(io, raw"\begin{minipage}[t][0.225\textheight][t]{0.485\textwidth}")
    println(io, raw"\centering")
    println(io, "{\\normalsize\\bfseries Triangle $triangle_index}\\par")
    println(io, "{\\small\\color{muted} Rays \$\\{$rays\\}\$}\\par\\vspace{0.3mm}")
    println(
        io,
        "{\\scriptsize $(triangle_status_tex(analysis))\\color{muted}\\quad " *
        "\$\\dim/\\deg=$(analysis.dimension)/$(analysis.degree)\$\\quad " *
        "GB $(analysis.groebner_generators)/$(analysis.maximum_terms)}\\par",
    )
    body_summary = body_summary_tex(analysis)
    if body_summary === nothing
        println(io, raw"\vspace{1.9mm}")
    else
        println(io, "{\\scriptsize\\color{bodyviolet} $body_summary}\\par\\vspace{0.3mm}")
    end
    println(io, raw"\begin{tikzpicture}[x=1cm,y=1cm]")
    write_tree_tikz(io, record; radius=1.95)
    println(io, raw"\end{tikzpicture}")
    println(io, raw"\end{minipage}")
end

function write_triangle_atlas_preamble(io::IO, orbit_sizes, analyses)
    print(
        io,
        raw"""\documentclass[10pt]{article}
\usepackage[a4paper,left=11mm,right=11mm,top=17mm,bottom=5mm,
            includehead,headheight=14pt,headsep=3mm]{geometry}
\usepackage{amsmath}
\usepackage{amssymb}
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
\definecolor{toricgreen}{HTML}{167A5B}
\definecolor{warningamber}{HTML}{A96513}
\definecolor{nonprimered}{HTML}{B33A3A}
\definecolor{bodyviolet}{HTML}{654A87}
\tikzset{
  pendant edge/.style={draw=pendantgray,line width=0.75pt},
  internal edge/.style={draw=internalblue,line width=1.8pt},
  edge label/.style={font=\tiny,fill=white,fill opacity=.88,text opacity=1,
                     inner sep=.8pt,sloped,above},
  internal vertex/.style={circle,fill=internalblue,draw=white,line width=.35pt,
                          minimum size=3.8pt,inner sep=0pt},
  leaf/.style={circle,fill=white,draw=ink,line width=.7pt,minimum size=4.4mm,
               inner sep=0pt,font=\scriptsize\bfseries,text=ink},
  special leaf/.style={leaf,fill=specialseven,draw=specialseven,text=white},
}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small\color{muted} Primitive integral trees of the maximal $G_2$ triangles}
\fancyhead[R]{\small\color{muted}\thepage}
\renewcommand{\headrulewidth}{0pt}
\setlength{\parindent}{0pt}
\color{ink}
\begin{document}
\begin{center}
  {\Huge\bfseries Trees of the maximal triangles}\par\vspace{3mm}
  {\large A symmetric atlas of the 251 maximal cones of $\operatorname{Trop}(G_2/P)$}\par\vspace{7mm}
  \begin{tikzpicture}
    \node[rounded corners=3pt,fill=paperblue,inner xsep=10mm,inner ysep=5mm,
          text width=.80\textwidth,align=center] {
      For a maximal triangle $\tau=\{i,j,k\}$, first take the relative-interior
      representative and its pendant lengths
      \[
        w_\tau=w_i+w_j+w_k,\qquad
        \ell_\tau=\ell_i+\ell_j+\ell_k.
      \]
      The chosen ray representatives make $\ell_\tau$ integral and orthogonal
      to all three Coxeter roots. The only remaining normalization is the
      diagonal-torus translation followed by primitive lattice rescaling
      \[
        a_\tau=\min_s\ell_{\tau,s},\qquad
        \ell_\tau^0=\frac{\ell_\tau-a_\tau\mathbf 1}{g_\tau}.
      \]
      Subtracting $a_\tau\mathbf 1$ is a lineality translation; dividing by
      the edge-length gcd $g_\tau$ is a positive scaling, performed afterwards.
      Thus every displayed edge length is a nonnegative integer, the complete
      edge vector is primitive, and the shortest pendant edge is zero. Thin
      grey edges are pendant; thick blue edges are the lineality-invariant
      internal splits. Each card also records a concise initial-ideal
      classification at the displayed triangle's relative-interior point.
    };
  \end{tikzpicture}
\end{center}
\vfill
\begin{center}\color{muted}
The 251 triangles form 24 combinatorial facet orbits. Taxon $7$ is the
rust-colored top vertex of each regular heptagon; internal vertices use a
harmonic layout.\\[1mm]
""",
    )
    println(io, "Orbit sizes: \\($(join(orbit_sizes, ",\\,"))\\).")
    toric = count(analysis -> analysis.toric, analyses)
    prime_nonbinomial = count(analysis -> analysis.prime && !analysis.binomial, analyses)
    nonprime = count(analysis -> !analysis.prime, analyses)
    println(
        io,
        "\\\\[1mm]{\\textcolor{toricgreen}{\\bfseries $toric prime toric}}\\quad" *
        "{\\textcolor{warningamber}{\\bfseries $prime_nonbinomial prime non-binomial}}\\quad" *
        "{\\textcolor{nonprimered}{\\bfseries $nonprime non-prime}}.",
    )
    println(io, "\\\\[1mm]{\\scriptsize On each card, `GB \\(a/b\\)' means \\(a\\) reduced Gr\\\"obner-basis elements and at most \\(b\\) terms per element; for a Newton--Okounkov body, `degree-one/lattice' compares distinct Pl\\\"ucker values with all lattice points.}")
    println(io, raw"\end{center}")
end

function write_triangle_orbit_header(io::IO, orbit_number, orbit, continuation, analyses)
    println(io, raw"\clearpage")
    suffix = continuation ? " (continued)" : ""
    println(io, raw"\begin{tikzpicture}")
    println(
        io,
        "\\node[rounded corners=2pt,fill=paperblue,inner xsep=5mm,inner ysep=2mm,text width=.94\\textwidth] " *
        "{\\Large\\bfseries Orbit $orbit_number$suffix\\hfill\\normalsize $(length(orbit)) triangles};",
    )
    println(io, raw"\end{tikzpicture}\par\vspace{0.5mm}")
    orbit_analyses = analyses[orbit]
    toric = count(analysis -> analysis.toric, orbit_analyses)
    prime_nonbinomial = count(analysis -> analysis.prime && !analysis.binomial, orbit_analyses)
    nonprime = count(analysis -> !analysis.prime, orbit_analyses)
    println(
        io,
        "{\\scriptsize\\color{muted} Triangle indices: $(join(orbit, ", "))}\\par" *
        "{\\scriptsize\\color{muted} Algebra: \\textcolor{toricgreen}{$toric toric}, " *
        "\\textcolor{warningamber}{$prime_nonbinomial prime non-binomial}, " *
        "\\textcolor{nonprimered}{$nonprime non-prime}.}\\par\\vspace{1mm}",
    )
end

function write_triangle_atlas(io::IO, records, orbit_list, analyses)
    write_triangle_atlas_preamble(io, length.(orbit_list), analyses)
    per_page = 8
    for (orbit_number, orbit) in enumerate(orbit_list)
        for page_start in 1:per_page:length(orbit)
            write_triangle_orbit_header(io, orbit_number, orbit, page_start > 1, analyses)
            page_triangles = orbit[page_start:min(page_start + per_page - 1, length(orbit))]
            for (position, triangle_index) in enumerate(page_triangles)
                lone_last_card = isodd(length(page_triangles)) && position == length(page_triangles)
                lone_last_card && println(io, raw"\hfill")
                write_triangle_tree_card(
                    io,
                    triangle_index,
                    records[triangle_index],
                    analyses[triangle_index],
                )
                if lone_last_card
                    println(io, raw"\hfill\par")
                elseif isodd(position)
                    println(io, raw"\hfill")
                else
                    println(io, raw"\par")
                end
            end
        end
    end
    println(io, raw"\end{document}")
end

function build_balanced_triangle_atlas(
    ;
    fan_path::AbstractString=DEFAULT_BALANCED_FAN_PATH,
    ray_matrix_path::AbstractString=DEFAULT_BALANCED_RAYS_PATH,
    triangle_matrix_path::AbstractString=DEFAULT_BALANCED_TRIANGLES_PATH,
    triangle_indices_path::AbstractString=DEFAULT_TRIANGLE_INDICES_PATH,
    analysis_csv_path::AbstractString=DEFAULT_TRIANGLE_ANALYSIS_CSV_PATH,
    tex_path::AbstractString=DEFAULT_TRIANGLE_ATLAS_TEX_PATH,
    compile_pdf::Bool=true,
    data_only::Bool=false,
)
    fan = load(fan_path)
    balanced_ray_matrix = load(ray_matrix_path)
    records = balanced_triangle_tree_data(fan, balanced_ray_matrix)
    triangles = [record.rays for record in records]

    triangle_matrix = matrix(
        QQ,
        vectors_as_matrix([record.balanced_generator for record in records]),
    )
    triangle_indices = matrix(ZZ, reduce(vcat, permutedims.(triangles)))
    save(triangle_matrix_path, triangle_matrix)
    save(triangle_indices_path, triangle_indices)
    @assert load(triangle_matrix_path) == triangle_matrix
    @assert load(triangle_indices_path) == triangle_indices

    if data_only
        println("primitive integral triangle matrix: ", triangle_matrix_path)
        println("triangle ray-index triples: ", triangle_indices_path)
        println("triangles: ", length(records))
        return (
            fan=fan,
            triangle_matrix=triangle_matrix,
            triangle_indices=triangle_indices,
            records=records,
            analyses=nothing,
            orbits=nothing,
            tex_path=nothing,
            pdf_path=nothing,
        )
    end

    analyses = load_triangle_analysis_csv(analysis_csv_path)
    orbit_list = triangle_orbits(fan)
    @assert length(orbit_list) == 24
    @assert sort!(reduce(vcat, orbit_list)) == collect(1:length(records))
    open(tex_path, "w") do io
        write_triangle_atlas(io, records, orbit_list, analyses)
    end

    pdf_path = compile_pdf ? compile_atlas(tex_path) : nothing
    println("primitive integral triangle matrix: ", triangle_matrix_path)
    println("triangle ray-index triples: ", triangle_indices_path)
    println("triangles: ", length(records))
    println("triangle orbits: ", length.(orbit_list))
    println("triangle analysis CSV: ", analysis_csv_path)
    println("LaTeX atlas: ", tex_path)
    pdf_path === nothing || println("PDF atlas: ", pdf_path)
    return (
        fan=fan,
        triangle_matrix=triangle_matrix,
        triangle_indices=triangle_indices,
        records=records,
        analyses=analyses,
        orbits=orbit_list,
        tex_path=tex_path,
        pdf_path=pdf_path,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    data_only = "--data-only" in ARGS
    build_balanced_triangle_atlas(
        ;
        compile_pdf=!data_only && !("--no-pdf" in ARGS),
        data_only=data_only,
    )
end
