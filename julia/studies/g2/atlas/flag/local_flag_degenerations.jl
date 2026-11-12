# Stage-6 census data shown by the local G2 flag cell atlas: the initial ideal
# of every support cell with the type of its Newton--Okounkov polytope, the two
# Minkowski summands of every prime flag cell, and polymake drawings of the
# three- and four-dimensional smallest pieces of these polytopes.

import Polymake

const CENSUS_TYPES = TropicalCensus.RecursivePolytopeTypes
const CENSUS_DEGENERATIONS = "06_canonical_interior_degenerations.csv"

"The Grassmannian census of a field, `:qq` or `:gf2`."
census_grassmannian_csv(field::Symbol) = TropicalCensus.data_path(
    "g2", field == :qq ? "grassmannian" : "grassmannian_gf2", CENSUS_DEGENERATIONS,
)

# The atlas reads stage-2 flag cell indices; stage 3 keeps their incidence.
"The flag census of a field, `:qq` or `:gf2`."
census_flag_csv(field::Symbol) = TropicalCensus.data_path(
    "g2", field == :qq ? "flag" : "flag_gf2", CENSUS_DEGENERATIONS,
)

# The census tests every Grassmannian body against the Chevalley polytope
# P_{G2/P2, varpi_2, s} of Spacek--Wang, where s = s_2 s_1 s_2 s_1 s_2 is the
# only reduced word of w^P.  That simplex is the lecture hall simplex of this
# sequence with reversed coordinates; the atlas certifies the equivalence for
# every body that the census marks as Chevalley.
const CENSUS_LECTURE_HALL_SEQUENCE = [1, 3, 2, 3, 1]

"""The s-lecture hall simplex conv(0, (0,...,0,s_n), ..., (s_1,...,s_n))."""
function census_lecture_hall_simplex(sequence)
    n = length(sequence)
    return convex_hull(QQ, [
        [index > n - step ? sequence[index] : 0 for index in 1:n] for step in 0:n
    ])
end

census_rotation_x(c, s) = QQ[1 0 0; 0 c -s; 0 s c]
census_rotation_y(c, s) = QQ[c 0 s; 0 1 0; -s 0 c]

# Symmetric realizations drawn for the named three-dimensional pieces, keyed by
# classifier name.  Rational rotations (Pythagorean cosines and sines) turn
# them towards the fixed camera of polymake's TikZ backend.  Other pieces are
# drawn in a rounded linear image of their first occurrence.  Every drawing is
# certified combinatorially isomorphic to the piece in the census.
const CENSUS_PIECE_REALIZATIONS = Dict(
    # The two quadrilaterals share the vertical front edge, like an open book.
    "multiplex M(3, 5)" =>
        QQ[-2 2 0; 0 2 2; 3 2 -1; -3 -2 -1; 0 -2 2; 2 -2 0] * transpose(
            census_rotation_x(12 // 13, 5 // 13) * census_rotation_y(12 // 13, 5 // 13),
        ),
    "4-gonal antiprism" =>
        QQ[2 0 2; -2 0 2; -2 0 -2; 2 0 -2;
           14//5 12//5 0; 0 12//5 14//5; -14//5 12//5 0; 0 12//5 -14//5] *
        transpose(census_rotation_y(4 // 5, -3 // 5)),
)

# Facet fills by number of vertices, with the atlas edge and vertex color.
const CENSUS_FACET_COLORS = Dict(3 => "A9CADB", 4 => "E3A15A")
const CENSUS_OTHER_FACET_COLOR = "8CC7B9"
const CENSUS_EDGE_COLOR = "205B7A"

"""Split one census CSV line, honouring the double quotes of `csv_escape`."""
function census_csv_fields(line)
    characters = collect(line)
    fields = String[]
    buffer = IOBuffer()
    quoted = false
    index = 1
    while index <= length(characters)
        character = characters[index]
        if character == '"'
            if quoted && index < length(characters) && characters[index + 1] == '"'
                print(buffer, '"')
                index += 1
            else
                quoted = !quoted
            end
        elseif character == ',' && !quoted
            push!(fields, String(take!(buffer)))
        else
            print(buffer, character)
        end
        index += 1
    end
    quoted && error("unterminated CSV quote")
    push!(fields, String(take!(buffer)))
    return fields
end

function census_csv_rows(path)
    lines = readlines(path)
    isempty(lines) && error("empty census CSV: $path")
    header = census_csv_fields(first(lines))
    header == TropicalCensus.DEGENERATION_COLUMNS || error("unexpected census header: $path")
    rows = Dict{String,String}[]
    for line in lines[2:end]
        fields = census_csv_fields(line)
        length(fields) == length(header) || error("unexpected census column count: $path")
        push!(rows, Dict(zip(header, fields)))
    end
    return rows
end

"""Classify a census polytope and certify the description stored in the census."""
function census_classified_polytope(path, expected_type)
    polytope = load(path)
    description = CENSUS_TYPES.determine_polytope(polytope)
    recomputed = CENSUS_TYPES.describe_type(description)
    recomputed == expected_type || error(
        "$path is now classified as $recomputed, but the census records $expected_type",
    )
    return (polytope=polytope, type=description)
end

"""The innermost base of a description, below its pyramids, wedges, and prisms."""
function census_smallest_piece(description)
    depth = 0
    while description.kind in (:pyramid, :wedge)
        description = description.base
        depth += 1
    end
    return (type=description, depth=depth)
end

"""Vertices of the smallest piece, as a face of the polytope in its own lattice."""
function census_piece_vertices(polytope, description)
    vertices = Matrix{Rational{BigInt}}(polytope.pm_polytope.VERTICES)
    points = vertices[:, 2:end] ./ vertices[:, 1]
    while description.kind in (:pyramid, :wedge)
        points = points[description.witness.base_vertices, :]
        description = description.base
    end
    return CENSUS_TYPES._full_dimensional(points)
end

function census_named_tex(name)
    name == "point" && return "\\Delta_0"
    name == "segment" && return "\\Delta_1"
    name == "3-gon" && return "\\Delta_2"
    name == "4-gon" && return "\\square"
    simplex = match(r"^(\d+)-simplex$", name)
    simplex === nothing || return "\\Delta_{$(simplex[1])}"
    multiplex = match(r"^multiplex M\((\d+), (\d+)\)$", name)
    multiplex === nothing || return "\\mathrm{M}($(multiplex[1]),$(multiplex[2]))"
    antiprism = match(r"^(\d+)-gonal antiprism$", name)
    antiprism === nothing || return "\\mathrm{AP}_{$(antiprism[1])}"
    text = replace(name, "Δ" => "\\(\\Delta\\)", "×" => "\\(\\times\\)")
    return "\\text{$text}"
end

"""Compact notation for a classifier description; page 3 of the atlas defines it."""
function census_type_tex(description)
    if description.kind == :pyramid
        power = 0
        while description.kind == :pyramid
            description = description.base
            power += 1
        end
        # A k-fold pyramid is the join with a (k-1)-simplex.
        base = census_type_tex(description)
        startswith(base, "\\text{") && (base = "($base)")
        return "$base\\ast\\Delta_{$(power - 1)}"
    elseif description.kind == :wedge
        base = census_type_tex(description.base)
        description.wedged_face.dimension == -1 && return "\\mathrm{Prism}($base)"
        return "\\mathrm{Wed}($base;\\," * census_type_tex(description.wedged_face) * ")"
    elseif description.kind == :dual
        return "(" * census_type_tex(description.base) * ")^{*}"
    elseif description.kind == :unknown
        vertices, facets = size(description.incidence, 2), size(description.incidence, 1)
        return "\\mathrm{U}^{$(description.dimension)}_{$vertices,$facets}"
    end
    return census_named_tex(description.name)
end

census_f_vector_tex(values) = "(" * join(values, ",") * ")"

"""The Newton--Okounkov polytope `P_O` of the support orbit labelled `O`."""
census_polytope_tex(label) = "P_{$label}"

"""Dimensions of the wedged faces anywhere in a description."""
function census_wedge_faces(description)
    description === nothing && return Int[]
    own = description.kind == :wedge ? [description.wedged_face.dimension] : Int[]
    return vcat(own, census_wedge_faces(description.base),
                census_wedge_faces(description.wedged_face))
end

"""
Sort key of a toric support polytope: the Chevalley simplex, then the other
simplices, then by the smallest piece (dimension, recognized before not,
vertices) and by the description (wedges, their faces, operations).
"""
function census_complexity(description, chevalley::Bool)
    piece = census_smallest_piece(description)
    simplex = description.kind == :known && occursin(r"^\d+-simplex$", description.name)
    faces = census_wedge_faces(description)
    return (
        chevalley ? 0 : simplex ? 1 : 2,
        piece.type.dimension,
        piece.type.kind == :unknown ? 1 : 0,
        size(piece.type.incidence, 2),
        length(faces),
        sum(faces; init=0),
        piece.depth,
    )
end

"""
The support orbits of both fields as atlas pages, in order and labelled by
their position. An orbit present in both fields is one page. First come the orbits
with a toric degeneration: in both fields, only over QQ, only over GF(2), each
by increasing complexity. Then the others: present in both fields, only over
QQ, only over GF(2), the prime ones first.
"""
function census_pages(qq_specs, gf2_specs, qq_supports, gf2_supports)
    pages = Any[(spec=spec, qq=orbit, gf2=findfirst(other -> other.qq_orbit == orbit, gf2_specs))
                for (orbit, spec) in enumerate(qq_specs)]
    for (orbit, spec) in enumerate(gf2_specs)
        spec.qq_orbit === nothing && push!(pages, (spec=spec, qq=nothing, gf2=orbit))
    end
    function key(page)
        records = (page.qq === nothing ? nothing : qq_supports[page.qq],
                   page.gf2 === nothing ? nothing : gf2_supports[page.gf2])
        toric = [record !== nothing && record.toric for record in records]
        census = page.qq === nothing ? length(qq_specs) + page.gf2 : page.qq
        if any(toric)
            group = all(toric) ? 1 : toric[1] ? 2 : 3
            return (group, records[findfirst(toric)].complexity..., census)
        end
        presence = page.qq === nothing ? 6 : page.gf2 === nothing ? 5 : 4
        prime = any(record !== nothing && record.prime for record in records)
        return (presence, prime ? 0 : 1, 0, 0, 0, 0, 0, 0, census)
    end
    sort!(pages; by=key)
    return [merge(page, (label=string(position), position=position))
            for (position, page) in enumerate(pages)]
end

"""Status line of a support cell: only the status unless it is prime toric."""
function census_support_summary_tex(support, label)
    support.prime || return raw"\textcolor{nonprimered}{\bfseries NON-PRIME}"
    support.toric || return support.binomial ?
        raw"\textcolor{warningamber}{\bfseries PRIME, NOT TORIC}" :
        raw"\textcolor{warningamber}{\bfseries PRIME NONBINOMIAL}"
    text = raw"\textcolor{toricgreen}{\bfseries PRIME TORIC}" *
           "\\quad \\($(census_polytope_tex(label))\\):" *
           "\\enspace \\($(support.type_tex)\\)" *
           "\\quad \\(f=$(census_f_vector_tex(support.f_vector))\\)" *
           "\\quad $(support.lattice_points) lattice points"
    if support.piece_f_vector !== nothing
        text *= "\\quad \\(f($(support.piece_tex))=" *
                "$(census_f_vector_tex(support.piece_f_vector))\\)"
    end
    if support.chevalley
        sequence = join(CENSUS_LECTURE_HALL_SEQUENCE, ",")
        text *= "\\\\[1mm]" * raw"\textcolor{bodyviolet}{\bfseries CHEVALLEY POLYTOPE}" *
                "\\quad \\(P_{G_2/P_2,\\,\\varpi_2,\\,s}\\) with \\(s=s_2s_1s_2s_1s_2\\), " *
                "unimodularly equivalent to the lecture hall simplex \\(P_{($sequence)}\\)"
    end
    return text
end

function census_piece_label(description)
    description.kind == :unknown && return "unrecognized"
    startswith(description.name, "multiplex ") && return "multiplex"
    description.name == "4-gonal antiprism" && return "square antiprism"
    return description.name
end

"Record the smallest piece of a support polytope (`kind == :support`) or flag summand of `orbit = (field, number)`."
function census_record_piece!(pieces, polytope, description, kind::Symbol, orbit::Tuple{Symbol,Int})
    piece = census_smallest_piece(description)
    piece.type.dimension in (3, 4) || return pieces
    name = CENSUS_TYPES.describe_type(piece.type)
    entry = get!(pieces, name) do
        (
            type=piece.type,
            vertices=census_piece_vertices(polytope, description),
            support_orbits=Set{Tuple{Symbol,Int}}(),
            flag_orbits=Set{Tuple{Symbol,Int}}(),
        )
    end
    push!(kind == :support ? entry.support_orbits : entry.flag_orbits, orbit)
    return pieces
end

census_tikz_rgb(hex) = join(
    (@sprintf("%.4f", parse(Int, hex[index:(index + 1)]; base=16) / 255)
     for index in (1, 3, 5)),
    " ",
)

"""
An exact rational linear image of a point configuration with roughly isotropic
spread: the centred coordinate columns are orthogonalized and scaled to unit
mean square.  It keeps the combinatorial type and keeps projections legible.
"""
function census_rounded_realization(points)
    centred = points .- sum(points; dims=1) ./ size(points, 1)
    columns = [centred[:, column] for column in axes(centred, 2)]
    for later in eachindex(columns), earlier in 1:(later - 1)
        columns[later] -= sum(columns[earlier] .* columns[later]) /
                          sum(columns[earlier] .* columns[earlier]) .* columns[earlier]
    end
    unit(column) = rationalize(
        BigInt, sqrt(length(column) / Float64(sum(column .* column))); tol=1e-3,
    )
    return reduce(hcat, (unit(column) .* column for column in columns))
end

"""
Write a polymake visual as TikZ, scaled to fit `width` by `height` centimetres.
The header with user name and date is dropped, so that the atlas stays
reproducible.
"""
function census_tikz_code(visual, name::String, vertices::Int; width::Float64, height::Float64)
    code = mktempdir() do directory
        path = joinpath(directory, "$name.tikz")
        Polymake.call_function(Nothing, :common, :tikz, visual; File=path)
        read(path, String)
    end
    start = findfirst("\\begin{tikzpicture}", code)
    start === nothing && error("polymake wrote no TikZ picture for $name")
    code = code[first(start):end]
    axis = Dict(
        String(found[1]) => parse.(Float64, (found[2], found[3]))
        for found in eachmatch(r"\b([xyz])\s*=\s*\{\(([-0-9.e]+)cm,\s*([-0-9.e]+)cm\)\}", code)
    )
    Set(keys(axis)) == Set(["x", "y", "z"]) || error("polymake TikZ axes not found for $name")
    projected = [
        let point = parse.(Float64, split(found[1], ','))
            Tuple(sum(point[index] * axis[label][component]
                      for (index, label) in enumerate(("x", "y", "z")))
                  for component in 1:2)
        end
        for found in eachmatch(r"\\coordinate \([^)]*\) at \(([^)]*)\);", code)
    ]
    length(projected) == vertices || error("polymake TikZ vertices not found for $name")
    extent(component) = maximum(point[component] for point in projected) -
                        minimum(point[component] for point in projected)
    scale = min(width / extent(1), height / extent(2))
    occursin("scale = 1,", code) || error("polymake TikZ scale not found for $name")
    return replace(code, "scale = 1," => @sprintf("scale = %.4f,", scale); count=1)
end

"""polymake's solid drawing of a 3-polytope, its facets filled by vertex count."""
function census_solid_tikz(polytope, name::String; width::Float64, height::Float64)
    pm = polytope.pm_polytope
    incidence = Matrix{Bool}(pm.VERTICES_IN_FACETS)
    facet_colors = [
        census_tikz_rgb(get(CENSUS_FACET_COLORS, count(incidence[facet, :]), CENSUS_OTHER_FACET_COLOR))
        for facet in axes(incidence, 1)
    ]
    visual = Polymake.visual(
        Polymake.Visual, pm;
        Name=name,
        FacetColor=facet_colors,
        FacetTransparency=0.8,
        EdgeColor=census_tikz_rgb(CENSUS_EDGE_COLOR),
        VertexColor=census_tikz_rgb(CENSUS_EDGE_COLOR),
        VertexThickness=2.5,
        VertexLabels="hidden",
    )
    return census_tikz_code(visual, name, size(incidence, 2); width, height)
end

"""
polymake's Schlegel diagram of a 4-polytope, projected into the first facet
with the most vertices; polymake draws that facet's edges in black.
"""
function census_schlegel_tikz(polytope, name::String; width::Float64, height::Float64)
    pm = polytope.pm_polytope
    incidence = Matrix{Bool}(pm.VERTICES_IN_FACETS)
    facet = argmax([count(incidence[row, :]) for row in axes(incidence, 1)]) - 1
    # SCHLEGEL takes its diagram, projection-facet, and edge options separately.
    visual = Polymake.call_method(
        pm, :SCHLEGEL,
        Polymake.OptionSet(Dict(:FACET => facet)),
        Polymake.OptionSet(Dict{Symbol,String}()),
        Polymake.OptionSet(Dict(
            :Name => name,
            :EdgeColor => census_tikz_rgb(CENSUS_EDGE_COLOR),
            :VertexColor => census_tikz_rgb(CENSUS_EDGE_COLOR),
            :VertexThickness => 2.5,
            :VertexLabels => "hidden",
        )),
    )
    return census_tikz_code(visual, name, size(incidence, 2); width, height)
end

"""
Draw the recorded pieces. `position` gives the page position of a support
orbit `(field, number)`; each piece lists the positions of its pages.
"""
function census_drawn_pieces(pieces, position)
    pages(name) = sort!(unique(position[orbit] for orbit in pieces[name].support_orbits))
    order = sort!(collect(keys(pieces)); by=name -> (
        pieces[name].type.dimension,
        isempty(pages(name)) ? typemax(Int) : first(pages(name)),
        name,
    ))
    result = Any[]
    for (index, name) in enumerate(order)
        entry = pieces[name]
        dimension = entry.type.dimension
        vertices = get(CENSUS_PIECE_REALIZATIONS, name) do
            census_rounded_realization(entry.vertices)
        end
        polytope = convex_hull(QQ, vertices)
        census = convex_hull(QQ, entry.vertices)
        drawn = CENSUS_TYPES.describe_type(CENSUS_TYPES.determine_polytope(polytope))
        drawn == name &&
            Bool(Polymake.polytope.isomorphic(polytope.pm_polytope, census.pm_polytope)) ||
            error("the drawn realization of $name is not combinatorially the census piece")
        tikz = dimension == 3 ?
            census_solid_tikz(polytope, "piece$index"; width=12.5, height=12.5) :
            census_schlegel_tikz(polytope, "piece$index"; width=7.6, height=7.6)
        push!(result, (
            name=name,
            dimension=dimension,
            label=census_piece_label(entry.type),
            tex=census_type_tex(entry.type),
            f_vector=Int.(f_vector(polytope)),
            pages=pages(name),
            tikz=tikz,
        ))
    end
    return result
end

"""
Attach the stage-6 censuses of one field to its support specs and flag fan:
the support record of every support orbit, its Newton--Okounkov polytope when
it is toric, and the Minkowski summands of every prime flag cell. `ray_index`
gives the atlas index of each census Grassmannian ray. The smallest pieces of
dimension three and four are added to `pieces`. The Pluecker summand of a
prime flag cell is certified to be the polytope of its support orbit.
"""
function census_field_data(field::Symbol, specs, flag_data, pieces; ray_index=nothing)
    atlas_rays(rays) = Tuple(sort!(ray_index === nothing ? rays : ray_index[rays]))
    supports = Any[]
    bodies = Dict{Int,Any}()
    grassmannian_csv = census_grassmannian_csv(field)
    grassmannian_rows = census_csv_rows(grassmannian_csv)
    length(grassmannian_rows) == length(specs) ||
        error("the $field Grassmannian census does not have one row per support orbit")
    for (orbit_number, (spec, row)) in enumerate(zip(specs, grassmannian_rows))
        members = Set(Tuple(sort!(collect(member))) for member in spec.members)
        parse(Int, row["orbit"]) == orbit_number &&
            atlas_rays(parse_csv_ints(row["rays"])) in members &&
            parse(Int, row["orbit_size"]) == length(spec.members) || error(
            "$field Grassmannian census orbit $orbit_number is not support orbit $(spec.representative)",
        )
        support = (
            orbit=orbit_number,
            prime=parse_csv_bool(row["prime"]),
            binomial=parse_csv_bool(row["binomial"]),
            toric=parse_csv_bool(row["toric"]),
            chevalley=row["chevalley_equivalent"] == "true",
        )
        support.chevalley && !support.toric &&
            error("$field Grassmannian census orbit $orbit_number is Chevalley but not toric")
        if support.toric
            body = census_classified_polytope(
                joinpath(dirname(grassmannian_csv), row["newton_okounkov_body"]),
                row["body_type"],
            )
            support.chevalley && !TropicalCensus.affine_unimodular_equivalent(
                body.polytope, census_lecture_hall_simplex(CENSUS_LECTURE_HALL_SEQUENCE),
            ) && error(
                "the Chevalley polytope of $field support orbit $orbit_number is not the " *
                "lecture hall simplex",
            )
            piece = census_smallest_piece(body.type)
            piece_f_vector = piece.type.kind == :unknown ?
                Int.(f_vector(convex_hull(QQ, census_piece_vertices(body.polytope, body.type)))) :
                nothing
            census_record_piece!(pieces, body.polytope, body.type, :support, (field, orbit_number))
            bodies[orbit_number] = body.polytope
            support = merge(support, (
                f_vector=parse_csv_ints(row["body_f_vector"]),
                lattice_points=parse(Int, row["body_lattice_points"]),
                type_tex=census_type_tex(body.type),
                piece_tex=census_type_tex(piece.type),
                piece_f_vector=piece_f_vector,
                complexity=census_complexity(body.type, support.chevalley),
            ))
        end
        push!(supports, support)
    end

    support_of = Dict(
        Tuple(sort!(collect(member))) => orbit_number
        for (orbit_number, spec) in enumerate(specs) for member in spec.members
    )
    summands = Dict{Int,Any}()
    covered = Int[]
    flag_csv = census_flag_csv(field)
    for row in census_csv_rows(flag_csv)
        orbit = parse(Int, row["orbit"])
        cells = parse_csv_ints(row["cells"])
        representative = parse(Int, row["representative_cell"])
        findall(flag_data.incidence[representative, :]) == parse_csv_ints(row["rays"]) ||
            error("$field flag census orbit $orbit uses another cell numbering than the atlas fan")
        append!(covered, cells)
        parse_csv_bool(row["prime"]) || continue
        support = only(unique(support_of[Tuple(flag_data.projected[cell])] for cell in cells))
        summand = (slice_tex=nothing, support=support)
        if parse_csv_bool(row["toric"])
            body = joinpath(dirname(flag_csv), row["newton_okounkov_body"])
            types = split(row["body_slice_types"], " + ")
            length(types) == 2 || error("$field flag census orbit $orbit does not have two summands")
            slices = [
                census_classified_polytope(
                    replace(body, r"\.poly$" => "_slice_$index.poly"), types[index],
                )
                for index in 1:2
            ]
            for slice in slices
                census_record_piece!(pieces, slice.polytope, slice.type, :flag, (field, orbit))
            end
            haskey(bodies, support) && TropicalCensus.affine_unimodular_equivalent(
                slices[2].polytope, bodies[support],
            ) || error(
                "the Pluecker summand of $field flag orbit $orbit is not the polytope of " *
                "its support orbit $support",
            )
            summand = (slice_tex=census_type_tex(slices[1].type), support=support)
        end
        foreach(cell -> summands[cell] = summand, cells)
    end
    sort(covered) == collect(axes(flag_data.incidence, 1)) ||
        error("the $field flag census orbits do not partition the atlas cells")
    return (supports=supports, bodies=bodies, summands=summands)
end

"A toric support orbit present in both fields has the same polytope in both."
function census_verify_common_polytopes(pages, qq, gf2)
    for page in pages
        (page.qq === nothing || page.gf2 === nothing) && continue
        haskey(qq.bodies, page.qq) && haskey(gf2.bodies, page.gf2) || continue
        TropicalCensus.affine_unimodular_equivalent(qq.bodies[page.qq], gf2.bodies[page.gf2]) ||
            error("support orbit $(page.label) has different polytopes over QQ and over GF(2)")
    end
    return true
end
