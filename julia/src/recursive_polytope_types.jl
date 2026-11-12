module RecursivePolytopeTypes

import Polymake

export determine_5_polytope, determine_polytope, describe_type,
       print_type, unrecognized_parts, PolytopeType

"""
A certified combinatorial description. `base` and `wedged_face` are recursive
descriptions; `incidence` records the complete facet-by-vertex incidence matrix.
For kind `:dual`, `base` describes the polar, whose incidence is the transpose:
vertex `j` of the polar is facet `j` here, facet `i` of the polar is vertex `i`.
All witness indices are Julia indices (one-based), local to this node or its
base as specified by their field names. The description need not be unique.
"""
struct PolytopeType
    kind::Symbol
    name::String
    dimension::Int
    incidence::BitMatrix
    base::Union{Nothing,PolytopeType}
    wedged_face::Union{Nothing,PolytopeType}
    witness::NamedTuple
    fully_named::Bool
end

_incidence(pm) = BitMatrix(Matrix{Bool}(pm.VERTICES_IN_FACETS))

function _pm_polytope(I::BitMatrix, d::Int)
    Polymake.polytope.Polytope{Polymake.Rational}(
        VERTICES_IN_FACETS=Polymake.IncidenceMatrix(I),
        COMBINATORIAL_DIM=d, N_VERTICES=size(I, 2), N_FACETS=size(I, 1),
    )
end

_signature(I) = (size(I), Tuple(vec(I)))
_row_signature(I) = sort([Tuple(I[i, :]) for i in axes(I, 1)])
_column_degrees(I) = sort(vec(sum(I; dims=1)))
_row_degrees(I) = sort(vec(sum(I; dims=2)))

function _same_type(I::BitMatrix, J::BitMatrix, d::Int)
    size(I) == size(J) || return false
    _row_degrees(I) == _row_degrees(J) || return false
    _column_degrees(I) == _column_degrees(J) || return false
    Bool(Polymake.polytope.isomorphic(_pm_polytope(I, d), _pm_polytope(J, d)))
end

const _CATALOG_CACHE = Dict{Any,BitMatrix}()

function _cached(builder, key)
    get!(_CATALOG_CACHE, key) do
        BitMatrix(builder())
    end
end

function _simplex(d)
    I = trues(d + 1, d + 1)
    for i in 1:(d + 1)
        I[i, i] = false
    end
    I
end

function _cube(d)
    I = falses(2d, 2^d)
    for v in 0:(2^d - 1), j in 1:d
        bit = (v >> (j - 1)) & 1
        I[2j - 1 + bit, v + 1] = true
    end
    I
end

function _cross(d)
    I = falses(2^d, 2d)
    for f in 0:(2^d - 1), j in 1:d
        bit = (f >> (j - 1)) & 1
        I[f + 1, 2j - 1 + bit] = true
    end
    I
end

function _multiplex(d, nv)
    I = falses(nv, nv)
    for i in 1:nv, k in 1:(d - 1)
        I[i, clamp(i - k, 1, nv)] = true
        I[i, clamp(i + k, 1, nv)] = true
    end
    I
end

function _partitions(n, least=1)
    n == 0 && return [Int[]]
    result = Vector{Vector{Int}}()
    for k in least:n, tail in _partitions(n - k, k)
        push!(result, [k; tail])
    end
    result
end

function _simplex_product(parts)
    vertices = collect(Iterators.product((0:k for k in parts)...))
    I = falses(sum(parts .+ 1), length(vertices))
    row = 0
    for (factor, k) in enumerate(parts), missing in 0:k
        row += 1
        for (v, labels) in enumerate(vertices)
            I[row, v] = labels[factor] != missing
        end
    end
    I
end

function _polygon(n)
    I = falses(n, n)
    for i in 1:n
        I[i, i] = I[i, mod1(i + 1, n)] = true
    end
    I
end

function _polygon_bipyramid(n)
    polygon = _polygon(n)
    I = falses(2n, n + 2)
    for i in 1:n, apex in 1:2
        row = 2i - 2 + apex
        I[row, 1:n] = polygon[i, :]
        I[row, n + apex] = true
    end
    I
end

function _polygon_prism(n)
    polygon = _polygon(n)
    I = falses(n + 2, 2n)
    I[1, 1:n] .= true
    I[2, (n + 1):2n] .= true
    for i in 1:n
        I[i + 2, 1:n] = polygon[i, :]
        I[i + 2, (n + 1):2n] = polygon[i, :]
    end
    I
end

"""Antiprism over an n-gon; columns: top n-gon, then bottom n-gon."""
function _antiprism(n)
    I = falses(2n + 2, 2n)
    I[1, 1:n] .= true
    I[2, (n + 1):2n] .= true
    for i in 1:n
        j = mod1(i + 1, n)
        I[2i + 1, [i, j, n + i]] .= true
        I[2i + 2, [n + i, n + j, j]] .= true
    end
    I
end

"""Product incidence: rows `F × B`, then `A × G`; vertex `(i, j)` is column `(i - 1)nb + j`."""
_product(A, B) = BitMatrix([repeat(A; inner=(1, size(B, 2)));
                            repeat(B; outer=(1, size(A, 2)))])

"""Zonotope of the four cube diagonals and one coordinate axis."""
function _elongated_dodecahedron()
    generators = [1 1 1; 1 -1 1; -1 1 1; -1 -1 1; 0 0 1]
    points = reduce(vcat, [[1 collect(s)' * generators]
                           for s in Iterators.product(fill((-1, 1), 5)...)])
    Polymake.polytope.Polytope{Polymake.Rational}(POINTS=points)
end

"""
Exact points of the 600-cell of circumradius 2: first the 24-cell
`(±2, 0, 0, 0)`, `(±1, ±1, ±1, ±1)`, then even permutations of `(±φ, ±1, ±1/φ, 0)`.
"""
function _600_cell_points()
    QE = Polymake.QuadraticExtension{Polymake.Rational}
    φ = QE(1//2, 1//2, 5)
    signs(k) = Iterators.product(fill((-1, 1), k)...)
    even_permutations = [collect(p) for p in Iterators.product(ntuple(_ -> 1:4, 4)...)
        if allunique(p) && iseven(count(p[i] > p[j] for i in 1:4 for j in (i + 1):4))]
    points = Vector{Vector{QE}}()
    for i in 1:4, s in (-2, 2)
        push!(points, [QE(j == i ? s : 0) for j in 1:4])
    end
    for s in signs(4)
        push!(points, QE.(collect(s)))
    end
    for p in even_permutations, s in signs(3)
        push!(points, [QE(s[1]) * φ, QE(s[2]), QE(s[3]) * (φ - 1), QE(0)][p])
    end
    points
end

function _convex_hull(points)
    QE = Polymake.QuadraticExtension{Polymake.Rational}
    M = reduce(vcat, [permutedims([one(QE); x]) for x in points])
    Polymake.polytope.Polytope{QE}(POINTS=M)
end

"""The 600-cell without the vertices of an inscribed 24-cell."""
_snub_24_cell() = _convex_hull(_600_cell_points()[25:end])

"""The 600-cell without two decagons in completely orthogonal planes."""
function _grand_antiprism()
    φ = Polymake.QuadraticExtension{Polymake.Rational}(1//2, 1//2, 5)
    # The planes span(e1, (0, φ, 1, 0)) and its orthogonal complement.
    on_planes(x) = (iszero(x[4]) && x[2] == φ * x[3]) ||
                   (iszero(x[1]) && iszero(φ * x[2] + x[3]))
    _convex_hull(filter(!on_planes, _600_cell_points()))
end

# Archimedean solids by polymake name, with their duals, the Catalan solids.
const _ARCHIMEDEAN_SOLIDS = (
    ("truncated_tetrahedron", "triakis tetrahedron"),
    ("cuboctahedron", "rhombic dodecahedron"),
    ("truncated_cube", "triakis octahedron"),
    ("truncated_octahedron", "tetrakis hexahedron"),
    ("rhombicuboctahedron", "deltoidal icositetrahedron"),
    ("truncated_cuboctahedron", "disdyakis dodecahedron"),
    ("snub_cube", "pentagonal icositetrahedron"),
    ("icosidodecahedron", "rhombic triacontahedron"),
    ("truncated_dodecahedron", "triakis icosahedron"),
    ("truncated_icosahedron", "pentakis dodecahedron"),
    ("rhombicosidodecahedron", "deltoidal hexecontahedron"),
    ("truncated_icosidodecahedron", "disdyakis triacontahedron"),
    ("snub_dodecahedron", "pentagonal hexecontahedron"),
)

# Johnson solids J1, ..., J92, in the numbering of polymake's `johnson_solid`.
const _JOHNSON_NAMES = (
    "square pyramid", "pentagonal pyramid", "triangular cupola",
    "square cupola", "pentagonal cupola", "pentagonal rotunda",
    "elongated triangular pyramid", "elongated square pyramid",
    "elongated pentagonal pyramid", "gyroelongated square pyramid",
    "gyroelongated pentagonal pyramid", "triangular bipyramid",
    "pentagonal bipyramid", "elongated triangular bipyramid",
    "elongated square bipyramid", "elongated pentagonal bipyramid",
    "gyroelongated square bipyramid", "elongated triangular cupola",
    "elongated square cupola", "elongated pentagonal cupola",
    "elongated pentagonal rotunda", "gyroelongated triangular cupola",
    "gyroelongated square cupola", "gyroelongated pentagonal cupola",
    "gyroelongated pentagonal rotunda", "gyrobifastigium",
    "triangular orthobicupola", "square orthobicupola", "square gyrobicupola",
    "pentagonal orthobicupola", "pentagonal gyrobicupola",
    "pentagonal orthocupolarotunda", "pentagonal gyrocupolarotunda",
    "pentagonal orthobirotunda", "elongated triangular orthobicupola",
    "elongated triangular gyrobicupola", "elongated square gyrobicupola",
    "elongated pentagonal orthobicupola", "elongated pentagonal gyrobicupola",
    "elongated pentagonal orthocupolarotunda",
    "elongated pentagonal gyrocupolarotunda",
    "elongated pentagonal orthobirotunda", "elongated pentagonal gyrobirotunda",
    "gyroelongated triangular bicupola", "gyroelongated square bicupola",
    "gyroelongated pentagonal bicupola", "gyroelongated pentagonal cupolarotunda",
    "gyroelongated pentagonal birotunda", "augmented triangular prism",
    "biaugmented triangular prism", "triaugmented triangular prism",
    "augmented pentagonal prism", "biaugmented pentagonal prism",
    "augmented hexagonal prism", "parabiaugmented hexagonal prism",
    "metabiaugmented hexagonal prism", "triaugmented hexagonal prism",
    "augmented dodecahedron", "parabiaugmented dodecahedron",
    "metabiaugmented dodecahedron", "triaugmented dodecahedron",
    "metabidiminished icosahedron", "tridiminished icosahedron",
    "augmented tridiminished icosahedron", "augmented truncated tetrahedron",
    "augmented truncated cube", "biaugmented truncated cube",
    "augmented truncated dodecahedron", "parabiaugmented truncated dodecahedron",
    "metabiaugmented truncated dodecahedron",
    "triaugmented truncated dodecahedron", "gyrate rhombicosidodecahedron",
    "parabigyrate rhombicosidodecahedron", "metabigyrate rhombicosidodecahedron",
    "trigyrate rhombicosidodecahedron", "diminished rhombicosidodecahedron",
    "paragyrate diminished rhombicosidodecahedron",
    "metagyrate diminished rhombicosidodecahedron",
    "bigyrate diminished rhombicosidodecahedron",
    "parabidiminished rhombicosidodecahedron",
    "metabidiminished rhombicosidodecahedron",
    "gyrate bidiminished rhombicosidodecahedron",
    "tridiminished rhombicosidodecahedron", "snub disphenoid",
    "snub square antiprism", "sphenocorona", "augmented sphenocorona",
    "sphenomegacorona", "hebesphenomegacorona", "disphenocingulum",
    "bilunabirotunda", "triangular hebesphenorotunda",
)

# Wythoffian uniform 4-polytopes not named elsewhere: (name, Coxeter type,
# rings, nv, nf). Ring 0 is the 5-cell, 16-cell, 24-cell, 120-cell end.
const _UNIFORM_4_POLYTOPES = (
    ("truncated 5-cell", "A4", (0, 1), 20, 10),
    ("cantellated 5-cell", "A4", (0, 2), 30, 20),
    ("runcinated 5-cell", "A4", (0, 3), 20, 30),
    ("bitruncated 5-cell", "A4", (1, 2), 30, 10),
    ("cantitruncated 5-cell", "A4", (0, 1, 2), 60, 20),
    ("runcitruncated 5-cell", "A4", (0, 1, 3), 60, 30),
    ("rectified tesseract", "B4", (2,), 32, 24),
    ("truncated tesseract", "B4", (2, 3), 64, 24),
    ("cantellated tesseract", "B4", (1, 3), 96, 56),
    ("runcinated tesseract", "B4", (0, 3), 64, 80),
    ("bitruncated tesseract", "B4", (1, 2), 96, 24),
    ("cantitruncated tesseract", "B4", (1, 2, 3), 192, 56),
    ("runcitruncated tesseract", "B4", (0, 2, 3), 192, 80),
    ("omnitruncated tesseract", "B4", (0, 1, 2, 3), 384, 80),
    ("truncated 16-cell", "B4", (0, 1), 48, 24),
    ("runcitruncated 16-cell", "B4", (0, 1, 3), 192, 80),
    ("rectified 24-cell", "F4", (1,), 96, 48),
    ("truncated 24-cell", "F4", (0, 1), 192, 48),
    ("cantellated 24-cell", "F4", (0, 2), 288, 144),
    ("runcinated 24-cell", "F4", (0, 3), 144, 240),
    ("bitruncated 24-cell", "F4", (1, 2), 288, 48),
    ("cantitruncated 24-cell", "F4", (0, 1, 2), 576, 144),
    ("runcitruncated 24-cell", "F4", (0, 1, 3), 576, 240),
    ("omnitruncated 24-cell", "F4", (0, 1, 2, 3), 1152, 240),
    ("rectified 120-cell", "H4", (1,), 1200, 720),
    ("rectified 600-cell", "H4", (2,), 720, 720),
    ("truncated 120-cell", "H4", (0, 1), 2400, 720),
    ("truncated 600-cell", "H4", (2, 3), 1440, 720),
    ("cantellated 120-cell", "H4", (0, 2), 3600, 1920),
    ("cantellated 600-cell", "H4", (1, 3), 3600, 1440),
    ("runcinated 120-cell", "H4", (0, 3), 2400, 2640),
    ("bitruncated 120-cell", "H4", (1, 2), 3600, 720),
    ("cantitruncated 120-cell", "H4", (0, 1, 2), 7200, 1920),
    ("cantitruncated 600-cell", "H4", (1, 2, 3), 7200, 1440),
    ("runcitruncated 120-cell", "H4", (0, 1, 3), 7200, 2640),
    ("runcitruncated 600-cell", "H4", (0, 2, 3), 7200, 2640),
    ("omnitruncated 120-cell", "H4", (0, 1, 2, 3), 14400, 2640),
)

function _known_name(I, d, extra_known; extended_catalog=true)
    nf, nv = size(I)
    d == -1 && return "empty face"
    d == 0 && return "point"
    d == 1 && return "segment"
    # Every bounded convex 2-polytope is combinatorially an nv-gon.
    d == 2 && return "$(nv)-gon"

    # Count filters avoid generating candidates which cannot match.
    nv == d + 1 && nf == d + 1 && _same_type(I, _simplex(d), d) &&
        return "$(d)-simplex"
    nv == 2^d && nf == 2d && _same_type(I, _cube(d), d) &&
        return "$(d)-cube"
    nv == 2d && nf == 2^d && _same_type(I, _cross(d), d) &&
        return "$(d)-cross-polytope"
    nv == nf && nv >= d + 1 && _same_type(I, _multiplex(d, nv), d) &&
        return "multiplex M($(d), $(nv - 1))"

    for parts in _partitions(d)
        length(parts) >= 2 || continue
        prod(parts .+ 1) == nv && sum(parts .+ 1) == nf || continue
        J = _cached(() -> _simplex_product(parts), (:product, Tuple(parts)))
        _same_type(I, J, d) && return join(["Δ$(k)" for k in parts], " × ")
    end

    # A cyclic d-polytope is simplicial and has the maximal facet count.
    m = d ÷ 2
    cyclic_facets = iseven(d) ?
        binomial(big(nv - m), m) + binomial(big(nv - m - 1), m - 1) :
        2binomial(big(nv - m - 1), m)
    if all(==(d), vec(sum(I; dims=2))) && nf == cyclic_facets
        J = _cached(() -> _incidence(Polymake.polytope.cyclic(d, nv)),
                    (:cyclic, d, nv))
        _same_type(I, J, d) && return "cyclic polytope C($(d), $(nv))"
    end

    if d == 3
        for (name, constructor, expected_nv, expected_nf) in (
            ("dodecahedron", () -> Polymake.polytope.dodecahedron(), 20, 12),
            ("icosahedron", () -> Polymake.polytope.icosahedron(), 12, 20),
        )
            nv == expected_nv && nf == expected_nf || continue
            J = _cached(() -> _incidence(constructor()), (name,))
            _same_type(I, J, d) && return name
        end
        if iseven(nv) && nv >= 6 && nf == nv ÷ 2 + 2
            _same_type(I, _polygon_prism(nv ÷ 2), d) &&
                return "$(nv ÷ 2)-gonal prism"
        end
        if nv >= 5 && nf == 2(nv - 2)
            _same_type(I, _polygon_bipyramid(nv - 2), d) &&
                return "$(nv - 2)-gonal bipyramid"
        end
        if iseven(nv) && nv >= 6 && nf == nv + 2
            _same_type(I, _antiprism(nv ÷ 2), d) &&
                return "$(nv ÷ 2)-gonal antiprism"
        end
        if iseven(nf) && nf >= 6 && nv == nf + 2
            _same_type(I, BitMatrix(transpose(_antiprism(nf ÷ 2))), d) &&
                return "$(nf ÷ 2)-gonal trapezohedron"
        end
        if extended_catalog
            for (name, dual_name) in _ARCHIMEDEAN_SOLIDS
                J = _cached(() -> _incidence(Polymake.polytope.archimedean_solid(name)),
                            (:archimedean, name))
                _same_type(I, J, d) && return replace(name, "_" => " ")
                _same_type(I, BitMatrix(transpose(J)), d) && return dual_name
            end
            for (k, name) in enumerate(_JOHNSON_NAMES)
                J = _cached(() -> _incidence(Polymake.polytope.johnson_solid(k)),
                            (:johnson, k))
                _same_type(I, J, d) && return "$(name) (J$(k))"
            end
            if nv == 18 && nf == 12
                J = _cached(() -> _incidence(_elongated_dodecahedron()),
                            (:elongated_dodecahedron,))
                _same_type(I, J, d) && return "elongated dodecahedron"
            end
        end
    end

    if d == 4
        # Small cases are Δ2 × Δ2, Δ1 × Δ1 × Δ2, the 4-cube, C(4, 6), the 16-cell.
        for p in 3:(nf ÷ 2)
            q = nf - p
            p * q == nv || continue
            _same_type(I, _product(_polygon(p), _polygon(q)), d) &&
                return "$(p),$(q)-duoprism"
        end
        for p in 3:(nv ÷ 2)
            q = nv - p
            p * q == nf || continue
            _same_type(I, BitMatrix(transpose(_product(_polygon(p), _polygon(q)))), d) &&
                return "$(p),$(q)-duopyramid"
        end
    end

    if extended_catalog
        if d == 4
            for (name, constructor, expected_nv, expected_nf) in (
                ("24-cell", () -> Polymake.polytope.regular_24_cell(), 24, 24),
                ("120-cell", () -> Polymake.polytope.regular_120_cell(), 600, 120),
                ("600-cell", () -> Polymake.polytope.regular_600_cell(), 120, 600),
                ("snub 24-cell", _snub_24_cell, 96, 144),
                ("grand antiprism", _grand_antiprism, 100, 320),
            )
                nv == expected_nv && nf == expected_nf || continue
                J = _cached(() -> _incidence(constructor()), (:exceptional4, name))
                _same_type(I, J, d) && return name
            end
            for (name, type, rings, expected_nv, expected_nf) in _UNIFORM_4_POLYTOPES
                nv == expected_nv && nf == expected_nf || continue
                J = _cached(() -> _incidence(Polymake.polytope.wythoff(
                                type, Polymake.Set{Int}(collect(rings)))),
                            (:wythoff, type, rings))
                _same_type(I, J, d) && return name
            end
        end
        for k in 2:((d + 1) ÷ 2)
            nv == binomial(d + 1, k) || continue
            J = _cached(() -> _incidence(Polymake.polytope.hypersimplex(k, d + 1)),
                        (:hypersimplex, k, d + 1))
            _same_type(I, J, d) && return "hypersimplex Δ($(k), $(d + 1))"
        end
        if nv == factorial(big(d + 1)) && nf == 2^(d + 1) - 2
            J = _cached(() -> _incidence(Polymake.polytope.permutahedron(d)),
                        (:permutahedron, d))
            _same_type(I, J, d) && return "$(d)-permutahedron"
        end
        catalan_number = binomial(big(2d + 2), d + 1) ÷ (d + 2)
        if nv == catalan_number && nf == d * (d + 3) ÷ 2
            J = _cached(() -> _incidence(Polymake.polytope.associahedron(d)),
                        (:associahedron, d))
            _same_type(I, J, d) && return "$(d)-associahedron"
        end
    end

    for (name, J, candidate_dimension) in extra_known
        candidate_dimension == d && _same_type(I, J, d) && return name
    end
    nothing
end

"""Incidence of a face, with vertices ordered as in `vertices`."""
function _face_incidence(I, vertices)
    isempty(vertices) && return falses(0, 0)
    length(vertices) == 1 && return falses(0, 1)
    rows = unique([Tuple(I[i, vertices]) for i in axes(I, 1)])
    filter!(r -> any(r) && !all(r), rows)
    maximal = [r for r in rows if !any(
        s -> r != s && all(!r[j] || s[j] for j in eachindex(r)), rows)]
    J = falses(length(maximal), length(vertices))
    for (i, row) in enumerate(maximal)
        J[i, :] = collect(row)
    end
    J
end

"""Exact dimension from a maximal chain of faces in a valid incidence matrix."""
function _incidence_dimension(I)
    nv = size(I, 2)
    nv == 0 && return -1
    nv == 1 && return 0
    size(I, 1) > 0 || throw(ArgumentError("Invalid face incidence."))
    f = argmin(vec(sum(I; dims=2)))
    1 + _incidence_dimension(_face_incidence(I, findall(I[f, :])))
end

"""Containing facets of a proper face mask; `nothing` if it is not a face."""
function _proper_face_facets(Q, fixed)
    any(!, fixed) || return nothing # exclude the entire base
    containing = [i for i in axes(Q, 1) if all(.!fixed .| Q[i, :])]
    # The empty face is allowed, including for a point base.
    any(fixed) || return containing
    closure = trues(size(Q, 2))
    for i in containing
        closure .&= Q[i, :]
    end
    closure == fixed ? containing : nothing
end

"""Face wedge; columns: old vertices, then clones outside the marked face."""
function _wedge_incidence(Q, face_vertices::AbstractVector{<:Integer})
    nf, nv = size(Q)
    all(v -> 1 <= v <= nv, face_vertices) ||
        throw(ArgumentError("Face vertex indices must be in 1:$nv."))
    length(unique(face_vertices)) == length(face_vertices) ||
        throw(ArgumentError("Face vertex indices must be distinct."))
    fixed = falses(nv)
    fixed[face_vertices] .= true
    _proper_face_facets(Q, fixed) !== nothing ||
        throw(ArgumentError("The marked vertices must form a proper face."))
    outside = findall(!, fixed)
    sides = [i for i in 1:nf if any(Q[i, outside])]
    W = falses(2 + length(sides), nv + length(outside))
    W[1, 1:nv] .= true
    W[2, 1:nv] = fixed
    W[2, (nv + 1):end] .= true
    for (j, i) in enumerate(sides)
        W[j + 2, 1:nv] = Q[i, :]
        W[j + 2, (nv + 1):end] = Q[i, outside]
    end
    W
end

_wedge_incidence(Q, facet_index::Integer) =
    _wedge_incidence(Q, findall(Q[facet_index, :]))

function _inverse_wedge(I, a, b)
    nf, nv = size(I)
    all(I[a, :] .| I[b, :]) || return nothing
    shared = findall(I[a, :] .& I[b, :])
    U = findall(I[a, :])
    left = findall(I[a, :] .& .!I[b, :])
    right = findall(I[b, :] .& .!I[a, :])
    isempty(left) && return nothing
    length(left) == length(right) || return nothing
    side_facets = [i for i in 1:nf if i != a && i != b]

    left_signatures = [Tuple(I[side_facets, v]) for v in left]
    right_signatures = [Tuple(I[side_facets, v]) for v in right]
    length(unique(left_signatures)) == length(left) || return nothing
    length(unique(right_signatures)) == length(right) || return nothing
    right_by_signature = Dict(s => v for (s, v) in zip(right_signatures, right))
    all(haskey(right_by_signature, s) for s in left_signatures) || return nothing
    clones = [right_by_signature[s] for s in left_signatures]

    Q = _face_incidence(I, U)
    shared_in_base = [i for (i, v) in enumerate(U) if v in shared]
    fixed = falses(length(U))
    fixed[shared_in_base] .= true
    containing = _proper_face_facets(Q, fixed)
    containing === nothing && return nothing
    f = findfirst(i -> Q[i, :] == fixed, axes(Q, 1))

    # Verify the entire incidence, under the explicit vertex correspondence.
    W = _wedge_incidence(Q, shared_in_base)
    vertex_order = [U; clones]
    length(vertex_order) == nv || return nothing
    _row_signature(W) == _row_signature(I[:, vertex_order]) || return nothing
    face = _face_incidence(Q, shared_in_base)
    face_dimension = _incidence_dimension(face)
    (base=Q, face=face, face_dimension=face_dimension,
     witness=(base_facets=(a, b), base_vertices=U,
              paired_vertices=collect(zip(left, clones)),
              fixed_vertices=shared, base_wedge_facet=f,
              base_facets_containing_face=containing,
              face_vertices_in_base=shared_in_base,
              face_dimension=face_dimension))
end

function _node(kind, name, d, I; base=nothing, face=nothing, witness=(;))
    fully_named = kind == :known ||
        (kind in (:pyramid, :dual) && base.fully_named) ||
        (kind == :wedge && base.fully_named && face.fully_named)
    PolytopeType(kind, name, d, copy(I), base, face, witness, fully_named)
end

_named_count(r) = (r.kind == :known ? 1 : 0) +
    (r.base === nothing ? 0 : _named_count(r.base)) +
    (r.wedged_face === nothing ? 0 : _named_count(r.wedged_face))

function _recognize(I, d, extra_known, extended_catalog, memo; try_dual=true)
    key = (d, _signature(I), try_dual)
    haskey(memo, key) && return memo[key]
    nf, nv = size(I)
    # Destack pyramids before the catalog, since small multiplexes are pyramids
    # too (M(3, 4) is the square pyramid). Simplices (nv == d + 1) stay named.
    if nv > d + 1
        # A facet with all but one vertex is exactly a pyramid base.
        for f in 1:nf
            sum(I[f, :]) == nv - 1 || continue
            U = findall(I[f, :])
            apex = only(findall(!, I[f, :]))
            base = _recognize(_face_incidence(I, U), d - 1,
                              extra_known, extended_catalog, memo)
            result = _node(:pyramid, "pyramid", d, I; base,
                           witness=(apex=apex, base_facet=f, base_vertices=U))
            memo[key] = result
            return result
        end
    end

    name = _known_name(I, d, extra_known; extended_catalog)
    if name !== nothing
        result = _node(:known, name, d, I)
        memo[key] = result
        return result
    end

    best = nothing
    # Try all candidate base-facet pairs until a fully named description exists.
    for a in 1:(nf - 1), b in (a + 1):nf
        inverse = _inverse_wedge(I, a, b)
        inverse === nothing && continue
        base = _recognize(inverse.base, d - 1,
                          extra_known, extended_catalog, memo)
        face = _recognize(inverse.face, inverse.face_dimension,
                          extra_known, extended_catalog, memo)
        operation = inverse.face_dimension == -1 ? "prism" :
            inverse.face_dimension == d - 2 ? "facet wedge" : "face wedge"
        candidate = _node(:wedge, operation, d, I;
                          base, face, witness=inverse.witness)
        if candidate.fully_named
            memo[key] = candidate
            return candidate
        end
        if best === nothing || _named_count(candidate) > _named_count(best)
            best = candidate
        end
    end
    # Fall back to the polar, which may not fall back to its own polar in turn.
    if try_dual
        polar = _recognize(BitMatrix(transpose(I)), d, extra_known,
                           extended_catalog, memo; try_dual=false)
        if polar.kind != :unknown
            candidate = _node(:dual, "dual", d, I; base=polar)
            if best === nothing || candidate.fully_named ||
                    _named_count(candidate) > _named_count(best)
                best = candidate
            end
        end
    end
    result = best === nothing ? _node(:unknown, "unrecognized", d, I) : best
    memo[key] = result
    result
end

function _extra_catalog(extra_known)
    [(String(name), _incidence(Q.pm_polytope), Int(Q.pm_polytope.COMBINATORIAL_DIM))
     for (name, Q) in extra_known]
end

"""
    determine_polytope(P; extra_known=[], extended_catalog=true)

Recognize the combinatorial type of a nonempty bounded Oscar polytope.
Priority: pyramid (unless a simplex) -> named catalog -> wedge over a proper
face -> dual -> unrecognized. Pyramids, wedges and duals are recognized
recursively; a dual is only used when no wedge is fully named, and the polar
itself is not dualized again. Wedges retain both the base
type and wedged-face type, together with the embedded face and vertex pairing.
Here a wedge means `{(x,t): x in Q, 0 <= t <= h(x)}`, where `h` is affine,
nonnegative on `Q`, and vanishes exactly on the marked face. The empty face
is allowed and gives a prism. Every decomposition is checked by reconstructing
the complete vertex-facet incidence.

The catalog is finite apart from parameterized families; unrecognized is not
a classification theorem. Supply additional named Oscar polytopes as
`extra_known=["my type" => Q]`. The main families are simplices, cubes,
cross-polytopes, multiplexes, products of simplices, cyclic polytopes, the
Platonic solids, prisms, bipyramids, antiprisms, trapezohedra, duoprisms and
duopyramids. `extended_catalog=false` skips the secondary catalog (Archimedean
and Catalan solids, Johnson solids, the elongated dodecahedron, the uniform
4-polytopes, hypersimplices, permutahedra, associahedra).
"""
function determine_polytope(P; extra_known=[], extended_catalog=true)
    pm = P.pm_polytope
    Bool(pm.FEASIBLE) && Bool(pm.BOUNDED) ||
        throw(ArgumentError("P must be a nonempty bounded polytope."))
    d = Int(pm.COMBINATORIAL_DIM)
    _recognize(_incidence(pm), d, _extra_catalog(extra_known),
               extended_catalog, Dict{Any,PolytopeType}())
end

"""
    determine_polytope(incidence, d; extra_known=[], extended_catalog=true)

As `determine_polytope(P)`, for a `d`-polytope given by its facet-by-vertex
incidence matrix.
"""
determine_polytope(incidence::AbstractMatrix{Bool}, d::Integer; extra_known=[], extended_catalog=true) =
    _recognize(BitMatrix(incidence), Int(d), _extra_catalog(extra_known),
               extended_catalog, Dict{Any,PolytopeType}())

"""
    determine_5_polytope(P; extra_known=[], extended_catalog=true)

As `determine_polytope`, requiring intrinsic dimension five. Uses
`P.pm_polytope` throughout; no `pm_object` calls.
"""
function determine_5_polytope(P; kwargs...)
    Int(P.pm_polytope.COMBINATORIAL_DIM) == 5 ||
        throw(ArgumentError("P must be 5-dimensional."))
    determine_polytope(P; kwargs...)
end

function describe_type(r::PolytopeType)
    r.kind == :pyramid && return "Pyramid(" * describe_type(r.base) * ")"
    r.kind == :dual && return "Dual(" * describe_type(r.base) * ")"
    r.kind == :wedge && r.wedged_face.dimension == -1 &&
        return "Prism(" * describe_type(r.base) * ")"
    r.kind == :wedge && return "Wedge(" * describe_type(r.base) *
        "; face=" * describe_type(r.wedged_face) * ")"
    r.kind == :unknown && return "Unrecognized($(r.dimension)D; " *
        "$(size(r.incidence, 2)) vertices, $(size(r.incidence, 1)) facets)"
    r.name
end

function print_type(r::PolytopeType; io=stdout, indent="", witnesses=false)
    println(io, indent, describe_type(r))
    if witnesses && !isempty(r.witness)
        println(io, indent, "  witness: ", r.witness)
    end
    if r.base !== nothing
        println(io, indent, r.kind == :dual ? "  polar:" : "  base:")
        print_type(r.base; io, indent=indent * "    ", witnesses)
    end
    if r.wedged_face !== nothing
        println(io, indent, "  wedged face:")
        print_type(r.wedged_face; io, indent=indent * "    ", witnesses)
    end
    nothing
end

"""Pivot columns of a row echelon form, by exact Gaussian elimination."""
function _pivot_columns(Y)
    A = copy(Y)
    pivots = Int[]
    for j in axes(A, 2)
        row = length(pivots) + 1
        row > size(A, 1) && break
        p = findfirst(!iszero, A[row:end, j])
        p === nothing && continue
        A[[row, row + p - 1], :] = A[[row + p - 1, row], :]
        for i in (row + 1):size(A, 1)
            A[i, :] -= A[i, j] // A[row, j] * A[row, :]
        end
        push!(pivots, j)
    end
    pivots
end

"""
Vertices of a polar of `conv(X)` (vertex rows ordered as the columns of `I`),
ordered as the facets of `I`. The polar is taken in a coordinate projection of
the linear span of `X`, centered at its vertex barycenter.
"""
function _polar_vertices(X, I)
    Y = X .- sum(X; dims=1) / size(X, 1)
    Z = Y[:, _pivot_columns(Y)]
    pm = Polymake.polytope.Polytope{Polymake.Rational}(POINTS=[ones(Int, size(Z, 1)) Z])
    F = Matrix{Rational{BigInt}}(pm.FACETS)
    polar = similar(Z, size(I, 1), size(Z, 2))
    for f in axes(F, 1)
        tight = [iszero(F[f, 1] + sum(F[f, 2:end] .* Z[v, :])) for v in axes(Z, 1)]
        polar[only(i for i in axes(I, 1) if I[i, :] == tight), :] = -F[f, 2:end] / F[f, 1]
    end
    polar
end

"""
Vertices of `conv(X)` in coordinates of its affine hull: unchanged if already
full-dimensional; if the vertices are integral, the lattice `aff ∩ Zⁿ` becomes
`Zᵏ` (polymake's `ambient_lattice_normalization`); otherwise a coordinate
projection, which keeps the combinatorial type but not the lattice.
"""
function _full_dimensional(X)
    pivots = _pivot_columns(X .- X[1:1, :])
    length(pivots) == size(X, 2) && return X
    all(isinteger, X) || return X[:, pivots]
    pm = Polymake.polytope.Polytope{Polymake.Rational}(POINTS=[ones(Int, size(X, 1)) X])
    V = Matrix{Rational{BigInt}}(Polymake.polytope.ambient_lattice_normalization(pm).VERTICES)
    V[:, 2:end] ./ V[:, 1]
end

function _unrecognized_parts!(parts, r, X, path)
    if r.kind == :unknown
        push!(parts, join(path, "_") => _full_dimensional(X))
    elseif r.kind == :pyramid
        _unrecognized_parts!(parts, r.base, X[r.witness.base_vertices, :], [path; "base"])
    elseif r.kind == :wedge
        U = r.witness.base_vertices
        _unrecognized_parts!(parts, r.base, X[U, :], [path; "base"])
        _unrecognized_parts!(parts, r.wedged_face,
                             X[U[r.witness.face_vertices_in_base], :], [path; "face"])
    elseif r.kind == :dual
        _unrecognized_parts!(parts, r.base, _polar_vertices(X, r.incidence), [path; "polar"])
    end
    parts
end

"""
    unrecognized_parts(P, r)

Rational vertex coordinates (one vertex per row) of every unrecognized node of
`r = determine_polytope(P)`, as `path => vertices` pairs. The `path` names the
position below the root, such as `"base"` or `"polar_base_face"`, and is empty
for the root. Pyramid bases and wedge pieces are faces of their parent, a polar
is realized from its parent's facets. Each part is returned full-dimensional,
in coordinates of its affine hull, lattice-preserving when its vertices are
integral; the vertex order is then not that of the node's incidence.
"""
function unrecognized_parts(P, r::PolytopeType)
    V = Matrix{Rational{BigInt}}(P.pm_polytope.VERTICES)
    _unrecognized_parts!(Pair{String,Matrix{Rational{BigInt}}}[], r,
                         V[:, 2:end] ./ V[:, 1], String[])
end

end # module
