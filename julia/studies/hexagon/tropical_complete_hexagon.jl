using Oscar
using Random

# Completing a tropical hexagon: the tropicalization trop(σ) of the section σ of Section 9 of main.tex. Min convention.
#
# Given tropical flags (x^(1), p^(12)) and (x^(4), p^(45)), every coordinate f of the polynomial functions of Lemma 9.1
# (polynomial_expressions.jl) is replaced by a tropically faithful representative at the weight w = (p^(12), x^(4)),
# respectively w = (p^(45), x^(1)), modulo the ideal I of linespace × pointspace. Its tropical value at w is the valuation
# of f at a generic lift. By default the representative is looked up in the certificate certificate/hexagon_certificate.mrdi
# (`certificate_representative`); with `method=:algorithm` it is computed at w by Algorithm 2
# (tropically_faithful_representatives.jl).
#
# Steps 1-3 of Algorithm 2 only depend on the Gröbner cone of w, so they are done once per Gröbner cone. These cones are
# the products of a face of the raw fan of trop(linespace), data/g2/grassmannian/01_raw_fan.mrdi, which is
# the Gröbner fan restricted to the tropical variety, and a cone of the tropical quadric. The Gröbner bases are
# computed for a small representative weight of the cone, as Singular cannot handle large weights.
#
# The tropical variety depends on the characteristic (Section 4); all of this is over QQ.

include(joinpath(@__DIR__, "complete_hexagon.jl"))
include(joinpath(@__DIR__, "tropically_faithful_representatives.jl"))

const RAW_FAN_PATH = joinpath(@__DIR__, "..", "..", "data", "g2", "grassmannian", "01_raw_fan.mrdi")
const FLAG_FAN_PATH = joinpath(@__DIR__, "..", "..", "data", "g2", "flag", "03_normalized_coarsest_fan.mrdi")

# --- The Gröbner cones ------------------------------------------------------------------------------------

const _RAW_FAN_FACES = Ref{Any}(nothing)

"""
The faces of the raw fan of trop(linespace): for each a small integer weight `point` in its relative interior (the sum
of its rays), and rational matrices `E`, `A` with `E*w = 0`, `A*w < 0` exactly on the relative interior. Cached.
"""
function raw_fan_faces()
    if _RAW_FAN_FACES[] === nothing
        F = load(RAW_FAN_PATH)
        lineality = matrix(QQ, reduce(vcat, [permutedims(integer_weight(l)) for l in lineality_space(F)]))
        L = positive_hull(zeros(QQ, 0, 21), lineality)
        faces = Any[(point=zeros(Int, 21), E=linear_equation_matrix(linear_span(L)), A=zero_matrix(QQ, 0, 21))]
        for d in (lineality_dim(F) + 1):dim(F), C in cones(F, d)
            push!(faces, (point=sum(integer_weight(r) for r in rays_modulo_lineality(C)[1]),
                          E=linear_equation_matrix(linear_span(C)), A=linear_inequality_matrix(facets(C))))
        end
        _RAW_FAN_FACES[] = faces
    end
    return _RAW_FAN_FACES[]
end

"The index of the face of the raw fan whose relative interior contains the line weight `w`, or `nothing`."
function locate_face(w)
    v = matrix(QQ, 21, 1, QQ.(w))
    return findfirst(f -> all(iszero, f.E * v) && all(<(0), f.A * v), raw_fan_faces())
end

"The cone of the tropical quadric containing the point weight `x`: the terms of q of minimal weight, among x1x4, x2x5, x3x6, x7^2."
function quadric_cone(x)
    weights = [x[1] + x[4], x[2] + x[5], x[3] + x[6], 2 * x[7]]
    return findall(==(minimum(weights)), weights)
end

"A weight in the relative interior of the quadric cone `S`: the terms in `S` have weight 0, the others weight 1 or 2."
function quadric_cone_point(S)
    w = zeros(Int, 7)
    for k in setdiff(1:4, S)
        k <= 3 ? (w[k + 3] = 1) : (w[7] = 1)
    end
    return w
end

const _CONE_SETUPS = Dict{Any,FaithfulSetup}()

"""
    cone_setup(w)

The `FaithfulSetup` of the Gröbner cone containing `w = vcat(p, x)`, with `w` a point of trop(linespace) × trop(pointspace).
Computed once per Gröbner cone, which takes some seconds. The ideal of Remark 9.2 is prime, hence saturated.
"""
function cone_setup(w)
    k = locate_face(w[1:21])
    k === nothing && throw(ArgumentError("the line $(w[1:21]) is not a point of the tropical G2-Grassmannian"))
    S = quadric_cone(w[22:28])
    length(S) >= 2 || throw(ArgumentError("the point $(w[22:28]) is not a point of the tropical quadric"))
    return get!(_CONE_SETUPS, (k, S)) do
        faithful_setup(line_point_ideal(), vcat(raw_fan_faces()[k].point, quadric_cone_point(S)); saturate=false)
    end
end

# --- Tropically faithful representatives from the certificate ---------------------------------------------

const HEXAGON_CERTIFICATE_PATH = joinpath(@__DIR__, "certificate", "hexagon_certificate.mrdi")
const _HEXAGON_CERTIFICATE = Ref{Any}(nothing)

"A polynomial with the variables of `EXPRESSION_RING` as an element of it (each `load` creates a new ring)."
function _certificate_polynomial(p)
    B = MPolyBuildCtx(EXPRESSION_RING)
    foreach((c, e) -> push_term!(B, QQ(c), e), AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
    return finish(B)
end

"""
The certificate certificate/hexagon_certificate.mrdi (see certificate/README.md), loaded on first use: for each of the 56
coordinates its representatives h_j = x^(-u_j) g_j, with the sorted exponent vectors of g_j as the rows of a matrix, and
for each maximal cone of Trop the 56 trees, parsed when first used.
"""
function hexagon_certificate()
    if _HEXAGON_CERTIFICATE[] === nothing
        d = load(HEXAGON_CERTIFICATE_PATH)
        reps = map(zip(d["g"], d["u"])) do (gs, us)
            map(zip(gs, us)) do (g0, u)
                g = _certificate_polynomial(g0)
                (g=g, u=Int.(u), exponents=reduce(vcat, permutedims.(sort(collect(AbstractAlgebra.exponent_vectors(g))))))
            end
        end
        cones = [(k=Int(k), S=[parse(Int, ch) for ch in S]) for (k, S) in zip(d["face"], d["quadric"])]
        _HEXAGON_CERTIFICATE[] = (names=collect(String, d["coordinate"]), reps=reps, cones=cones,
                                  trees=[collect(String, t) for t in d["tree"]], parsed=Dict{Tuple{Int,Int},Any}())
    end
    return _HEXAGON_CERTIFICATE[]
end

"A tree of the certificate as nested named tuples (j, children), children a vector of (t, node), a leaf being `nothing`."
function parse_certificate_tree(s)
    pos = Ref(1)
    function number()
        start = pos[]
        while pos[] <= ncodeunits(s) && isdigit(s[pos[]])
            pos[] += 1
        end
        return parse(Int, s[start:pos[] - 1])
    end
    function node()
        s[pos[]] == '.' && (pos[] += 1; return nothing)
        j = number()
        pos[] += 1                                  # '('
        children = Tuple{Int,Any}[]
        while true
            t = number()
            pos[] += 1                              # ':'
            push!(children, (t, node()))
            s[pos[]] == ',' || break
            pos[] += 1
        end
        pos[] += 1                                  # ')'
        return (j=j, children=children)
    end
    return node()
end

"The index in `hexagon_certificate().cones` of a maximal cone of Trop containing `w = vcat(p, x)`."
function certificate_cone(w)
    C = hexagon_certificate()
    v = matrix(QQ, 21, 1, QQ.(w[1:21]))
    faces = raw_fan_faces()
    in_face = Dict{Int,Bool}()
    in_closed_face(k) = get!(() -> all(iszero, faces[k].E * v) && all(<=(0), faces[k].A * v), in_face, k)
    minimal_terms = quadric_cone(w[22:28])
    n = findfirst(c -> issubset(c.S, minimal_terms) && in_closed_face(c.k), C.cones)
    n === nothing && throw(ArgumentError("$w is not a point of trop(linespace) × trop(pointspace)"))
    return n
end

"""
    certificate_representative(i, w)

The tropically faithful representative f_σ = x^(-u) g of the coordinate `i` (1..56: x2_1..7, x3_1..7, p23_1..21,
p34_1..21) at a point `w = vcat(p, x)` of Trop, from the certificate: in the tree of a maximal cone containing `w`, the
child of each node is one whose vertex v attains min <e, w> over the exponents e of g_j, i.e. whose cell σ ∩ N_v contains
`w`. Returns `(j, g, u, value)` with `value` = trop(f_σ)(w) = val(f(z)) for a generic lift z of `w`. On a wall several
cells contain `w`; their representatives are all faithful at `w`.
"""
function certificate_representative(i, w)
    C = hexagon_certificate()
    n = certificate_cone(w)
    node = get!(() -> parse_certificate_tree(C.trees[n][i]), C.parsed, (n, i))
    q = QQ.(w)
    den = reduce(lcm, denominator.(q); init=ZZ(1))
    wi = [BigInt(numerator(c * den)) for c in q]
    while true
        r = C.reps[i][node.j]
        at = r.exponents * wi
        m = minimum(at)
        k = findfirst(((t, _),) -> at[t] == m, node.children)
        k === nothing && error("coordinate $i, cone $n: no cell of the tree contains $w")
        child = node.children[k][2]
        child === nothing && return (j=node.j, g=r.g, u=r.u, value=QQ(m - sum(r.u .* wi)) // den)
        node = child
    end
end

# --- Tropical hexagons ------------------------------------------------------------------------------------

"""
    tropical_path(p, x; method=:certificate)

The tropicalization of Lemma lem:complete_hexagon at a tropical line `p` and a tropical point `x`: the values of tropically
faithful representatives of all coordinates. Returns `(x2, x3, p23, p34, no_representative)`.

With `method=:certificate` the representatives are those of certificate/hexagon_certificate.mrdi
(`certificate_representative`), which are faithful on all of Trop, so `no_representative` is empty. With
`method=:algorithm` they are computed at `(p, x)` by the algorithm of tropically_faithful_representatives.jl;
`no_representative` lists the coordinates, as `(name, c)`, for which it finds that none exists. For these the value is a
lower bound of the valuation at a lift, and the valuation at a generic lift depends on the component of V(in_w I) the
lift reduces to.
"""
function tropical_path(p, x; method::Symbol=:certificate)
    w = vcat(QQ.(p), QQ.(x))
    if method == :certificate
        values = [certificate_representative(i, w).value for i in 1:56]
        return (x2=values[1:7], x3=values[8:14], p23=values[15:35], p34=values[36:56], no_representative=Tuple{Symbol,Int}[])
    end
    method == :algorithm || throw(ArgumentError("method must be :certificate or :algorithm"))
    T = cone_setup(w)
    H = hexagon_expressions()
    no_representative = Tuple{Symbol,Int}[]
    function values(name)
        return map(enumerate(H[name])) do (c, f)
            r = tropically_faithful_representative(f, T; w)
            r.status == :no_representative && push!(no_representative, (name, c))
            r.status == :in_ideal && error("coordinate $c of $name vanishes on linespace × pointspace")
            return r.value
        end
    end
    x2, x3, p23, p34 = values(:x2), values(:x3), values(:p23), values(:p34)
    return (x2=x2, x3=x3, p23=p23, p34=p34, no_representative=no_representative)
end

"""
Whether (x, p) satisfies the relations of the flag variety tropically: the minimum is attained twice in the quadric, the
7 linear forms, and the Plücker and incidence relations. This is a necessary condition only.
"""
function is_tropical_flag_candidate(x, p)
    nu = tropical_semiring_map(QQ)
    w = vcat(QQ.(p), QQ.(x))
    relations = vcat(linear_forms(P_VARS), plucker_relations(P_VARS), [quadric(X_VARS)], incidence_relations(X_VARS, P_VARS))
    return all(g -> length(initial(g, nu, w)) >= 2, relations)
end

"A tropical vector shifted to minimum 0."
shift_to_zero(v) = v .- minimum(v)

"""
    tropical_complete_hexagon(x1, p12, x4, p45; check=true, normalize=true, method=:certificate)
    tropical_complete_hexagon((x1, p12), (x4, p45); ...)

The tropical hexagon through the tropical flags (x^(1), p^(12)) and (x^(4), p^(45)): the valuations of the hexagon of a
generic lift. Input vectors have rational entries, the output has the fields of `complete_hexagon`, with vectors over
QQ, and `no_representative`: the coordinates of each of the two completions without a faithful representative, see
`tropical_path` (always empty for `method=:certificate`, the default). With `normalize` every vector is shifted to
minimum 0. With `check` the inputs are tested by `is_tropical_flag_candidate`.
"""
function tropical_complete_hexagon(x1, p12, x4, p45; check::Bool=true, normalize::Bool=true, method::Symbol=:certificate)
    x1, p12, x4, p45 = QQ.(x1), QQ.(p12), QQ.(x4), QQ.(p45)
    if check
        is_tropical_flag_candidate(x1, p12) || throw(ArgumentError("(x1, p12) does not satisfy the tropical flag relations"))
        is_tropical_flag_candidate(x4, p45) || throw(ArgumentError("(x4, p45) does not satisfy the tropical flag relations"))
    end
    A = tropical_path(p12, x4; method)
    B = tropical_path(p45, x1; method)
    points = [x1, A.x2, A.x3, x4, B.x2, B.x3]
    lines = [p12, A.p23, A.p34, p45, B.p23, B.p34]
    if normalize
        points, lines = shift_to_zero.(points), shift_to_zero.(lines)
    end
    return (points=points, lines=lines, flags=[(points[k], lines[k]) for k in 1:6],
            no_representative=(A.no_representative, B.no_representative))
end

tropical_complete_hexagon(f1::Tuple, f4::Tuple; kwargs...) = tropical_complete_hexagon(f1..., f4...; kwargs...)

# --- From classical to tropical ---------------------------------------------------------------------------

"The t-adic valuations of a vector over K = QQ(t) with nonzero entries."
function tropicalize(v)
    nu = tropical_semiring_map(parent(v[1]), gen(parent(v[1])))
    return [QQ(nu(c)) for c in v]
end

# --- Random tropical flags ---------------------------------------------------------------------------------

const _TROPICAL_FLAG_FAN = Ref{Any}(nothing)

"The fan of trop(flagspace) in R^7 x R^21 from `FLAG_FAN_PATH`. Cached."
tropical_flag_fan() = _TROPICAL_FLAG_FAN[] === nothing ? (_TROPICAL_FLAG_FAN[] = load(FLAG_FAN_PATH)) : _TROPICAL_FLAG_FAN[]

"""
    random_tropical_flag(; cone=nothing, coeff=1000, rng=Random.default_rng())

A random point `(x, p)` of trop(flagspace): a random positive combination, with integer coefficients in `1:coeff`, of
the rays of a maximal cone of the fan, chosen uniformly or given by its number `cone`.
"""
function random_tropical_flag(; cone=nothing, coeff::Integer=1000, rng=Random.default_rng())
    C = maximal_cones(tropical_flag_fan())
    i = cone === nothing ? rand(rng, 1:length(C)) : cone
    v = sum(rand(rng, 1:coeff) * QQFieldElem[c for c in r] for r in rays_modulo_lineality(C[i])[1])
    return (v[1:7], v[8:28])
end
