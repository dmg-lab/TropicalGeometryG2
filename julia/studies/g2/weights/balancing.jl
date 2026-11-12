using Oscar

# Balancing check for the multiplicities in weights.csv (computed by weights.jl).
#
# For a wall τ (codimension-one cone of the raw fan) and a maximal cone σ ⊃ τ, let L_σ = Z^21 ∩ span(σ) and
# L_τ = Z^21 ∩ span(τ). Then L_σ/L_τ ≅ Z, and u_{σ/τ} is its generator pointing into σ. The fan is balanced at τ
# with weights m if Σ_{σ ⊃ τ} m(σ) u_{σ/τ} ∈ span(τ) (MS15, Definition 3.3.1). This is checked for both columns
# of weights.csv: mult_QQ and mult_QQbar.
#
# From the `julia` folder:  julia --project=. studies/g2/weights/balancing.jl

include(joinpath(@__DIR__, "weights.jl"))

"Integer row vectors spanning the cone (its primitive rays and the lineality basis)."
cone_generators(C, lineality) = vcat([integer_ray(Vector(r)) for r in rays_modulo_lineality(C)[1]], lineality)

"Row basis of the saturated lattice Z^21 ∩ span(rows)."
saturated_lattice(rows) = saturate(hnf(matrix(ZZ, reduce(vcat, permutedims.(rows)))))

"""
The lattice normal vector u_{σ/τ} modulo span(τ), as the rational vector g·r, where r is a ray of σ outside τ and
g generates the image of L_σ under the coordinate along r with respect to span(τ).
"""
function lattice_normal(sigma_gens, tau_gens, r)
    Lsigma = saturated_lattice(sigma_gens)
    T = matrix(QQ, reduce(vcat, permutedims.(vcat([r], tau_gens))))
    coords = QQFieldElem[]
    for i in 1:nrows(Lsigma)
        b = matrix(QQ, 1, ncols(Lsigma), [QQ(Lsigma[i, j]) for j in 1:ncols(Lsigma)])
        ok, x = can_solve_with_solution(T, b; side=:left)
        @assert ok "lattice vector not in span(σ)"
        push!(coords, x[1, 1])
    end
    g = reduce(gcd, coords)
    return g .* QQ.(r)
end

"Read the multiplicities (mult_QQ, mult_QQbar) from weights.csv, indexed by maximal cone."
function read_weights(path=WEIGHTS_PATH)
    rows = [split(l, ",") for l in Iterators.drop(eachline(path), 1) if !isempty(l)]
    mQQ = Dict(parse(Int, r[1]) => parse(Int, r[6]) for r in rows)
    mbar = Dict(parse(Int, r[1]) => parse(Int, r[8]) for r in rows if !isempty(r[8]))
    return mQQ, mbar
end

"""
    check_balancing()

For every wall of the raw fan: the maximal cones containing it, and whether the fan is balanced there with the
QQ and the QQbar multiplicities. Returns the walls where balancing fails, for each weighting.
"""
function check_balancing()
    F = load(RAW_FAN_PATH)
    lineality = [integer_ray(Vector(l)) for l in lineality_space(F)]
    maxcones = collect(maximal_cones(F))
    ray_sets(C) = Set(map(r -> integer_ray(Vector(r)), rays_modulo_lineality(C)[1]))
    max_rays = ray_sets.(maxcones)
    mQQ, mbar = read_weights()
    walls = collect(cones(F, dim(F) - 1))
    fail_QQ, fail_bar = Int[], Int[]
    for (k, tau) in enumerate(walls)
        tau_rays = ray_sets(tau)
        tau_gens = vcat(collect(tau_rays), lineality)
        adjacent = [i for i in eachindex(maxcones) if issubset(tau_rays, max_rays[i])]
        normals = [lattice_normal(cone_generators(maxcones[i], lineality), tau_gens,
                                  only(setdiff(max_rays[i], tau_rays))) for i in adjacent]
        span_tau = matrix(QQ, reduce(vcat, permutedims.(tau_gens)))
        in_span(v) = can_solve_with_solution(span_tau, matrix(QQ, 1, length(v), v); side=:left)[1]
        in_span(sum(mQQ[i] .* normals[n] for (n, i) in enumerate(adjacent))) || push!(fail_QQ, k)
        isempty(mbar) || in_span(sum(mbar[i] .* normals[n] for (n, i) in enumerate(adjacent))) || push!(fail_bar, k)
    end
    println("$(length(walls)) walls; balancing fails with QQ multiplicities at $(length(fail_QQ)) walls, ",
            "with QQbar multiplicities at $(length(fail_bar)) walls")
    for (name, fails) in (("QQ", fail_QQ), ("QQbar", fail_bar)), k in fails
        tau_rays = ray_sets(walls[k])
        adjacent = [i for i in eachindex(maxcones) if issubset(tau_rays, max_rays[i])]
        println("  $name fails at wall $k, adjacent maximal cones ", adjacent)
    end
    return (QQ=fail_QQ, QQbar=fail_bar)
end

if abspath(PROGRAM_FILE) == @__FILE__
    check_balancing()
end
