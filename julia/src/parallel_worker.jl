using Oscar
using LinearAlgebra


const PARALLEL_TROP_STATE = Ref{Any}(nothing)
const PARALLEL_TROP_SYMMETRY = Ref{Any}(nothing)
const PARALLEL_TROP_PURITY = Ref{Any}(nothing)

"The degree-first weighted order refines weights only on a homogeneous ideal."
function _parallel_trop_validate_homogeneous_state(state)
    for (source, polynomials) in (("defining ideal", gens(state.ideal)),
                                  ("exported generators", state.generators))
        for polynomial in polynomials
            iszero(polynomial) && continue
            all(sum(exponent) == total_degree(polynomial) for exponent in exponents(polynomial)) ||
                error("parallel tropicalization requires homogeneous $source")
        end
    end
    return state
end

"""
Complete a provider's `(ring, ideal, generators)` with its torus coordinates.
`variable_names` are always the ring's names. A provider working in a reduced
torus chart also returns `original_variable_names` and `lift_indices`, the
index of the ring coordinate lifting each original coordinate; otherwise the
original coordinates are the ring's own. A provider may declare the Krull
`dimension` of its ideal when computing it is out of reach.
"""
function _parallel_trop_provider_state(data)
    for field in (:ring, :ideal, :generators)
        hasproperty(data, field) || error("ideal provider must return field `$field`")
    end
    variable_names = string.(symbols(data.ring))
    hasproperty(data, :variable_names) && data.variable_names != variable_names &&
        error("provider variable_names differ from the names of its ring")
    reduced = hasproperty(data, :lift_indices)
    original = hasproperty(data, :original_variable_names) ?
        collect(String, data.original_variable_names) : copy(variable_names)
    reduced || original == variable_names ||
        error("a provider in a reduced torus chart must also return lift_indices")
    lift_indices = reduced ? collect(Int, data.lift_indices) : collect(eachindex(variable_names))
    length(lift_indices) == length(original) ||
        error("lift_indices must have one entry per original coordinate")
    all(in(eachindex(variable_names)), lift_indices) ||
        error("lift_indices must index the ring coordinates")
    !hasproperty(data, :dimension) || data.dimension isa Integer && data.dimension >= 0 ||
        error("a declared dimension must be a nonnegative integer")
    return merge(data, (variable_names=variable_names, original_variable_names=original,
                        lift_indices=lift_indices))
end

"The sign of each image relative to its expected variable."
function _parallel_trop_signs(image, variables)
    length(image) == length(variables) || error("signed action is not on the full coordinate ring")
    return [image[k] == variables[k] ? 1 : image[k] == -variables[k] ? -1 :
            error("signed coordinate action does not match its permutation") for k in eachindex(image)]
end

"Signed images of the variables under the group generated gfan-style by `generators`."
function _parallel_trop_variable_group(state, generators)
    variables = gens(state.ring)
    group = Main.TropicalCensus.signed_permutation_group(generators, length(variables))
    images = Main.TropicalCensus.signed_coordinate_images(group)
    return group.permutations, [images(p, variables) for p in group.permutations]
end

"""Load an ideal provider and initialize this worker's private Oscar state."""
function initialize_parallel_trop_worker(
    provider_file::String,
    provider_module_name::Symbol,
    builder_name::Symbol,
    candidates_name::Union{Nothing, Symbol}=nothing,
    coordinate_images_name::Union{Nothing, Symbol}=nothing,
    certification_strategy::Symbol=:groebner,
    purity_verifier_name::Union{Nothing, Symbol}=nothing,
    symmetry_action::Symbol=:taxa,
)
    isdefined(Main, provider_module_name) || Base.include(Main, provider_file)
    provider = Base.invokelatest(getfield, Main, provider_module_name)
    builder = Base.invokelatest(getfield, provider, builder_name)
    PARALLEL_TROP_STATE[] = _parallel_trop_validate_homogeneous_state(
        _parallel_trop_provider_state(Base.invokelatest(builder)))
    PARALLEL_TROP_SYMMETRY[] = nothing
    PARALLEL_TROP_PURITY[] = nothing
    if !isnothing(candidates_name)
        symmetry_action in (:taxa, :variables) || error("unknown symmetry action $symmetry_action")
        generators = Base.invokelatest(getfield, provider, candidates_name)
        state = PARALLEL_TROP_STATE[]
        if symmetry_action == :variables
            isnothing(coordinate_images_name) ||
                error("variable symmetries take their signs as generators `p => signs`")
            candidates, images = _parallel_trop_variable_group(state, generators)
            Main.TropicalCensus.validate_tropical_permutation_group(candidates)
        else
            isnothing(coordinate_images_name) && error("taxon symmetry requires a signed coordinate action")
            candidates = [Int.(collect(permutation)) for permutation in generators]
            coordinate_images = Base.invokelatest(getfield, provider, coordinate_images_name)
            Main.TropicalCensus.validate_tropical_permutation_group(candidates)
            images = [Base.invokelatest(coordinate_images, permutation, gens(state.ring))
                      for permutation in candidates]
            for (permutation, image) in zip(candidates, images)
                indices = Main.TropicalCensus.tropical_coordinate_permutation(invperm(permutation), ngens(state.ring))
                _parallel_trop_signs(image, gens(state.ring)[indices])
            end
        end
        actions = [hom(state.ring, state.ring, image) for image in images]
        normalized = Set(_parallel_trop_normalize_witness.(state.generators))
        for action in actions
            transformed = action.(state.generators)
            Set(_parallel_trop_normalize_witness.(transformed)) == normalized ||
                ideal(state.ring, transformed) == state.ideal ||
                error("supplied signed $(symmetry_action == :taxa ? "taxon" : "variable") action does not preserve the defining ideal")
        end
        PARALLEL_TROP_SYMMETRY[] = (permutations=candidates, actions=actions, action=symmetry_action)
    end
    certification_strategy in (:groebner, :balanced_support) || error("unknown tropical certification strategy")
    if certification_strategy == :balanced_support
        isnothing(PARALLEL_TROP_SYMMETRY[]) && error("balanced support requires verified signed symmetry")
        isnothing(purity_verifier_name) && error("balanced support requires an exact provider purity verifier")
        verifier = Base.invokelatest(getfield, provider, purity_verifier_name)
        certificate = Base.invokelatest(verifier, PARALLEL_TROP_STATE[])
        Main.TropicalCensus.validate_support_purity_certificate(certificate, dim(PARALLEL_TROP_STATE[].ideal))
        PARALLEL_TROP_PURITY[] = certificate
    end
    return nothing
end


function parallel_trop_state()
    isnothing(PARALLEL_TROP_STATE[]) && error(
        "parallel tropical worker has not been initialized",
    )
    return PARALLEL_TROP_STATE[]
end


"""Render an Oscar polynomial in the syntax accepted by gfan."""
function parallel_trop_gfan_polynomial(f)
    iszero(f) && return "0"
    names = string.(symbols(parent(f)))
    coefficient_field = coefficient_ring(f)
    finite_coefficients = coefficient_field isa FinField
    finite_coefficients && degree(coefficient_field) != 1 &&
        throw(ArgumentError("gfan support export requires a prime coefficient field"))
    pieces = String[]

    for (coefficient, exponent) in zip(coefficients(f), exponents(f))
        # gfan computes the trivial-valuation prevariety from polynomial
        # supports. Integer lifts preserve those supports; every algebraic
        # check and round-trip parse still uses the worker's coefficient field.
        representative = finite_coefficients ? lift(ZZ, coefficient) : coefficient
        negative = representative < 0
        absolute = abs(representative)
        numerator_string = string(numerator(absolute))
        denominator_value = denominator(absolute)
        coefficient_string = isone(denominator_value) ?
            numerator_string : numerator_string * "/" * string(denominator_value)

        factors = String[]
        for (name, power) in zip(names, exponent)
            iszero(power) && continue
            push!(factors, power == 1 ? name : name * "^" * string(power))
        end
        monomial = join(factors, "*")
        body = if isempty(monomial)
            coefficient_string
        elseif isone(absolute)
            monomial
        else
            coefficient_string * "*" * monomial
        end

        if isempty(pieces)
            push!(pieces, negative ? "-" * body : body)
        else
            push!(pieces, negative ? " - " * body : " + " * body)
        end
    end
    return join(pieces)
end


function parallel_trop_export_data()
    state = parallel_trop_state()
    return (
        variable_names=state.variable_names,
        original_variable_names=state.original_variable_names,
        lift_indices=state.lift_indices,
        generators=parallel_trop_gfan_polynomial.(state.generators),
        # Under :groebner certification the dimension only orders the checks;
        # balanced support compares it with an exact purity certificate.
        ideal_dimension=hasproperty(state, :dimension) ? state.dimension : dim(state.ideal),
        declared_dimension=hasproperty(state, :dimension),
        homogeneous=true,
        symmetry_permutations=isnothing(PARALLEL_TROP_SYMMETRY[]) ? nothing :
            PARALLEL_TROP_SYMMETRY[].permutations,
        symmetry_action=isnothing(PARALLEL_TROP_SYMMETRY[]) ? nothing :
            PARALLEL_TROP_SYMMETRY[].action,
        purity_certificate=PARALLEL_TROP_PURITY[],
    )
end

"""Exact torus nonemptiness, with cyclic variable saturation for homogeneous ideals."""
function _parallel_trop_monomial_free(initial_ideal)
    ring = base_ring(initial_ideal)
    polynomials = collect(gens(initial_ideal))
    homogeneous = all(polynomial -> iszero(polynomial) ||
        all(sum(exponent) == total_degree(polynomial) for exponent in exponents(polynomial)),
        polynomials)
    if homogeneous && isdefined(Oscar, :GITFans) &&
            isdefined(Oscar.GITFans, :is_monomial_free)
        # A single satstd call only partially saturates. The GITFans routine
        # cycles every variable into the last position in a dp ring, where
        # homogeneity makes saturation complete (BKR20, Prop. 3.1,
        # Remark 3.2 and Algorithm 3.3). Together the cycles compute
        # saturation by the product of all coordinates exactly.
        # Use a fresh, explicit dp ring so no weighted-order cache or custom
        # default order can invalidate the last-variable argument.
        dp_ring, dp_variables = polynomial_ring(base_ring(ring), string.(symbols(ring));
            internal_ordering=:degrevlex, cached=false)
        default_ordering(dp_ring) == degrevlex(dp_ring) ||
            error("cyclic variable saturation requires degree-reverse-lexicographic ordering")
        to_dp = hom(ring, dp_ring, dp_variables)
        dp_ideal = ideal(dp_ring, to_dp.(polynomials))
        return (monomial_free=Oscar.GITFans.is_monomial_free(dp_ideal), method=:cyclic_satstd)
    end
    return (monomial_free=!isone(saturation(initial_ideal, ideal(ring, prod(gens(ring))))),
            method=:product_saturation)
end

"Exact interior membership without constructing an unnecessary excluded-point witness."
function parallel_trop_support_membership(job)
    start_time = time()
    state = parallel_trop_state()
    length(job.weight) == ngens(state.ring) || throw(ArgumentError("weight has the wrong length"))
    basis = _parallel_trop_weighted_groebner_basis(state, job.weight)
    basis_seconds = time() - start_time
    empty_directions = Vector{Vector{Int}}()
    direct = _parallel_trop_select_witness(basis,
        (weight=job.weight, rays=empty_directions, lineality=empty_directions))
    saturation_seconds = 0.0
    method = :direct_monomial
    tropical = if !isnothing(direct)
        false
    else
        initial_forms = _parallel_trop_normalized_initial_form.(basis, Ref(job.weight))
        if all(polynomial -> length(polynomial) == 2, initial_forms)
            # Weight-initial forms of a weight-refining Groebner basis are a
            # Groebner basis of the weight-initial ideal. Dividing a monomial
            # by a binomial basis replaces it by one nonzero monomial at each
            # step. A monomial therefore cannot reduce to zero when this basis
            # contains no monomial, proving torus nonemptiness exactly.
            method = :binomial_groebner_basis
            true
        else
            initial_ideal = ideal(state.ring, initial_forms)
            saturation_started = time()
            membership = _parallel_trop_monomial_free(initial_ideal)
            saturation_seconds = time() - saturation_started
            method = membership.method
            membership.monomial_free
        end
    end
    return (valid=tropical, reason=tropical ? :tropical_interior : :monomial_contained,
            basis_size=length(basis), seconds=time() - start_time,
            basis_seconds=basis_seconds, saturation_seconds=saturation_seconds,
            membership_method=method)
end

"Parse only the arithmetic syntax emitted by parallel_trop_gfan_polynomial."
function _parallel_trop_parse_polynomial(text::AbstractString, ring)
    variables = Dict(Symbol(name) => variable for (name, variable) in
                     zip(string.(symbols(ring)), gens(ring)))
    function convert(expression)
        expression isa Integer && return ring(expression)
        expression isa Symbol && return variables[expression]
        expression isa Expr && expression.head == :call || error("unsupported polynomial syntax")
        operation, arguments = expression.args[1], expression.args[2:end]
        operation == :+ && return sum(convert.(arguments); init=zero(ring))
        operation == :* && return prod(convert.(arguments); init=one(ring))
        operation == :- && return length(arguments) == 1 ? -convert(only(arguments)) :
            convert(arguments[1]) - convert(arguments[2])
        operation == :/ && length(arguments) == 2 && arguments[2] isa Integer &&
            return convert(arguments[1]) / QQ(arguments[2])
        operation == :^ && length(arguments) == 2 && arguments[2] isa Integer &&
            arguments[2] >= 0 && return convert(arguments[1]) ^ arguments[2]
        error("unsupported polynomial operation $operation")
    end
    return convert(Meta.parse(text))
end

_parallel_trop_support_key(polynomial) = Tuple(sort!(Tuple.(collect(exponents(polynomial)))))

function _parallel_trop_polynomial_orbit(polynomial, actions)
    representatives = Dict{Tuple, String}()
    for action in actions
        image = _parallel_trop_normalize_witness(action(polynomial))
        key = _parallel_trop_support_key(image)
        text = parallel_trop_gfan_polynomial(image)
        representatives[key] = min(get(representatives, key, text), text)
    end
    return sort!(collect(values(representatives)))
end

function parallel_trop_symmetrize_polynomials(polynomials::AbstractVector{<:AbstractString})
    state = parallel_trop_state()
    symmetry = PARALLEL_TROP_SYMMETRY[]
    isnothing(symmetry) && return polynomials
    result = Dict{Tuple, String}()
    for text in polynomials
        polynomial = _parallel_trop_parse_polynomial(text, state.ring)
        for text in _parallel_trop_polynomial_orbit(polynomial, symmetry.actions)
            image = _parallel_trop_parse_polynomial(text, state.ring)
            key = _parallel_trop_support_key(image)
            result[key] = min(get(result, key, text), text)
        end
    end
    return sort!(collect(values(result)))
end

function parallel_trop_basis_support_signature(polynomials::AbstractVector{<:AbstractString})
    ring = parallel_trop_state().ring
    supports = sort!(unique(string(_parallel_trop_support_key(
        _parallel_trop_parse_polynomial(text, ring),
    )) for text in polynomials))
    return join(vcat(string.(symbols(ring)), supports), '\n')
end

function _parallel_trop_witness_orbit(polynomial)
    symmetry = PARALLEL_TROP_SYMMETRY[]
    return isnothing(symmetry) ? [parallel_trop_gfan_polynomial(polynomial)] :
        _parallel_trop_polynomial_orbit(polynomial, symmetry.actions)
end


function _parallel_trop_normalize_witness(f)
    iszero(f) && return f
    return inv(leading_coefficient(f)) * f
end


function _parallel_trop_normalized_initial_form(polynomial, weight)
    polynomial_terms = collect(terms(polynomial))
    term_weights = [dot(weight, first(exponents(term))) for term in polynomial_terms]
    maximum_weight = maximum(term_weights)
    initial_form = sum(polynomial_terms[findall(==(maximum_weight), term_weights)])
    return _parallel_trop_normalize_witness(initial_form)
end


function _parallel_trop_weighted_groebner_basis(state, weight)
    ordering = _parallel_trop_weighted_ordering(state, weight)
    return collect(groebner_basis(
        state.ideal;
        ordering=ordering,
        complete_reduction=true,
    ))
end

function _parallel_trop_weighted_ordering(state, weight)
    R = state.ring
    return weight_ordering(ones(Int, ngens(R)), weight_ordering(weight, lex(R)))
end

"Test the complete initial ideal, lifting a monomial when its torus is empty."
function _parallel_trop_exact_initial_witness(state, basis, weight)
    R = state.ring
    initial_ideal = ideal(R, _parallel_trop_normalized_initial_form.(basis, Ref(weight)))
    variables = gens(R)
    product = prod(variables)
    torus_ideal = ideal(R, product)
    isone(saturation(initial_ideal, torus_ideal)) || return nothing
    saturated, exponent = saturation_with_index(initial_ideal, torus_ideal)
    isone(saturated) || error("inconsistent monomial-containment saturation")
    exponent >= 0 || error("invalid saturation index")
    monomial = product ^ Int(exponent)

    # Removing factors while retaining exact ideal membership lowers the
    # lifted witness degree. Bound this optional reduction work independently
    # of the saturation index; an unreduced monomial is still a valid witness.
    attempts = 0
    for variable in variables
        for _ in 1:Int(exponent)
            attempts >= 128 && break
            attempts += 1
            trial = divexact(monomial, variable)
            iszero(normal_form(trial, initial_ideal)) || break
            monomial = trial
        end
        attempts >= 128 && break
    end
    ordering = _parallel_trop_weighted_ordering(state, weight)
    witness = monomial - normal_form(monomial, state.ideal; ordering=ordering)
    iszero(witness) && error("monomial lifting produced a zero witness")
    iszero(normal_form(witness, state.ideal; ordering=ordering)) ||
        error("lifted monomial witness does not belong to the defining ideal")
    empty_directions = Vector{Vector{Int}}()
    _parallel_trop_cone_leaves_weight_cell(witness, weight, empty_directions, empty_directions) == :monomial ||
        error("lifted witness does not have a monomial weight-initial form")
    return witness
end


function _parallel_trop_cone_leaves_weight_cell(
    polynomial,
    weight,
    cone_rays,
    lineality,
)
    exponent_vectors = collect(exponents(polynomial))
    term_weights = [dot(weight, exponent) for exponent in exponent_vectors]
    maximum_weight = maximum(term_weights)
    maximal_terms = findall(==(maximum_weight), term_weights)
    length(maximal_terms) == 1 && return :monomial

    reference = exponent_vectors[first(maximal_terms)]
    for term_index in Iterators.drop(maximal_terms, 1)
        difference = reference - exponent_vectors[term_index]
        any(direction -> !iszero(dot(direction, difference)), cone_rays) &&
            return :crosses_groebner_cone
        any(direction -> !iszero(dot(direction, difference)), lineality) &&
            return :crosses_groebner_cone
    end

    for term_index in eachindex(exponent_vectors)
        term_index in maximal_terms && continue
        difference = reference - exponent_vectors[term_index]
        any(direction -> dot(direction, difference) < 0, cone_rays) &&
            return :crosses_groebner_cone
        any(direction -> !iszero(dot(direction, difference)), lineality) &&
            return :crosses_groebner_cone
    end

    return :contained
end


"""Prefer interior-excluding witnesses, then lower degree and fewer terms."""
function _parallel_trop_select_witness(basis, job)
    best = nothing
    best_key = nothing
    for (index, polynomial) in enumerate(basis)
        reason = _parallel_trop_cone_leaves_weight_cell(
            polynomial,
            job.weight,
            job.rays,
            job.lineality,
        )
        reason == :contained && continue
        # A monomial initial form removes the tested interior from the
        # prevariety. A crossing witness may only subdivide a tropical cone.
        # Keep the original basis order as a deterministic final tie-break.
        key = (
            reason == :monomial ? 0 : 1,
            total_degree(polynomial),
            length(polynomial),
            index,
        )
        if isnothing(best_key) || key < best_key
            best = (polynomial=polynomial, reason=reason)
            best_key = key
        end
    end
    return best
end

"Prove one term is uniquely maximal on the entire relative interior, in MAX convention."
function _parallel_trop_relative_interior_dominance(polynomial, job)
    exponent_vectors = [BigInt.(exponent) for exponent in exponents(polynomial)]
    isempty(exponent_vectors) && return nothing
    nvariables = length(first(exponent_vectors))
    length(job.weight) == nvariables && all(length(direction) == nvariables for
        direction in Iterators.flatten((job.rays, job.lineality))) ||
        throw(ArgumentError("term dominance directions have the wrong ambient dimension"))
    pairing(direction, difference) = sum(BigInt(value) * power for
        (value, power) in zip(direction, difference); init=BigInt(0))
    weights = [pairing(job.weight, exponent) for exponent in exponent_vectors]
    maximal = findall(==(maximum(weights)), weights)
    length(maximal) == 1 || return nothing
    dominant = only(maximal)
    comparisons = NamedTuple[]
    for index in eachindex(exponent_vectors)
        index == dominant && continue
        difference = exponent_vectors[dominant] - exponent_vectors[index]
        ray_differences = [pairing(ray, difference) for ray in job.rays]
        lineality_differences = [pairing(direction, difference) for direction in job.lineality]
        all(>=(0), ray_differences) && any(>(0), ray_differences) || return nothing
        all(iszero, lineality_differences) || return nothing
        push!(comparisons, (competitor_index=index, competitor_exponent=exponent_vectors[index],
                            ray_differences=ray_differences,
                            lineality_differences=lineality_differences))
    end
    return (dominant_index=dominant, dominant_exponent=exponent_vectors[dominant],
            comparisons=comparisons)
end

"Exclude a whole relative interior only after exact ray and lineality comparisons."
function parallel_trop_whole_interior_exclusion(job)
    start_time = time()
    state = parallel_trop_state()
    length(job.weight) == ngens(state.ring) || throw(ArgumentError("weight has the wrong length"))
    basis = _parallel_trop_weighted_groebner_basis(state, job.weight)
    candidates = sort!(collect(basis); by=polynomial -> (total_degree(polynomial), length(polynomial), string(polynomial)))
    for polynomial in candidates
        witness = _parallel_trop_normalize_witness(polynomial)
        proof = _parallel_trop_relative_interior_dominance(witness, job)
        isnothing(proof) && continue
        return (excluded=true, reason=:whole_relative_interior_monomial,
                witness=parallel_trop_gfan_polynomial(witness),
                witness_orbit=_parallel_trop_witness_orbit(witness), proof=proof,
                basis_size=length(basis), seconds=time() - start_time)
    end
    lifted = _parallel_trop_exact_initial_witness(state, basis, job.weight)
    if !isnothing(lifted)
        witness = _parallel_trop_normalize_witness(lifted)
        proof = _parallel_trop_relative_interior_dominance(witness, job)
        if !isnothing(proof)
            return (excluded=true, reason=:whole_relative_interior_lifted_monomial,
                    witness=parallel_trop_gfan_polynomial(witness),
                    witness_orbit=_parallel_trop_witness_orbit(witness), proof=proof,
                    basis_size=length(basis), seconds=time() - start_time)
        end
    end
    return (excluded=false, reason=:no_whole_relative_interior_proof,
            witness=nothing, proof=nothing, basis_size=length(basis), seconds=time() - start_time)
end


"""Certify a complete candidate cone or return a polynomial witness."""
function parallel_trop_cone_certificate(job)
    start_time = time()
    state = parallel_trop_state()
    weight = job.weight
    length(weight) == ngens(state.ring) ||
        throw(ArgumentError("weight has the wrong length"))

    basis = _parallel_trop_weighted_groebner_basis(state, weight)
    candidate = _parallel_trop_select_witness(basis, job)
    if isnothing(candidate) || candidate.reason != :monomial
        lifted = _parallel_trop_exact_initial_witness(state, basis, weight)
        isnothing(lifted) || (candidate = (polynomial=lifted, reason=:monomial_saturation))
    end
    if !isnothing(candidate)
        witness = _parallel_trop_normalize_witness(candidate.polynomial)
        return (
            valid=false,
            reason=candidate.reason,
            witness=parallel_trop_gfan_polynomial(witness),
            witness_orbit=_parallel_trop_witness_orbit(witness),
            basis_size=length(basis),
            seconds=time() - start_time,
        )
    end

    return (
        valid=true,
        reason=:contained,
        witness=nothing,
        basis_size=length(basis),
        seconds=time() - start_time,
    )
end


"""Return a canonical textual signature of the weight-initial ideal."""
function parallel_trop_initial_ideal_signature(weight::Vector{Int})
    start_time = time()
    state = parallel_trop_state()
    length(weight) == ngens(state.ring) ||
        throw(ArgumentError("weight has the wrong length"))
    basis = _parallel_trop_weighted_groebner_basis(state, weight)
    initial_forms = parallel_trop_gfan_polynomial.(
        _parallel_trop_normalized_initial_form.(basis, Ref(weight)),
    )
    sort!(initial_forms)
    return (
        signature=join(initial_forms, '\n'),
        basis_size=length(basis),
        seconds=time() - start_time,
    )
end


"""Membership-only helper for diagnostics and unit checks."""
function parallel_trop_monomial_witness(weight::Vector{Int})
    empty_directions = Vector{Vector{Int}}()
    return parallel_trop_cone_certificate((
        weight=weight,
        rays=empty_directions,
        lineality=empty_directions,
    ))
end

"Probe bounded, strictly positive ray sums without certifying a whole cone."
function parallel_trop_support_probe(job)
    weights = Vector{Vector{Int}}([job.weight])
    for sample in 1:3
        coefficients = [1 + mod(5 * index^2 + (sample + 2) * index +
                               11 * job.representative, 13) for index in eachindex(job.rays)]
        weight = zeros(Int, length(job.weight))
        for (coefficient, ray) in zip(coefficients, job.rays)
            weight .+= coefficient .* ray
        end
        weight in weights || push!(weights, weight)
    end
    elapsed = 0.0
    tested = Vector{Vector{Int}}()
    last_result = nothing
    for (sample, weight) in enumerate(weights)
        push!(tested, weight)
        result = parallel_trop_monomial_witness(weight)
        elapsed += result.seconds
        last_result = result
        result.valid || return merge(result, (
            seconds=elapsed, sample_index=sample, probe_weights=tested,
        ))
    end
    return merge(last_result, (seconds=elapsed, sample_index=0, probe_weights=tested))
end
