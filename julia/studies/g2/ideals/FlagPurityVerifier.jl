using SHA
using TOML

export verify_flag_tropical_purity

const FLAG_PURITY_CERTIFICATE = TropicalCensus.work_path("g2", "flag", "purity", "flag_ideal_certificate.toml")
const FLAG_GF2_PURITY_CERTIFICATE =
    TropicalCensus.work_path("g2", "flag_gf2", "purity", "flag_ideal_certificate.toml")
const FLAG_PURITY_VERIFICATION_POLICY = "flag-prime-binomial-weight-initial-v1"
const FLAG_PURITY_WEIGHT_MAX = -[6, 0, -2, 0, 5, 2, 0, -6, -8, -6, -2, -4, -6, -6,
                               -4, -8, -6, -6, -2, -10, -8, -8, -8, -6, -6, -6, -8, -6]

function _flag_ideal_digest(state)
    polynomials = sort!([string(polynomial / leading_coefficient(polynomial))
                        for polynomial in gens(state.ideal) if !iszero(polynomial)])
    payload = join(string.(symbols(state.ring)), ',') * "\n" * join(polynomials, '\n') * "\n"
    return bytes2hex(sha256(payload))
end

function _write_flag_purity_certificate(path, certificate)
    mkpath(dirname(path))
    temporary = path * ".tmp.$(getpid())"
    open(temporary, "w") do io
        TOML.print(io, certificate; sorted=true)
    end
    mv(temporary, path; force=true)
    return path
end

function _flag_purity_result(certificate, path; reused=false)
    verified = get(certificate, "verified", false)
    return (
        verified=verified,
        exact=get(certificate, "exact", false),
        pure=verified && get(certificate, "geometrically_prime", false),
        dimension=get(certificate, "dimension", -1),
        prime=verified && get(certificate, "prime", false),
        geometrically_prime=verified && get(certificate, "geometrically_prime", false),
        torus_saturated=verified && get(certificate, "torus_saturated", false),
        method=get(certificate, "method", ""),
        details=get(certificate, "details", ""),
        defining_ideal_sha256=get(certificate, "defining_ideal_sha256", ""),
        oscar_version=get(certificate, "oscar_version", ""),
        verification_policy=get(certificate, "verification_policy", ""),
        reused_certificate=reused,
        certificate_path=path,
    )
end

"""
Prove flag purity from a freshly computed geometrically prime weight-initial
ideal. A monomial-free binomial ideal saturated by the variable product is a
lattice ideal; Smith index one makes that lattice saturated and its quotient
a domain over the algebraic closure. A positive filtration transfers this
primality to the homogeneous defining ideal. Coordinate nonmembership then
proves saturation of the defining prime ideal by the full variable product.
Only a completed certificate bound to the actual defining polynomials and
Oscar version can be reused.

The coefficients are in QQ or in GF(2). A lattice ideal of a saturated lattice
is prime over every field, so the same argument holds in characteristic two;
its certificate is kept apart and records the field.
"""
function verify_flag_tropical_purity(state; certificate_path=nothing, reuse_only::Bool=false)
    expected_names = vcat(["x$i" for i in 1:7],
                          ["p$i$j" for i in 1:6 for j in (i + 1):7])
    string.(symbols(state.ring)) == expected_names ||
        error("flag purity verification requires the full 28-coordinate ring")
    field = base_ring(state.ring)
    field_name = field == QQ ? "QQ" :
        field isa FinField && order(field) == 2 ? "GF(2)" :
        error("flag purity verification requires QQ or GF(2)")
    certificate_path = something(certificate_path,
        field == QQ ? FLAG_PURITY_CERTIFICATE : FLAG_GF2_PURITY_CERTIFICATE)
    polynomials = collect(gens(state.ideal))
    length(polynomials) == 78 || error("flag purity verification requires 78 defining equations")
    homogeneous = all(all(sum(exponent) == total_degree(polynomial) for exponent in
                          exponents(polynomial)) for polynomial in polynomials)
    homogeneous || error("flag initial-prime transfer requires homogeneous equations")
    digest = _flag_ideal_digest(state)
    version = string(Base.pkgversion(Oscar))
    if isfile(certificate_path)
        certificate = TOML.parsefile(certificate_path)
        current = get(certificate, "status", "") == "complete" &&
                  get(certificate, "coefficient_field", "QQ") == field_name &&
                  get(certificate, "defining_ideal_sha256", "") == digest &&
                  get(certificate, "oscar_version", "") == version &&
                  get(certificate, "verification_policy", "") == FLAG_PURITY_VERIFICATION_POLICY &&
                  get(certificate, "weight_max", Int[]) == FLAG_PURITY_WEIGHT_MAX
        recorded_path = get(certificate, "initial_ideal_path", "")
        initial_path = isabspath(recorded_path) ? recorded_path : joinpath(dirname(certificate_path), recorded_path)
        artifact_current = isfile(initial_path) &&
                           get(certificate, "initial_ideal_sha256", "") == bytes2hex(sha256(read(initial_path)))
        consistent = artifact_current && all(get(certificate, flag, false) === true for flag in (
            "verified", "exact", "homogeneous", "initial_binomial", "initial_monomial_free",
            "initial_torus_saturated", "initial_geometrically_prime", "initial_coordinates_nonzero",
            "prime", "geometrically_prime", "torus_saturated")) &&
            get(certificate, "initial_lattice_index", "") == "1" &&
            get(certificate, "dimension", -1) == get(certificate, "initial_dimension", -2) == 8 &&
            get(certificate, "variables", -1) == 28 && get(certificate, "generators", -1) == 78
        current && consistent && return _flag_purity_result(certificate, certificate_path; reused=true)
    end
    if reuse_only
        return _flag_purity_result(Dict{String, Any}(
            "verified" => false, "exact" => false,
            "defining_ideal_sha256" => digest, "oscar_version" => version,
            "verification_policy" => FLAG_PURITY_VERIFICATION_POLICY,
            "method" => "Completed exact prime-initial certificate lookup",
            "details" => "No matching completed certificate exists; no computation was started.",
        ), certificate_path)
    end
    certificate = Dict{String, Any}(
        "status" => "initial_running", "verified" => false, "exact" => false,
        "coefficient_field" => field_name,
        "homogeneous" => homogeneous, "variables" => 28, "generators" => 78,
        "defining_ideal_sha256" => digest, "oscar_version" => version,
        "verification_policy" => FLAG_PURITY_VERIFICATION_POLICY,
        "weight_max" => FLAG_PURITY_WEIGHT_MAX,
        "positive_filtration_weights" => FLAG_PURITY_WEIGHT_MAX .+ (1 - minimum(FLAG_PURITY_WEIGHT_MAX)),
        "method" => "Exact prime binomial weight-initial ideal: product saturation and Smith lattice index one; positive filtration transfers primality",
        "details" => "The initial ideal is a geometrically prime lattice ideal. Homogeneity permits an all-ones shift to positive filtration weights, preserving the weight-initial ideal. A domain associated graded ring implies a domain defining quotient. No coordinate belongs to the initial ideal, hence none belongs to the defining prime ideal, proving product saturation and pure dimension eight.",
    )
    _write_flag_purity_certificate(certificate_path, certificate)
    println("Flag purity: recomputing an exact 28-coordinate candidate prime initial ideal")
    flush(stdout)
    started = time()
    initial_ideal = initial(state.ideal, tropical_semiring_map(field, min), QQ.(-FLAG_PURITY_WEIGHT_MAX))
    reduced = collect(gens(groebner_basis(initial_ideal; complete_reduction=true)))
    certificate["initial_seconds"] = time() - started
    certificate["initial_generators"] = length(reduced)
    certificate["dimension"] = Int(dim(state.ideal))
    certificate["initial_dimension"] = Int(dim(initial_ideal))
    certificate["initial_binomial"] = all(length(polynomial) == 2 for polynomial in reduced)
    certificate["initial_lattice_index"] = certificate["initial_binomial"] ?
        string(TropicalCensus.binomial_lattice_index(reduced)) : "not-binomial"
    certificate["status"] = "saturation_running"
    _write_flag_purity_certificate(certificate_path, certificate)
    started = time()
    saturated = saturation(initial_ideal, ideal(state.ring, prod(gens(state.ring))))
    certificate["saturation_seconds"] = time() - started
    certificate["initial_monomial_free"] = !isone(saturated)
    certificate["initial_torus_saturated"] = saturated == initial_ideal
    certificate["initial_coordinates_nonzero"] = all(!(variable in initial_ideal) for variable in gens(state.ring))
    initial_prime = certificate["initial_binomial"] && certificate["initial_lattice_index"] == "1" &&
                    certificate["initial_monomial_free"] && certificate["initial_torus_saturated"]
    certificate["initial_geometrically_prime"] = initial_prime
    certificate["prime"] = initial_prime
    certificate["geometrically_prime"] = initial_prime
    certificate["torus_saturated"] = initial_prime && certificate["initial_coordinates_nonzero"]
    certificate["exact"] = true
    certificate["verified"] = initial_prime && certificate["torus_saturated"] &&
                              certificate["dimension"] == certificate["initial_dimension"] == 8
    certificate["status"] = "complete"
    # Workers may compute a fresh proof concurrently. Publish immutable,
    # content-addressed evidence before atomically publishing its TOML, so
    # another worker cannot replace the evidence behind a completed digest.
    initial_temporary = joinpath(dirname(certificate_path),
        "verified_prime_initial_ideal.tmp.$(getpid()).mrdi")
    save(initial_temporary, initial_ideal)
    initial_digest = bytes2hex(sha256(read(initial_temporary)))
    initial_path = joinpath(dirname(certificate_path),
        "verified_prime_initial_ideal_$initial_digest.mrdi")
    mv(initial_temporary, initial_path; force=true)
    certificate["initial_ideal_path"] = basename(initial_path)
    certificate["initial_ideal_sha256"] = initial_digest
    _write_flag_purity_certificate(certificate_path, certificate)
    println("Flag purity: verified=$(certificate["verified"]), exact initial+ saturation time=$(round(certificate["initial_seconds"] + certificate["saturation_seconds"]; digits=2))s")
    flush(stdout)
    return _flag_purity_result(certificate, certificate_path)
end
