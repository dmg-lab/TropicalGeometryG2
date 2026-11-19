# Verify the certificates of compute.jl WITHOUT recomputing any Groebner basis of J (as
# ../verify.jl for the Grassmannian).  For each orbit representative C (cone_NNNN.mrdi), the
# pieces P cover C:  cover = every ray of every P lies in C, the pieces have pairwise
# intersections of lower dimension, and their volumes sum to that of C (slices
# {c <= 1} with c positive on C modulo lineality).  For each piece, with w = sum of its rays:
#   in_I      G in J: every g lies in <G0>_ZZ (division by G0, identity checked), where G0
#             (certificates/base.mrdi) lies in J by its stored cofactors (checked once)
#   J_in_G    every generator of J lies in <G>_ZZ (division by G, identity checked)
#   spairs    every S-polynomial over ZZ of G (non-coprime leading monomials) has a standard
#             representation with all terms below the lcm (identity and order checked)
#   lcs       leading coefficients of G (powers of 2: G reduces to a Groebner basis mod odd p)
#   lin_hom   every g is homogeneous for the lineality grading of C
#   rays      the ray conditions of ray_violation (common_flag.jl): in_v(g) = in_w(g) on relint P
#   torus     the +-1 point is a common zero of all in_w(g), or the slice certificate:
#             slice_gens recomputed, slice_gens in <slice_G>_ZZ, S-polynomials of slice_G with
#             standard representations, no constant, leading coefficients powers of 2
# Also checks that the 12 coordinate permutations of the signed action (basis.jl) carry the
# representatives onto all 2714 maximal cones.
# usage: julia --project=julia --heap-size-hint=2500M verify.jl [--resume]   -> results/verify.csv
include("common_flag.jl")
mkpath(FRES)
R, gens_Z = flag_ideal_ZZ()
F = flag_fan()
cones = collect(maximal_cones(F))
basis = load(joinpath(FCERT, "basis.mrdi"))
perms = [coordinate_permutation(Vector{Int}(s)) for s in basis["weyl_sigma"]]
orbit, img, reps = cone_orbits(F, perms)
sizes = [count(==(o), orbit) for o in 1:length(reps)]
covered = Set(t[k] for k in reps for t in img)
println("orbits: ", length(reps), ", the representatives and their images cover all ", length(cones), " cones: ",
        covered == Set(1:length(cones))); flush(stdout)
canon = canonical_ray_map(F)
base = load(joinpath(FCERT, "base.mrdi"))
G0 = [evaluate(g, gens(R)) for g in base["G0"]]; o0 = weight_order(R, Vector{Int}(base["w0"]))
base_ok = check_cofactors(G0, [[evaluate(c, gens(R)) for c in row] for row in base["cofactors"]], gens_Z)
println("base: |G0| = ", length(G0), ", G0 in J over ZZ by its cofactors: ", base_ok); flush(stdout)
@assert base_ok
pow2(n) = n > 0 && count_ones(n) == 1
hrep(C) = (Matrix{QQFieldElem}(linear_inequality_matrix(facets(C))), Matrix{QQFieldElem}(linear_equation_matrix(linear_span(C))))
inside(H, v) = all(<=(0), H[1] * v) && (isempty(H[2]) || all(iszero, H[2] * v))
cone_of(rs, lin) = positive_hull(matrix(QQ, reduce(vcat, [permutedims(QQ.(r)) for r in rs])),
                                 matrix(QQ, reduce(vcat, [permutedims(v) for v in lin])))

"Do the pieces (ray lists) cover C?  Containment, interiors pairwise disjoint, volumes add up."
function covers(C, pieces_rays, lin)
  H = hrep(C)
  all(rs -> all(r -> inside(H, QQ.(r)), rs), pieces_rays) || return false
  Ps = [cone_of(rs, lin) for rs in pieces_rays]
  all(i -> all(j -> dim(intersect(Ps[i], Ps[j])) < dim(C), (i + 1):length(Ps)), eachindex(Ps)) || return false
  vol = slice_volume_function(C, canon)
  sum(vol(rs) for rs in pieces_rays) == vol(primitive_rays(C))
end

path = joinpath(FRES, "verify.csv")
resume = "--resume" in ARGS && isfile(path)
done = Set{Int}()
if resume                          # keep the rows of completely verified cones (piece == pieces)
  lines = filter(l -> !startswith(l, "#"), readlines(path))
  for l in lines[2:end]; f = split(l, ','); f[3] == f[5] && push!(done, parse(Int, f[1])); end
  open(path, "w") do io
    foreach(l -> println(io, l), [lines[1]; filter(l -> parse(Int, split(l, ',')[1]) in done, lines[2:end])])
  end
end
out = open(path, resume ? "a" : "w")
resume || println(out, "cone,orbit_size,pieces,cover,piece,gb_size,in_I,J_in_G,spairs_checked,spairs_failed,lcs,lcs_powers_of_2,lin_hom,rays,torus,torus_kind,slice_lcs,seconds")
flush(out)
for k in reps
  k in done && continue
  t0 = time()
  C = cones[k]; lin = lineality_generators(C)
  pieces = load(flag_cone_file(k))["pieces"]
  cover = covers(C, [[Vector{Int}(r) for r in p["rays"]] for p in pieces], lin)
  for (n, d) in enumerate(pieces)
    rs = [Vector{Int}(r) for r in d["rays"]]
    w = Vector{Int}(d["w"]); @assert w == sum(rs)
    G = [evaluate(g, gens(R)) for g in d["G"]]
    o = weight_order(R, w)
    in_I = check_in_span(G, G0, o0)
    J_in_G = check_in_span(gens_Z, G, o)
    nchk, nfail = check_standard_representations(G, o)
    lcs = sort(unique([abs(Int(leading_coefficient(g; ordering = o))) for g in G]))
    lin_hom = all(g -> all(v -> length(unique([sum(QQ(e[i]) * v[i] for i in 1:NV) for e in exponents(g)])) <= 1, lin), G)
    rays_ok = ray_violation(G, w, rs) === nothing
    Jw = [initial_form(g, w) for g in G]
    pt = Vector{Int}(d["pm1_point"])
    slice_lcs = ""
    if !isempty(pt)
      torus = all(f -> iszero(evaluate(f, ZZ.(pt))), Jw); kind = "pm1"
    else
      kind = "slice"
      T, z = polynomial_ring(ZZ, [NAMES; "y"])
      S = Vector{Int}(d["slice_S"]); rest = [j for j in 1:NV if !(j in S)]
      sg = [[evaluate(f, z[1:NV]) for f in Jw]; [z[j] - 1 for j in S]; z[NV + 1] * prod(z[rest]) - 1]
      stored_sg = [evaluate(f, z) for f in d["slice_gens"]]
      H = [evaluate(h, z) for h in d["slice_G"]]
      oT = degrevlex(T)
      sn, sf = check_standard_representations(H, oT)
      hl = sort(unique([abs(Int(leading_coefficient(h; ordering = oT))) for h in H]))
      slice_lcs = join(hl, ";")
      torus = sg == stored_sg && check_in_span(sg, H, oT) && sf == 0 && !any(is_constant, H) && all(pow2, hl)
    end
    println(out, k, ",", sizes[orbit[k]], ",", length(pieces), ",", cover, ",", n, ",", length(G), ",", in_I, ",", J_in_G, ",",
            nchk, ",", nfail, ",\"", join(lcs, ";"), "\",", all(pow2, lcs), ",", lin_hom, ",", rays_ok, ",", torus, ",", kind,
            ",\"", slice_lcs, "\",", round(time() - t0; digits = 1))
    flush(out)
    GC.gc()
  end
end
println(out, "# done"); close(out)
