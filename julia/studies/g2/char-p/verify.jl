# Verify the stored certificates of the second inclusion in the proof of Theorem thm:charp
# (main.tex) WITHOUT recomputing any Groebner basis.  For each cone file
# certificates/cones/cone_NNN.mrdi (see compute.jl) and the cone C of the raw fan:
#   in_I      G is contained in I: the stored cofactors give G[k] = sum_j c_kj gens_j over ZZ
#   I_in_G    every generator of I lies in <G>_ZZ (division by G, identity checked)
#   spairs    every S-polynomial of G (non-coprime leading monomials) has a standard
#             representation over ZZ with all terms below the lcm (identity + order checked)
#   lcs       leading coefficients of G (must be powers of 2 for every odd p to be good)
#   lin_hom   every g in G is homogeneous for the lineality grading of C
#   H, I      terms of in_w(g) tie at every ray of C, and are <= the other terms there (MIN)
#   torus     either the +-1 point is a common zero of all in_w(g), or the slice certificate:
#             slice_gens recomputed from in_w(G) and S, slice_gens in <slice_G>_ZZ,
#             S-polynomials of slice_G have standard representations, no constant in
#             slice_G, leading coefficients of slice_G (powers of 2 needed)
# These are exactly the hypotheses used with Lemma lem:groebnerZ in the proof.
# usage: julia verify.jl CHUNK NCHUNKS   -> results/verify-CHUNK.csv
include("common.jl")
chunk, nchunks = parse(Int, ARGS[1]), parse(Int, ARGS[2])
mkpath(joinpath(@__DIR__, "results"))
R, gens_Z = ideal_ZZ()
F = raw_fan()
cones = collect(maximal_cones(F))
pow2(n) = n > 0 && count_ones(n) == 1
out = open(joinpath(@__DIR__, "results", "verify-$chunk.csv"), "w")
println(out, "cone,gb_size,in_I,I_in_G,spairs_checked,spairs_failed,lcs,lcs_powers_of_2,lin_hom,H,I,torus,torus_kind,slice_lcs")
flush(out)
for k in chunk:nchunks:length(cones)
  d = load(cone_file(k))
  C = cones[k]; w = Vector{Int}(d["w"]); @assert w == interior_point(C)
  G = [evaluate(g, gens(R)) for g in d["G"]]
  o = weight_order(R, w)
  cof = [[evaluate(c, gens(R)) for c in row] for row in d["cofactors"]]
  in_I = check_cofactors(G, cof, gens_Z)
  I_in_G = check_in_span(gens_Z, G, o)
  nchk, nfail = check_standard_representations(G, o)
  lcs = sort(unique([abs(Int(leading_coefficient(g; ordering = o))) for g in G]))
  # homogeneity and ray conditions
  rs = primitive_rays(C); lin = lineality_generators(C)
  lin_hom = all(g -> all(v -> length(unique([sum(QQ(e[i]) * v[i] for i in 1:21) for e in exponents(g)])) <= 1, lin), G)
  Hok = true; Iok = true
  for g in G
    es = collect(exponents(g)); ws = [sum(e .* w) for e in es]; m = minimum(ws)
    ini = findall(==(m), ws); non = findall(>(m), ws)
    for r in rs
      a = sum(es[ini[1]] .* r)
      Hok &= all(i -> sum(es[i] .* r) == a, ini)
      Iok &= all(i -> sum(es[i] .* r) >= a, non)
    end
  end
  J = [initial_form(g, w) for g in G]
  pt = Vector{Int}(d["pm1_point"])
  slice_lcs = ""
  if !isempty(pt)
    torus = all(f -> iszero(evaluate(f, ZZ.(pt))), J); kind = "pm1"
  else
    kind = "slice"
    T, z = polynomial_ring(ZZ, [["p$i$j" for (i, j) in PAIRS]; "y"])
    S = Vector{Int}(d["slice_S"]); rest = [j for j in 1:21 if !(j in S)]
    sg = [[evaluate(f, z[1:21]) for f in J]; [z[j] - 1 for j in S]; z[22] * prod(z[rest]) - 1]
    stored_sg = [evaluate(f, z) for f in d["slice_gens"]]
    H = [evaluate(h, z) for h in d["slice_G"]]
    oT = degrevlex(T)
    sn, sf = check_standard_representations(H, oT)
    hl = sort(unique([abs(Int(leading_coefficient(h; ordering = oT))) for h in H]))
    slice_lcs = join(hl, ";")
    torus = sg == stored_sg && check_in_span(sg, H, oT) && sf == 0 &&
            !any(is_constant, H) && all(pow2, hl)
  end
  println(out, k, ",", length(G), ",", in_I, ",", I_in_G, ",", nchk, ",", nfail, ",\"", join(lcs, ";"), "\",",
          all(pow2, lcs), ",", lin_hom, ",", Hok, ",", Iok, ",", torus, ",", kind, ",\"", slice_lcs, "\"")
  flush(out)
end
println(out, "# done"); close(out)
