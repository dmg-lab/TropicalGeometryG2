###############################################################################
#
#   Generate the Chevalley polytopes of arXiv:2411.10276 and write them as
#   polymake files.
#
#   For each homogeneous space in Appendix A this walks ALL reduced expressions
#   for w^P up to commutation -- which is the full set of distinct Chevalley
#   polytopes for that space, since the polytope is invariant under commutation
#   moves -- and writes one .poly file per commutation class.
#
#   Each file is a native polymake object written with `Polymake.save` on the
#   `pm_polytope` subobject, and carries its provenance as polymake
#   attachments.  From polymake:
#
#       polytope > $p = load("data/chevalley-G2-varpi_0-1-word_2-1-2-1-2.poly");
#       polytope > print $p->get_attachment("CHEVALLEY_WORD");
#       2 1 2 1 2
#       polytope > print $p->N_VERTICES;
#       6
#
#   Run with:   julia generate_chevalley_polytopes.jl [outdir]
#
###############################################################################

using TropicalCensus
using Oscar

const OUTDIR = length(ARGS) >= 1 ? ARGS[1] : data_path("chevalley")

# The five homogeneous spaces of Appendix A.
const SPACES = [
  (:A, 3, [1, 0, 1], "Fl(1,3;4) = SL_4/P_{1,3}, adjoint     (A.1)"),
  (:B, 2, [0, 3], "OG(2,5) = Spin_5/P_2, varpi = 3w_2   (A.2)"),
  (:F, 4, [0, 0, 0, 1], "F_4/P_4, coadjoint                   (A.3)"),
  (:G, 2, [0, 1], "G_2/P_2, adjoint                     (A.4)"),
  (:C, 4, [0, 1, 0, 0], "SG(2,8) = Sp_8/P_2, coadjoint        (A.5)"),
]

mkpath(OUTDIR)
println("writing polymake files to ", OUTDIR, "\n")

total = 0
for (type, rk, varpi, label) in SPACES
  R = root_system(type, rk)
  W = weyl_group(R)
  wP = minimal_coset_representative(W, parabolic_indices(varpi))
  words = all_paper_words(wP; up_to_commutation=true)
  d, deg = degree_of_embedding(R, varpi)

  println("=== ", label)
  println("    dim X = ", d, ", deg X = ", deg,
    ", commutation classes = ", length(words))

  for w in words
    D = chevalley_data(type, rk, varpi; word=w)
    P = chevalley_polytope(D)
    vol = normalised_volume_of(P)
    path = save_chevalley_polytope(D, OUTDIR)
    global total += 1
    println("    word ", rpad(join(w, " "), 32),
      " vol ", rpad(vol, 5),
      vol == deg ? " = deg  " : " < deg  ",
      "verts ", rpad(n_vertices(P), 4),
      basename(path))
  end
  println()
end

println("wrote ", total, " polymake files to ", OUTDIR)
