using Oscar
using ProgressMeter

# Primality of the torus saturations P = in_w(I) : (x_1 ⋯ x_n)^∞, over QQ and over QQbar, for every face of the raw
# fan of trop(G) (G the G2-Grassmannian), and for the 11 initial forms of the quadric q.
#
# Over QQ:    `is_prime(P)`.
# Over QQbar: `absolute_primary_decomposition(P)` returns (Q, P_QQ, P_abs, d), one entry per QQ-primary component, where
#             d is the number of conjugate geometric components. P is prime over QQbar iff there is exactly one entry,
#             with Q == P_QQ (no embedded or non-reduced part) and d == 1.
#
# From the `julia` folder:  julia --project=. studies/g2/weights/saturation_primality.jl
# writes saturation_primality.csv (one row per face) and prints a summary.

include(joinpath(@__DIR__, "weights.jl"))

const PRIMALITY_PATH = joinpath(@__DIR__, "saturation_primality.csv")

"All faces of the raw fan with their dimension and an interior weight (the sum of the primitive rays)."
function raw_fan_faces_with_weights(path=RAW_FAN_PATH)
    F = load(path)
    faces = [(dim=lineality_dim(F), weight=zeros(Int, 21))]
    for d in (lineality_dim(F) + 1):dim(F), C in cones(F, d)
        push!(faces, (dim=d, weight=sum(integer_ray(Vector(r)) for r in rays_modulo_lineality(C)[1])))
    end
    return faces
end

"Primality of the ideal `P` over QQ and over QQbar."
function primality(P)
    prime_QQ = is_prime(P)
    components = absolute_primary_decomposition(P)
    prime_QQbar = length(components) == 1 && components[1][1] == components[1][2] && components[1][4] == 1
    return (prime_QQ=prime_QQ, prime_QQbar=prime_QQbar, components=length(components),
            conjugates=[Int(c[4]) for c in components])
end

function check_grassmannian(; path=PRIMALITY_PATH)
    I = load(IDEAL_PATH)
    faces = raw_fan_faces_with_weights()
    rows = String[]
    @showprogress desc = "faces: " for (k, face) in enumerate(faces)
        L = laurent_initial_ideal(I, face.weight)
        p = primality(L.saturation)
        push!(rows, join((k, face.dim, L.initial == L.saturation, p.prime_QQ, p.prime_QQbar, p.components,
                          join(p.conjugates, ";")), ","))
    end
    open(path, "w") do io
        println(io, "face,dim,already_saturated,prime_QQ,prime_QQbar,QQ_components,conjugates")
        foreach(r -> println(io, r), rows)
    end
    return path
end

"The initial forms of q at the 11 cones of the tropical quadric, and their primality over QQ and QQbar."
function check_quadric()
    R, x = polynomial_ring(QQ, ["x$i" for i in 1:7])
    terms = [x[1] * x[4], x[2] * x[5], x[3] * x[6], -x[7]^2]
    cones = [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4], [1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4], [1, 2, 3, 4]]
    return [(S, primality(ideal(R, [sum(terms[S])]))) for S in cones]
end

function summarize_primality(path=PRIMALITY_PATH)
    rows = [split(l, ",") for l in Iterators.drop(eachline(path), 1) if !isempty(l)]
    for d in sort(unique(r[2] for r in rows))
        rs = filter(r -> r[2] == d, rows)
        println("dim $d: $(length(rs)) faces; prime over QQ: $(count(r -> r[4] == "true", rs)); ",
                "prime over QQbar: $(count(r -> r[5] == "true", rs))")
    end
    bad = [(r[1], r[2], r[4], r[5], r[6], r[7]) for r in rows if r[4] != "true" || r[5] != "true"]
    println("faces not prime over QQ or QQbar (face, dim, QQ, QQbar, components, conjugates): ", bad)
end

if abspath(PROGRAM_FILE) == @__FILE__
    for (S, p) in check_quadric()
        println("quadric cone $S: prime over QQ $(p.prime_QQ), over QQbar $(p.prime_QQbar)")
    end
    summarize_primality(check_grassmannian())
end
