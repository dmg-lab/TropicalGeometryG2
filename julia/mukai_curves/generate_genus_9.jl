using Oscar

const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "data", "mukai_curves"))

R, (a,b,c,d,e,f,g,h,i,j,k,l,m,n) = polynomial_ring(QQ, ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"])

G = ideal(R, [
    f*k - e*l + j*n,
    e*k + f*m - i*n,
    f*j - b*k - d*l,
    e*j - d*k + b*m - a*n,
    f*i + d*k - c*l + a*n,
    e*i - c*k - d*m,
    g*h + 2*d*k - c*l - b*m + a*n,
    f*h - j*k + i*l,
    k^2 + l*m + h*n,
    e*h + i*k + j*m,
    d*h + i*j - a*k,
    c*h + i^2 + a*m,
    b*h - j^2 + a*l,
    f^2 + g*l + b*n,
    e*f + g*k - d*n,
    a*f - b*i - d*j,
    e^2 - g*m - c*n,
    d*e - c*f - g*i,
    b*e + d*f + g*j,
    a*e + d*i - c*j,
    b*c + d^2 + a*g
])


TV = tropical_variety(G;skip_saturation=false, skip_decomposition=true)
println(TV)
PC = polyhedral_complex(TV)
println("f-vector: ", f_vector(PC))


mkpath(OUTPUT_DIR)
save(joinpath(OUTPUT_DIR, "genus_9_fan.mrdi"), PC)
