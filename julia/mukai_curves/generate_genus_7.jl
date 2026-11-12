using Oscar

const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "data", "mukai_curves"))

R, (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) = polynomial_ring(QQ, ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p"])

G = ideal(R, [
    k*m - j*n + h*o - e*p,
    k*l - i*n + g*o - d*p,
    j*l - i*m + f*o - c*p,
    h*l - g*m + f*n - b*p,
    e*l - d*m + c*n - b*o,
    h*i - g*j + f*k - a*p,
    e*i - d*j + c*k - a*o,
    e*g - d*h + b*k - a*n,
    e*f - c*h + b*j - a*m,
    d*f - c*g + b*i - a*l
])

TV = tropical_variety(G; skip_decomposition=true)
println(TV)
PC = polyhedral_complex(TV)
println("f-vector: ", f_vector(PC))

mkpath(OUTPUT_DIR)
save(joinpath(OUTPUT_DIR, "genus_7_fan.mrdi"), PC)
