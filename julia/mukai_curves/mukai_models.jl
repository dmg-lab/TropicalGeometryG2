# The Mukai models of genus 7, 8, 9 and 10 (the last one is G2/P2 in PP^13, by the 4x4 Pfaffians of a 7x7 matrix),
# and the dimension and degree of a random linear section in PP^(genus-2) of the model of genus `genus`: for genus 10
# the 18 points of (d) in Section 10 of the paper (Krull dimension 1, degree 18).
#   julia --project=julia julia/mukai_curves/mukai_models.jl

using Oscar

# g=7: ring in 16 variables, ideal of the orthogonal Grassmannian OG(2,5) in PP^15 via 10 quadrics
R, (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) = polynomial_ring(QQ, ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p"])
G = ideal(R, [k*m-j*n+h*o-e*p, k*l-i*n+g*o-d*p, j*l-i*m+f*o-c*p,
              h*l-g*m+f*n-b*p, e*l-d*m+c*n-b*o, h*i-g*j+f*k-a*p, e*i-d*j+c*k-a*o,
              e*g-d*h+b*k-a*n, e*f-c*h+b*j-a*m, d*f-c*g+b*i-a*l])

# g=8: ring in 15 variables, Plucker ideal of the classical Grassmannian Gr(2,6) in PP^14
R, (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) = polynomial_ring(QQ, ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"])
G = grassmann_pluecker_ideal(R,2,6)

# g=9: ring in 14 variables, ideal of the Lagrangian Grassmannian LGr(3,6) in PP^13 via 21 quadrics
R, (a,b,c,d,e,f,g,h,i,j,k,l,m,n) = polynomial_ring(QQ, ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"])
G = ideal(R, [f*k-e*l+j*n, e*k+f*m-i*n, f*j-b*k-d*l, e*j-d*k+b*m-a*n,
              f*i+d*k-c*l+a*n, e*i-c*k-d*m, g*h+2*d*k-c*l-b*m+a*n, f*h-j*k+i*l, k^2+l*m+h*n,
              e*h+i*k+j*m, d*h+i*j-a*k, c*h+i^2+a*m, b*h-j^2+a*l, f^2+g*l+b*n, e*f+g*k-d*n,
              a*f-b*i-d*j, e^2-g*m-c*n, d*e-c*f-g*i, b*e+d*f+g*j, a*e+d*i-c*j, b*c+d^2+a*g])

# g=10: ring in 14 variables, ideal of the G2 orbit closure in PP^13 via 4-Pfaffians of a 7x7 matrix
R, (a,b,c,d,e,f,g,h,i,j,k,l,m,n) = polynomial_ring(QQ, ["a","b","c","d","e","f","g","h","i","j","k","l","m","n"])
M = matrix(R, [0 -f  e  g  h    i  a;
               f  0 -d  j  k    l  b;
              -e  d  0  m  n -g-k  c;
              -g -j -m  0  c   -b  d;
              -h -k -n -c  0    a  e;
              -i -l g+k  b -a   0  f;
              -a -b -c -d -e   -f  0])
G = ideal(R, pfaffians(M, 4))

# forward problem: set genus, define target ring S, and construct a random linear map L from R to S
genus = 10
S, x = polynomial_ring(QQ, ["x$i" for i in 0:genus-2])
L = hom(R, S, [sum(rand(-9:9) * xi for xi in gens(S)) for _ in 1:ngens(R)])
# compute the ideal of the linear section L(G) and display its dimension and degree
IX = L(G)
println("genus $genus: Krull dimension ", dim(IX), ", degree ", degree(IX))
