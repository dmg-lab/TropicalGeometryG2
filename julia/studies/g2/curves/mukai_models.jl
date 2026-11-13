# The Mukai model of genus 10, G2/P2 in PP^13 by the 4x4 Pfaffians of a 7x7 matrix, and the dimension and degree of a
# random linear section in PP^(genus-2): the 18 points of (d) in Section 10 of the paper (Krull dimension 1, degree 18).
#   julia --project=. studies/g2/curves/mukai_models.jl

using Oscar


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
