# Dihedral symmetries of the G2/P2 fan

## Result

The coarsest fan has a faithful dihedral taxon symmetry

```text
D6 = <r,s | r^6 = s^2 = 1, srs = r^(-1)>,    |D6| = 12,
```

where `r` rotates the taxon hexagon `1-5-3-4-2-6` and `s` reflects it
through the opposite vertices 1 and 4.  Both fix taxon 7.  The induced
permutations of the 21 Pluecker coordinates preserve the lineality space, all
51 quotient rays, and all 248 maximal cones.  Thus they are exact linear fan
automorphisms of R^21 and give the desired abstract dihedral action on taxa.

There is an important qualification: the six elements called *reflections* in
the abstract dihedral group are not Euclidean hyperplane reflections in this
natural 21-dimensional representation.  Their fixed spaces have codimension 8
or 9 (and not 1).  Their restrictions to the span of the fan, modulo
lineality, still have fixed-space codimension greater than 1, so changing only
their action on an ambient complement cannot repair this.  In fact, no
nonidentity element of this D6 image fixes a hyperplane.  The group is
therefore an orthogonal permutation representation of an abstract reflection
group, but it is not itself a real reflection group inside GL(21) in the
strict Coxeter sense.

Allowing signs in a monomial lift does not change this obstruction.  Each
2-cycle of coordinates in an involution contributes one `-1` eigenvalue,
independently of the compatible signs.  The two taxon-reflection types induce
8 and 9 such coordinate 2-cycles.

The fan rays span 12 of the 18
dimensions modulo the 3-dimensional lineality space.

## What remains outside lineality

The representation-theoretic decomposition explains precisely where the
literal G2 reflection representation lives.  Here `E1` is the faithful
2-dimensional representation in which `r` is rotation through 60 degrees and
the abstract reflections fix lines; `E2` is the 2-dimensional representation
that factors through a rotation of order 3.  The four `A/B` representations
are 1-dimensional.

| representation | 1 | r,+r^5 | r^2,+r^4 | r^3 | vertex refl. | edge refl. | decomposition |
|---|---:|---:|---:|---:|---:|---:|---|
| R^21 | 21 | 0 | 0 | 3 | 5 | 3 | `4A1 + 2B1 + B2 + 3E1 + 4E2` |
| lineality L | 3 | 2 | 0 | -1 | 1 | 1 | `A1 + E1` |
| fan span modulo L | 12 | -2 | 0 | 4 | 2 | 2 | `2A1 + B1 + B2 + E1 + 3E2` |
| R^21/L | 18 | -2 | 0 | 4 | 4 | 2 | `3A1 + 2B1 + B2 + 2E1 + 4E2` |

Consequently the lineality space is exactly `A1 + E1`: a fixed scaling
direction plus the expected rank-2 G2 reflection plane.  Outside lineality the
fan itself spans

```text
2A1 + B1 + B2 + E1 + 3E2,
```

and the full ambient quotient is

```text
3A1 + 2B1 + B2 + 2E1 + 4E2.
```

The entire complement is therefore not a reflection representation: on the
fan span each abstract reflection has a 5-dimensional `-1` eigenspace.
Nevertheless, its lone `E1` summand is canonical.  It is obtained by applying
the `E1` character projector and taking the orthogonal complement of
lineality.

This plane gives an exact second G2 picture.  In the rational basis selected
by row reduction, its intrinsic lattice is

```text
E intersect Z^21 = {(a,b) in Z^2 : a = b mod 3},
```

with basis `(1,1), (2,-1)`.  The primitive lattice vectors on the projected
directions are

```text
short: [(-2, 1), (-1, -1), (-1, 2), (1, -2), (1, 1), (2, -1)]
long:  [(-3, 0), (-3, 3), (0, -3), (0, 3), (3, -3), (3, 0)]
```

The inherited Gram matrix is

```text
[40/3  20/3]
[20/3  40/3],
```

so the squared long/short length ratio is exactly 3.  Taking the short root
`(-1,-1)` and long root `(3,0)` gives the G2 Cartan matrix

```text
[ 2 -3]
[-1  2].
```

Of the 51 fan rays, 36 have nonzero image and
15 vanish.  The nonzero images occupy 12 oriented
root directions, split into the short and long hexagons above, and every root
direction is hit by exactly three fan rays.  The ray-by-ray exact coordinates
are recorded in `projected_rays.csv`.

## Hidden affine G2 structure

The next incidence test gives a clean answer: there is an exact finite window
of the twisted affine root system `D4^(3)`, whose affine Weyl group is of G2
type.  It is a root pattern in a quotient of the fan's canonical isotypic
projection, not an affine Weyl symmetry of the finite fan itself.

The relevant four-dimensional projection is

```text
U = 2A1 + E1
```

where `E1` is the finite G2 root plane and `2A1` supplies a homogenizing
coordinate and an affine level.  After normalizing the `E1` part, the 36
nonzero rays give 24 distinct lifts: two lifts over each of the 12 finite G2
roots.

Cone incidence distinguishes three different roles among the fixed
projective directions:

| fixed direction | ray indices | D6 orbit sizes | degrees | maximal cells | zero-only | short/long directions | always two real roots? |
|---|---|---|---|---:|---:|---:|---|
| `(1, -6)` | `[8, 9, 10, 30, 33, 34, 37]` | `[1, 3, 3]` | `[9, 21, 21, 21, 21, 21, 21]` | 96 | 12 | 6/6 | false |
| `(1, 8)` | `[18, 35, 36]` | `[3]` | `[9, 9, 9]` | 27 | 0 | 6/6 | true |
| `(3, -4)` | `[14, 15, 21, 22, 29]` | `[2, 3]` | `[15, 15, 15, 15, 15]` | 57 | 6 | 6/0 | false |

Set

```text
c = (1,-6),      d = (1,8),      delta = (5/21)[d] in A1^2 / Qc.
```

The line `Qc` is canonically the homogenizing line: it contains the unique fan
ray fixed by all of D6, ray 37, and the nine maximal cones in that
ray's star contain only zero-finite rays.  By contrast, the three rays on
`Qd`, namely `[18, 35, 36]`, have 27 incident maximal cones and see all
12 finite roots.  Every such cone contains exactly two real-root rays.  Their
real-root links are exactly:

| delta ray | K3,3 bipartition of its real-root link |
|---:|---|
| 18 | `[6, 7, 17]` versus `[13, 40, 41]` |
| 35 | `[6, 12, 44]` versus `[11, 20, 40]` |
| 36 | `[7, 11, 45]` versus `[12, 19, 41]` |

Each link is `K3,3`.  Each side contains two long roots and one short root,
and the exact relation on each side is

```text
lambda_1 + lambda_2 = -3 alpha_short.
```

The three nonsimplicial cones in these stars are

```text
[13, 14, 17, 18]
[19, 21, 36, 45]
[20, 22, 35, 44]
```

and exact arithmetic gives, modulo lineality, the circuit relation

```text
[r_+] + [r_-] = 6([z_+] + [z_-]),
```

where `r_+` and `r_-` are the antipodal real-root rays and `z_+,z_-` are the
two zero-finite rays in the cone.  This is why the `Qd` incidence is an affine
interface, whereas ray 37 is a vertical apex.

The third direction is not independent after removing the homogenizing line:

```text
(3,-4) = d + 2c.
```

Thus `(3,-4)` and `(1,8)` define the same one-dimensional quotient direction.
In the chosen A1 coordinates, the intrinsic fixed-plane lattice is

```text
Lambda = {(a,b) in Z^2 : b is even}.
```

The primitive quotient coordinate on this lattice is

```text
q(a,b) = 3a + b/2,        ker(q) = Zc.
```

The normalized level functional

```text
n(a,b) = (3/5)q(a,b) = (3/10)(6a+b)
```

kills `c` and sends `delta` to 1.  The four normalized lift positions are

| finite roots | A1 coordinate | level | rays per root |
|---|---|---:|---:|
| short | `(10/21)c` | 0 | 1 |
| short | `(10/21)c + (5/21)d` | 1 | 2 |
| long | `(-5/14)c` | 0 | 2 |
| long | `(-5/7)c - (5/7)d` | -3 | 1 |

Consequently the 24 distinct lifts are exactly

```text
{alpha, alpha+delta : alpha short}
    union {beta, beta-3delta : beta long}.
```

Short-root levels are spaced by 1 and long-root levels by 3.  This is the
twisted affine system `D4^(3)`; untwisted `G2^(1)` would require a long-root
level `+1` or `-1`, and no such lift occurs.

A concrete affine simple system is

```text
alpha_0 = (-1,2) + delta,
alpha_1 = (-1,-1),
alpha_2 = (3,0).
```

It is represented by rays `[25, 46]`, `[45]`, and
`[1, 39]`, respectively.  Its Cartan matrix is

```text
[ 2 -1  0]
[-1  2 -3]
[ 0 -1  2],
```

and the Kac-label relation is exactly

```text
alpha_0 + 2 alpha_1 + alpha_2 = delta  modulo Qc.
```

On abstract level coordinates `(x,y,n)`, the three corresponding affine
reflections are

```text
s_0(x,y,n) = (x+y, -y, n-y),
s_1(x,y,n) = (-y, -x, n),
s_2(x,y,n) = (-x-y, y, n).
```

They are involutions with Coxeter orders `m_01=3`, `m_02=2`, `m_12=6`.
Moreover, `(s_0 s_1 s_2)^2` is a nontrivial unipotent translation, so they
generate the expected infinite affine Coxeter group.

The finite fan is only a window into that group.  No realization of the three
simple roots above is cofacial; in fact no level-zero short and level-zero
long lift share a maximal cone.  Also `s_0` preserves only
12 of the 24 displayed real roots.  Hence the finite fan
does not carry the full affine Weyl action, even though its projected rays and
incidence encode an exact `D4^(3)`/affine-G2 fragment.

## Exact certificate

- Rotation: `r = [5, 6, 4, 2, 3, 1, 7]` (taxon `i` maps to `r[i]`).
- Reflection: `s = [1, 3, 2, 4, 6, 5, 7]`.
- Exact relations checked: `order(r)=6`, `order(s)=2`, and `s*r*s=r^(-1)`.
- Generated group: 12 elements, exactly the stored fan-symmetry group.
- Fan action checked from the rational ray and maximal-cone data.

| kind | exponent | order | taxon cycle type | codim Fix in R^21 | codim Fix on fan span/L | codim Fix on R^21/L |
|---|---:|---:|---:|---:|---:|---:|
| identity | 0 | 1 | `1+1+1+1+1+1+1` | 0 | 0 | 0 |
| rotation | 1 | 6 | `6+1` | 17 | 10 | 15 |
| rotation | 2 | 3 | `3+3+1` | 14 | 8 | 12 |
| rotation | 3 | 2 | `2+2+2+1` | 9 | 4 | 7 |
| rotation | 4 | 3 | `3+3+1` | 14 | 8 | 12 |
| rotation | 5 | 6 | `6+1` | 17 | 10 | 15 |
| hexagon reflection | 0 | 2 | `2+2+1+1+1` | 8 | 5 | 7 |
| hexagon reflection | 1 | 2 | `2+2+2+1` | 9 | 5 | 8 |
| hexagon reflection | 2 | 2 | `2+2+1+1+1` | 8 | 5 | 7 |
| hexagon reflection | 3 | 2 | `2+2+2+1` | 9 | 5 | 8 |
| hexagon reflection | 4 | 2 | `2+2+1+1+1` | 8 | 5 | 7 |
| hexagon reflection | 5 | 2 | `2+2+2+1` | 9 | 5 | 8 |

Here `codim Fix in R^21` is `rank(I-g)` on the pair-coordinate permutation
representation.  A genuine linear reflection would have value 1.  The last
two columns perform the analogous exact computations on the 12-dimensional
span of the fan rays modulo lineality and on the whole 18-dimensional ambient
quotient, respectively.

## Reproduce

From the repository root:

```sh
julia --project=julia julia/studies/g2/fan_symmetries/analyze_dihedral_action.jl
```

