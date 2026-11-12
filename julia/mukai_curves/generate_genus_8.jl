using Oscar

const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "data", "mukai_curves"))

# ── 1. Ambient ring: Plücker-type variables ──
R, (x12, x13, x14, x15, x16, 
    x23, x24, x25, x26, 
    x34, x35, x36, 
    x45, x46, 
    x56, 
    ) = polynomial_ring(QQ, [
        "x12","x13","x14","x15","x16",
        "x23","x24","x25","x26",
        "x34","x35","x36",
        "x45","x46",
        "x56",
        ])



# ── 7×7 skew-symmetric matrix ──────────────────────────────────────────────────
X = R[  0    x12  x13  x14  x15  x16   ;
       -x12   0   x23  x24  x25  x26   ;
       -x13 -x23   0   x34  x35  x36   ;
       -x14 -x24 -x34   0   x45  x46   ;
       -x15 -x25 -x35 -x45   0   x56  ;
       -x16 -x26 -x36 -x46 -x56   0   ]

G = ideal(pfaffians(X,4))

TGrass = tropical_variety(G;skip_decomposition=true)

mkpath(OUTPUT_DIR)
save(joinpath(OUTPUT_DIR, "genus_8_fan.mrdi"), polyhedral_complex(TGrass))
