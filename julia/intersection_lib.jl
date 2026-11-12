using Oscar
using LinearAlgebra
using Combinatorics
using Random
using ProgressMeter
using HiGHS
using Graphs
using Printf

import Graphs: SimpleGraph, add_vertex!, add_edge!, has_edge,
               nv, ne, connected_components, vertices, degree


"""
Save an Oscar object to a temporary file beside `filename`, then atomically replace the destination. A killed process
can therefore leave at most a disposable temporary file, never a partially written checkpoint at `filename`.
"""
function atomic_oscar_save(filename::AbstractString, value)
    stem, extension = splitext(filename)
    temporary = "$(stem).tmp$(extension)"
    mkpath(dirname(filename))
    rm(temporary; force=true)
    try
        save(temporary, value)
        Base.Filesystem.rename(temporary, filename)
    finally
        rm(temporary; force=true)
    end
    return nothing
end

# Column indices for HiGHS rows (0-based, reused everywhere)

# ══════════════════════════════════════════════════════════════════════════════
# H-representation extraction
# ══════════════════════════════════════════════════════════════════════════════

function primitive(v::Vector{QQFieldElem})::Vector{QQFieldElem}
    k = lcm([vi.den for vi in v])
    return [k*vi for vi in v]
end

"""
    get_hrep_f64(P) -> (A, b)

Extract Float64 H-representation {x : Ax ≤ b} from a Polyhedron,
combining facets and affine hull (as pairs of inequalities).
"""
function get_hrep_f64(P::Polyhedron)
    n    = ambient_dim(P)
    pm_P = Oscar.pm_object(P)

    function to_f64_int_matrix(M::Matrix{QQFieldElem})
        nrows, ncols = size(M)
        isempty(M) && return Matrix{Float64}(undef, 0, ncols)
        result = Matrix{Float64}(undef, nrows, ncols)
        for i in 1:nrows
            result[i, :] = Float64.(primitive(M[i, :]))
        end
        return result
    end

    F = to_f64_int_matrix(Matrix{QQFieldElem}(pm_P.FACETS))
    E = to_f64_int_matrix(Matrix{QQFieldElem}(pm_P.AFFINE_HULL))

    A = isempty(F) ? zeros(0, n) : -F[:, 2:end]
    b = isempty(F) ? Float64[]   :  F[:, 1]

    if !isempty(E)
        A = vcat(A,  E[:, 2:end], -E[:, 2:end])
        b = vcat(b, -E[:, 1],      E[:, 1])
    end

    return A, b
end

"""
    get_hrep_qq(P) -> (A, b)

Extract the exact rational H-representation `{x : Ax ≤ b}` from a
`Polyhedron`, combining facets and affine-hull equations (the latter as pairs
of inequalities). Each row is scaled to have integer coefficients.
"""
function get_hrep_qq(P::Polyhedron)
    n    = ambient_dim(P)
    pm_P = Oscar.pm_object(P)

    function to_qq_int_matrix(M::Matrix{QQFieldElem})
        nrows, ncols = size(M)
        isempty(M) && return Matrix{QQFieldElem}(undef, 0, ncols)
        result = Matrix{QQFieldElem}(undef, nrows, ncols)
        for i in 1:nrows
            result[i, :] = primitive(M[i, :])
        end
        return result
    end

    F = to_qq_int_matrix(Matrix{QQFieldElem}(pm_P.FACETS))
    E = to_qq_int_matrix(Matrix{QQFieldElem}(pm_P.AFFINE_HULL))

    A = isempty(F) ? Matrix{QQFieldElem}(undef, 0, n) : -F[:, 2:end]
    b = isempty(F) ? QQFieldElem[]                    :  F[:, 1]

    if !isempty(E)
        A = vcat(A,  E[:, 2:end], -E[:, 2:end])
        b = vcat(b, -E[:, 1],      E[:, 1])
    end

    return A, b
end

"""
    translate_hrep(A, b, t) -> b_shifted

For a cone with H-rep {x : Ax ≤ b}, the translate C + t has H-rep
{x : A(x-t) ≤ b}, i.e. {x : Ax ≤ b + A*t}.  Returns the shifted RHS.
"""
translate_hrep(A::Matrix{Float64}, b::Vector{Float64}, t::Vector{Float64}) =
    b .+ A * t

translate_hrep(
    A::Matrix{QQFieldElem},
    b::Vector{QQFieldElem},
    t::Vector{QQFieldElem},
) = b .+ A * t

# ══════════════════════════════════════════════════════════════════════════════
# HiGHS direct API — hot-start feasibility
# ══════════════════════════════════════════════════════════════════════════════

const kHighsModelStatusOptimal = HiGHS.kHighsModelStatusOptimal
const kHighsModelStatusInfeasible = HiGHS.kHighsModelStatusInfeasible

"""
    make_base_highs(A_base, b_base) -> h

Create a HiGHS instance with N free variables and the rows of A_base.
The caller is responsible for calling Highs_destroy(h) when done.
"""
function make_base_highs(A_base::Matrix{Float64}, b_base::Vector{Float64})
    h = HiGHS.Highs_create()
    N = size(A_base, 2)
    COL_INDS = Cint.(collect(0:N-1))
    HiGHS.Highs_setBoolOptionValue(h, "output_flag",    false)
    HiGHS.Highs_setBoolOptionValue(h, "log_to_console", false)
    HiGHS.Highs_setIntOptionValue( h, "simplex_strategy", 1)    # dual simplex
    # HiGHS fixes its scheduler size at the first run of a process; the hybrid
    # models use 1 thread too, so both filters can run in one session.
    HiGHS.Highs_setIntOptionValue( h, "threads", 1)
    HiGHS.Highs_setStringOptionValue(h, "presolve", "on")
    for _ in 1:N
        HiGHS.Highs_addVar(h, -1e30, 1e30)
    end
    for i in axes(A_base, 1)
        HiGHS.Highs_addRow(h, -1e30, b_base[i], N, COL_INDS, A_base[i, :])
    end
    return h
end

"""
    test_and_remove!(h, n_base, A2, b2) -> Bool

Add rows of A2, run HiGHS, delete those rows, return feasibility.
`n_base` = number of rows already present in h.
"""
function test_and_remove!(h, n_base::Int,
                          A2::Matrix{Float64}, b2::Vector{Float64})::Bool
    N = size(A2, 2)
    n2= size(A2, 1)
    COL_INDS = Cint.(collect(0:N-1))
    for i in 1:n2
        HiGHS.Highs_addRow(h, -1e30, b2[i], N, COL_INDS, A2[i, :])
    end
    HiGHS.Highs_run(h)
    st = HiGHS.Highs_getModelStatus(h)
    if st != kHighsModelStatusOptimal
        # The basis kept from the previous pair, or presolve with right-hand
        # sides near 1e9, can call a feasible pair infeasible: solve it again
        # from scratch without presolve.
        HiGHS.Highs_clearSolver(h)
        HiGHS.Highs_setStringOptionValue(h, "presolve", "off")
        HiGHS.Highs_run(h)
        st = HiGHS.Highs_getModelStatus(h)
        HiGHS.Highs_setStringOptionValue(h, "presolve", "on")
    end
    HiGHS.Highs_deleteRowsByRange(h, Cint(n_base), Cint(n_base + n2 - 1))
    return st == kHighsModelStatusOptimal
end

# ══════════════════════════════════════════════════════════════════════════════
# Phase 1 — LP pre-filter
# ══════════════════════════════════════════════════════════════════════════════

"""
    lp_prefilter(tl_hreps, tg_hreps, norm_A, norm_b) -> feasible_pairs

Parallel LP pre-filter. For each i, builds a hot-start HiGHS instance
(linear-space cone + normalization rows), then tests all j (Mukai-fan cones).
"""
function lp_prefilter(
    tl_hreps :: Vector{Tuple{Matrix{Float64}, Vector{Float64}}},
    tg_hreps :: Vector{Tuple{Matrix{Float64}, Vector{Float64}}},
    norm_A   :: Matrix{Float64},
    norm_b   :: Vector{Float64},
)::Vector{Tuple{Int,Int}}

    feasible_pairs = Tuple{Int,Int}[]
    lk = ReentrantLock()

    @showprogress Threads.@threads for i in eachindex(tl_hreps)
        A1, b1 = tl_hreps[i]

        # Base model: linear-space cone rows + normalization rows
        A_base = vcat(A1, norm_A)
        b_base = vcat(b1, norm_b)
        n_base = size(A_base, 1)

        h        = make_base_highs(A_base, b_base)
        local_fp = Tuple{Int,Int}[]

        for j in eachindex(tg_hreps)
            A2, b2 = tg_hreps[j]
            test_and_remove!(h, n_base, A2, b2) && push!(local_fp, (i, j))
        end

        HiGHS.Highs_destroy(h)

        if !isempty(local_fp)
            lock(lk) do
                append!(feasible_pairs, local_fp)
            end
        end
    end

    return feasible_pairs
end

"""
    soplex_phase1_is_feasible(A, b) -> Bool

Reference implementation of the exact Phase-I test through Oscar/Polymake.
The production prefilter below calls the SoPlex C API directly instead.
"""
function soplex_phase1_is_feasible(
    A :: Matrix{QQFieldElem},
    b :: Vector{QQFieldElem},
)::Bool
    m, n = size(A)
    length(b) == m || throw(DimensionMismatch(
        "A has $m rows but b has length $(length(b))"
    ))

    # The last variable is the nonnegative artificial slack s. The Phase-I
    # problem is always feasible and bounded below, so SoPlex need only return
    # an exact optimum rather than classify an infeasible zero-objective LP.
    A_phase1 = hcat(A, fill(QQ(-1), m))
    slack_nonnegative = zeros(QQFieldElem, 1, n + 1)
    slack_nonnegative[1, end] = QQ(-1)
    A_phase1 = vcat(A_phase1, slack_nonnegative)
    b_phase1 = vcat(b, QQ(0))

    P = polyhedron(A_phase1, b_phase1)
    objective = zeros(QQFieldElem, n + 1)
    objective[end] = QQ(1)
    phase1_lp = linear_program(P, objective; convention=:min)

    optimum = optimal_value(phase1_lp)
    isnothing(optimum) && error("SoPlex Phase-I LP unexpectedly infeasible")
    return iszero(optimum)
end

# SoPlex is shipped by SCIP_jll, a transitive Oscar dependency. Loading it by
# PkgId avoids requiring users to add SCIP_jll as a direct project dependency.
const _SCIP_JLL = Base.require(Base.PkgId(
    Base.UUID("e5ac4fe4-a920-5659-9bf8-f9f73e9e79ce"), "SCIP_jll"
))
const _SOPLEX_LIBRARY = _SCIP_JLL.libsoplex
const _SOPLEX_INFINITY = 1.0e100
const _SOPLEX_STATUS_OPTIMAL = Cint(1)

mutable struct _SoPlexPhase1Model
    handle       :: Ptr{Cvoid}
    n_variables  :: Int
    n_rows       :: Int
    n_cols       :: Int
    zero_col     :: Vector{Cdouble}
    row_nums     :: Vector{Clong}
    row_denoms   :: Vector{Clong}
end

"""C-`long` encoding of rational rows, stored column-wise by LP row."""
struct _SoPlexRationalRows
    nums          :: Matrix{Clong}
    denoms        :: Matrix{Clong}
    rhs_nums      :: Vector{Clong}
    rhs_denoms    :: Vector{Clong}
    n_nonzeros    :: Vector{Cint}
end

function _make_soplex_phase1_model(
    n_variables::Int,
    max_rows::Int,
)::_SoPlexPhase1Model
    h = ccall((:SoPlex_create, _SOPLEX_LIBRARY), Ptr{Cvoid}, ())
    h == C_NULL && error("SoPlex_create failed")

    try
        # This enables rational read/solve/check modes and automatic
        # synchronization of the floating-point and rational LPs.
        ccall((:SoPlex_setRational, _SOPLEX_LIBRARY), Cvoid,
              (Ptr{Cvoid},), h)

        # IntParam codes from soplex.h: objective sense, verbosity, simplifier.
        ccall((:SoPlex_setIntParam, _SOPLEX_LIBRARY), Cvoid,
              (Ptr{Cvoid}, Cint, Cint), h, 0, -1)  # minimize
        ccall((:SoPlex_setIntParam, _SOPLEX_LIBRARY), Cvoid,
              (Ptr{Cvoid}, Cint, Cint), h, 9, 0)   # errors only
        ccall((:SoPlex_setIntParam, _SOPLEX_LIBRARY), Cvoid,
              (Ptr{Cvoid}, Cint, Cint), h, 10, 1)  # simplifier auto

        max_cols = n_variables + 1 + max_rows
        return _SoPlexPhase1Model(
            h, n_variables, 0, 0,
            zeros(Cdouble, max_rows),
            zeros(Clong, max_cols),
            ones(Clong, max_cols),
        )
    catch
        ccall((:SoPlex_free, _SOPLEX_LIBRARY), Cvoid,
              (Ptr{Cvoid},), h)
        rethrow()
    end
end

function _destroy_soplex!(model::_SoPlexPhase1Model)
    if model.handle != C_NULL
        ccall((:SoPlex_free, _SOPLEX_LIBRARY), Cvoid,
              (Ptr{Cvoid},), model.handle)
        model.handle = C_NULL
    end
    return nothing
end

function _clear_soplex_lp!(model::_SoPlexPhase1Model)
    # Despite its name, SoPlex_clearLPReal does not reliably clear the
    # synchronized rational LP in SoPlex 7.1. Remove entries explicitly; these
    # operations update both representations when SYNCMODE_AUTO is active.
    while model.n_rows > 0
        model.n_rows -= 1
        ccall((:SoPlex_removeRowReal, _SOPLEX_LIBRARY), Cvoid,
              (Ptr{Cvoid}, Cint), model.handle, Cint(model.n_rows))
    end
    while model.n_cols > 0
        model.n_cols -= 1
        ccall((:SoPlex_removeColReal, _SOPLEX_LIBRARY), Cvoid,
              (Ptr{Cvoid}, Cint), model.handle, Cint(model.n_cols))
    end
    return nothing
end

function _add_soplex_zero_col!(
    model::_SoPlexPhase1Model;
    objective::Cdouble=Cdouble(0),
    lower::Cdouble=Cdouble(-_SOPLEX_INFINITY),
    upper::Cdouble=Cdouble(_SOPLEX_INFINITY),
)
    model.n_rows <= length(model.zero_col) || error(
        "internal SoPlex row-buffer capacity exceeded"
    )
    ccall((:SoPlex_addColReal, _SOPLEX_LIBRARY), Cvoid,
          (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Cint,
           Cdouble, Cdouble, Cdouble),
          model.handle, model.zero_col, Cint(model.n_rows), Cint(0),
          objective, lower, upper)
    model.n_cols += 1
    return nothing
end

function _qq_clong_parts(x::QQFieldElem)::Tuple{Clong,Clong}
    try
        return Clong(x.num), Clong(x.den)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError(
            "SoPlex's C API only accepts rational numerators and " *
            "denominators that fit in a C long; coefficient $x is too large"
        ))
    end
end

function _encode_soplex_rows(
    A::Matrix{QQFieldElem},
    b::Vector{QQFieldElem},
)::_SoPlexRationalRows
    n_rows, n_variables = size(A)
    n_rows == length(b) || throw(DimensionMismatch(
        "A has $n_rows rows but b has length $(length(b))"
    ))

    # A column contains one original inequality row, making its coefficients
    # contiguous when copied into the dense SoPlex C row buffer.
    nums = Matrix{Clong}(undef, n_variables, n_rows)
    denoms = Matrix{Clong}(undef, n_variables, n_rows)
    rhs_nums = Vector{Clong}(undef, n_rows)
    rhs_denoms = Vector{Clong}(undef, n_rows)
    n_nonzeros = Vector{Cint}(undef, n_rows)

    for row in 1:n_rows
        nnz = 0
        for col in 1:n_variables
            num, den = _qq_clong_parts(A[row, col])
            nums[col, row] = num
            denoms[col, row] = den
            nnz += !iszero(num)
        end
        rhs_nums[row], rhs_denoms[row] = _qq_clong_parts(b[row])
        n_nonzeros[row] = Cint(nnz)
    end
    return _SoPlexRationalRows(
        nums, denoms, rhs_nums, rhs_denoms, n_nonzeros
    )
end

"""
Add `a*x - phase_slack ≤ rhs` exactly. The rational C row API has no
representation for an infinite row bound, so introduce a nonnegative row
slack and submit the equivalent rational equality

    a*x - phase_slack + row_slack = rhs.
"""
function _add_soplex_phase1_row!(
    model::_SoPlexPhase1Model,
    rows::_SoPlexRationalRows,
    row::Int,
)
    size(rows.nums, 1) == model.n_variables || throw(DimensionMismatch(
        "row has length $(size(rows.nums, 1)); expected $(model.n_variables)"
    ))
    checkbounds(rows.rhs_nums, row)

    # Add this row's nonnegative slack column before constructing the row.
    _add_soplex_zero_col!(
        model;
        lower=Cdouble(0), upper=Cdouble(_SOPLEX_INFINITY),
    )

    model.n_cols <= length(model.row_nums) || error(
        "internal SoPlex column-buffer capacity exceeded"
    )
    n = model.n_variables
    source_offset = (row - 1) * n + 1
    copyto!(model.row_nums, 1, rows.nums, source_offset, n)
    copyto!(model.row_denoms, 1, rows.denoms, source_offset, n)
    for col in n + 1:model.n_cols
        model.row_nums[col] = Clong(0)
        model.row_denoms[col] = Clong(1)
    end

    # Variable n+1 is the shared Phase-I slack; the newest variable is the
    # nonnegative slack belonging to this equality row.
    model.row_nums[n + 1] = Clong(-1)
    model.row_nums[model.n_cols] = Clong(1)
    rhs_num = rows.rhs_nums[row]
    rhs_den = rows.rhs_denoms[row]
    n_nonzero = rows.n_nonzeros[row] + Cint(2)

    ccall((:SoPlex_addRowRational, _SOPLEX_LIBRARY), Cvoid,
          (Ptr{Cvoid}, Ptr{Clong}, Ptr{Clong}, Cint, Cint,
           Clong, Clong, Clong, Clong),
          model.handle, model.row_nums, model.row_denoms,
          Cint(model.n_cols), n_nonzero,
          rhs_num, rhs_den, rhs_num, rhs_den)
    model.n_rows += 1
    return nothing
end

function _add_soplex_phase1_rows!(
    model::_SoPlexPhase1Model,
    rows::_SoPlexRationalRows,
)
    size(rows.nums, 1) == model.n_variables || throw(DimensionMismatch(
        "A has $(size(rows.nums, 1)) columns; expected $(model.n_variables)"
    ))
    for row in eachindex(rows.rhs_nums)
        _add_soplex_phase1_row!(model, rows, row)
    end
    return nothing
end

function _reset_soplex_phase1!(
    model::_SoPlexPhase1Model,
    base_row_sets::_SoPlexRationalRows...,
)
    _clear_soplex_lp!(model)

    # Original variables are free. The final common variable is the
    # nonnegative Phase-I slack with objective coefficient one.
    for _ in 1:model.n_variables
        _add_soplex_zero_col!(model)
    end
    _add_soplex_zero_col!(
        model;
        objective=Cdouble(1), lower=Cdouble(0),
        upper=Cdouble(_SOPLEX_INFINITY),
    )
    for rows in base_row_sets
        _add_soplex_phase1_rows!(model, rows)
    end
    return nothing
end

function _soplex_exact_objective_is_zero(model::_SoPlexPhase1Model)::Bool
    ptr = ccall((:SoPlex_objValueRationalString, _SOPLEX_LIBRARY),
                Ptr{Cchar}, (Ptr{Cvoid},), model.handle)
    ptr == C_NULL && error("SoPlex did not return a rational objective value")
    value = try
        unsafe_string(ptr)
    finally
        Libc.free(ptr)
    end

    # SoPlex returns either an integer or a numerator/denominator string.
    numerator_text = first(split(strip(value), '/'; limit=2))
    return iszero(parse(BigInt, numerator_text))
end

function _test_and_remove_soplex!(
    model::_SoPlexPhase1Model,
    dynamic_rows::_SoPlexRationalRows,
)::Bool
    n_base_rows = model.n_rows
    n_base_cols = model.n_cols

    try
        _add_soplex_phase1_rows!(model, dynamic_rows)
        status = ccall((:SoPlex_optimize, _SOPLEX_LIBRARY), Cint,
                       (Ptr{Cvoid},), model.handle)
        status == _SOPLEX_STATUS_OPTIMAL || error(
            "SoPlex Phase-I solve returned status $status instead of OPTIMAL"
        )
        return _soplex_exact_objective_is_zero(model)
    finally
        # Remove in reverse order. This preserves the base LP and its basis for
        # the next candidate tested by this thread.
        while model.n_rows > n_base_rows
            model.n_rows -= 1
            ccall((:SoPlex_removeRowReal, _SOPLEX_LIBRARY), Cvoid,
                  (Ptr{Cvoid}, Cint), model.handle, Cint(model.n_rows))
        end
        while model.n_cols > n_base_cols
            model.n_cols -= 1
            ccall((:SoPlex_removeColReal, _SOPLEX_LIBRARY), Cvoid,
                  (Ptr{Cvoid}, Cint), model.handle, Cint(model.n_cols))
        end
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Hybrid HiGHS + exact-certificate + SoPlex fallback
# ══════════════════════════════════════════════════════════════════════════════

mutable struct _HybridThreadModel
    highs       :: Ptr{Cvoid}
    soplex      :: _SoPlexPhase1Model
    col_indices :: Vector{Cint}
    primal      :: Vector{Cdouble}
    dual_ray    :: Vector{Cdouble}
end

"""Float64 rows stored as contiguous columns for the HiGHS C row API."""
struct _HighsRowSet
    columns :: Matrix{Float64}
    rhs     :: Vector{Float64}
end

function _encode_highs_rows(
    A::Matrix{QQFieldElem},
    b::Vector{QQFieldElem},
)::_HighsRowSet
    size(A, 1) == length(b) || throw(DimensionMismatch(
        "A has $(size(A, 1)) rows but b has length $(length(b))"
    ))
    columns = permutedims(Float64.(A))
    rhs = Float64.(b)
    all(isfinite, columns) && all(isfinite, rhs) || throw(ArgumentError(
        "exact H-representation cannot be represented finitely in Float64"
    ))
    return _HighsRowSet(columns, rhs)
end

function _make_hybrid_thread_model(
    n_variables::Int,
    max_rows::Int,
)::_HybridThreadModel
    highs = HiGHS.Highs_create()
    highs == C_NULL && error("Highs_create failed")
    try
        HiGHS.Highs_setBoolOptionValue(highs, "output_flag", false)
        HiGHS.Highs_setBoolOptionValue(highs, "log_to_console", false)
        HiGHS.Highs_setIntOptionValue(highs, "simplex_strategy", 1)
        HiGHS.Highs_setIntOptionValue(highs, "threads", 1)
        HiGHS.Highs_setStringOptionValue(highs, "presolve", "on")

        col_indices = Cint.(0:n_variables-1)
        primal = zeros(Cdouble, n_variables)
        dual_ray = zeros(Cdouble, max_rows)
        soplex = _make_soplex_phase1_model(n_variables, max_rows)
        return _HybridThreadModel(
            highs, soplex, col_indices, primal, dual_ray
        )
    catch
        HiGHS.Highs_destroy(highs)
        rethrow()
    end
end

function _destroy_hybrid_model!(model::_HybridThreadModel)
    if model.highs != C_NULL
        HiGHS.Highs_destroy(model.highs)
        model.highs = C_NULL
    end
    _destroy_soplex!(model.soplex)
    return nothing
end

function _add_highs_rows!(
    model::_HybridThreadModel,
    rows::_HighsRowSet,
)
    size(rows.columns, 1) == length(model.col_indices) ||
        throw(DimensionMismatch(
        "rows have $(size(rows.columns, 1)) columns; " *
        "expected $(length(model.col_indices))"
    ))
    size(rows.columns, 2) == length(rows.rhs) || throw(DimensionMismatch(
        "row matrix has $(size(rows.columns, 2)) rows but RHS has " *
        "length $(length(rows.rhs))"
    ))
    highs_infinity = HiGHS.Highs_getInfinity(model.highs)
    for row in eachindex(rows.rhs)
        ret = HiGHS.Highs_addRow(
            model.highs, -highs_infinity, rows.rhs[row],
            size(rows.columns, 1), model.col_indices,
            view(rows.columns, :, row),
        )
        ret == HiGHS.kHighsStatusOk || error(
            "Highs_addRow failed with status $ret"
        )
    end
    return nothing
end

function _reset_hybrid_highs!(
    model::_HybridThreadModel,
    base_row_sets::_HighsRowSet...,
)
    ret = HiGHS.Highs_clearModel(model.highs)
    ret == HiGHS.kHighsStatusOk || error(
        "Highs_clearModel failed with status $ret"
    )
    highs_infinity = HiGHS.Highs_getInfinity(model.highs)
    for _ in model.col_indices
        ret = HiGHS.Highs_addVar(
            model.highs, -highs_infinity, highs_infinity
        )
        ret == HiGHS.kHighsStatusOk || error(
            "Highs_addVar failed with status $ret"
        )
    end
    for rows in base_row_sets
        _add_highs_rows!(model, rows)
    end
    return nothing
end

function _remove_highs_dynamic_rows!(
    model::_HybridThreadModel,
    n_base::Int,
    n_dynamic::Int,
)
    n_dynamic == 0 && return nothing
    ret = HiGHS.Highs_deleteRowsByRange(
        model.highs, Cint(n_base), Cint(n_base + n_dynamic - 1)
    )
    ret == HiGHS.kHighsStatusOk || error(
        "Highs_deleteRowsByRange failed with status $ret"
    )
    return nothing
end

function _float_to_qq(x::Float64, tolerance::Float64)::QQFieldElem
    isfinite(x) || throw(ArgumentError(
        "cannot reconstruct non-finite floating-point value $x over QQ"
    ))
    r = rationalize(BigInt, x; tol=tolerance)
    return QQ(numerator(r), denominator(r))
end

function _satisfies_exactly(
    x::Vector{QQFieldElem},
    row_sets::Tuple{Matrix{QQFieldElem},Vector{QQFieldElem}}...,
)::Bool
    for (A, b) in row_sets
        for row in axes(A, 1)
            lhs = QQ(0)
            for col in axes(A, 2)
                lhs += A[row, col] * x[col]
            end
            lhs <= b[row] || return false
        end
    end
    return true
end

function _project_to_normalization!(
    x::Vector{QQFieldElem},
    norm_A::Matrix{QQFieldElem},
    norm_b::Vector{QQFieldElem},
)::Bool
    size(norm_A, 1) >= 2 || return false
    length(norm_b) >= 2 || return false
    all(norm_A[2, col] == -norm_A[1, col] for col in axes(norm_A, 2)) ||
        return false
    norm_b[2] == -norm_b[1] || return false

    pivot = findlast(!iszero, view(norm_A, 1, :))
    isnothing(pivot) && return false
    lhs = QQ(0)
    for col in axes(norm_A, 2)
        lhs += norm_A[1, col] * x[col]
    end
    x[pivot] += (norm_b[1] - lhs) / norm_A[1, pivot]
    return true
end

"""Try to turn a HiGHS primal point into an exactly verified QQ witness."""
function _highs_primal_is_exactly_feasible(
    primal::AbstractVector{Cdouble},
    norm_A::Matrix{QQFieldElem},
    norm_b::Vector{QQFieldElem},
    row_sets::Tuple{Matrix{QQFieldElem},Vector{QQFieldElem}}...,
)::Bool
    all(isfinite, primal) || return false
    for tolerance in (1.0e-7, 1.0e-9, 1.0e-11, 1.0e-13)
        x = [_float_to_qq(value, tolerance) for value in primal]
        _satisfies_exactly(x, row_sets...) && return true

        # The normalization is an exact equality represented by two opposite
        # inequalities. Repairing one coordinate often recovers the exact
        # rational basic solution from the floating-point point.
        projected = copy(x)
        if _project_to_normalization!(projected, norm_A, norm_b)
            _satisfies_exactly(projected, row_sets...) && return true
        end
    end
    return false
end

function _farkas_certificate_is_valid(
    lambda::Vector{QQFieldElem},
    row_sets::Tuple{Matrix{QQFieldElem},Vector{QQFieldElem}}...,
)::Bool
    total_rows = sum(size(A, 1) for (A, _) in row_sets)
    length(lambda) == total_rows || return false
    all(value -> value >= 0, lambda) || return false

    n_variables = size(first(row_sets)[1], 2)
    products = zeros(QQFieldElem, n_variables)
    rhs_product = QQ(0)
    offset = 0
    for (A, b) in row_sets
        for row in axes(A, 1)
            multiplier = lambda[offset + row]
            if !iszero(multiplier)
                for col in axes(A, 2)
                    products[col] += multiplier * A[row, col]
                end
                rhs_product += multiplier * b[row]
            end
        end
        offset += size(A, 1)
    end
    return all(iszero, products) && rhs_product < 0
end

function _collect_exact_rows(
    row_sets::Tuple{Matrix{QQFieldElem},Vector{QQFieldElem}}...,
)::Tuple{Matrix{QQFieldElem},Vector{QQFieldElem}}
    return (
        reduce(vcat, first.(row_sets)),
        reduce(vcat, last.(row_sets)),
    )
end

"""
Try to reconstruct and exactly verify a Farkas multiplier from a floating
HiGHS dual ray. For upper-bounded rows `A*x <= b`, HiGHS returns the opposite
sign, hence `lambda = -ray`.
"""
function _highs_ray_is_exact_infeasibility_certificate(
    ray::AbstractVector{Cdouble},
    row_sets::Tuple{Matrix{QQFieldElem},Vector{QQFieldElem}}...,
)::Bool
    isempty(ray) && return false
    all(isfinite, ray) || return false
    scale = maximum(abs, ray)
    iszero(scale) && return false

    # Recover the exact kernel on the numerical support of the ray. Simplex
    # normally returns an extreme ray, giving a one-dimensional kernel that can
    # be oriented and checked without any numerical reconstruction.
    A, b = _collect_exact_rows(row_sets...)
    for relative_cutoff in (1.0e-8, 1.0e-10, 1.0e-12)
        support = findall(value -> abs(value) > scale * relative_cutoff, ray)
        isempty(support) && continue
        length(support) > size(A, 2) + 2 && continue

        A_support = A[support, :]
        kernel_dim, kernel_matrix = nullspace(
            matrix(QQ, transpose(A_support))
        )
        kernel_dim >= 1 || continue
        exact_kernel = Matrix{QQFieldElem}(kernel_matrix)

        if kernel_dim == 1
            lambda_support = vec(exact_kernel)
            if all(value -> value <= 0, lambda_support)
                lambda_support .*= QQ(-1)
            end
            all(value -> value >= 0, lambda_support) || continue
            rhs_product = QQ(0)
            for k in eachindex(support)
                rhs_product += b[support[k]] * lambda_support[k]
            end
            rhs_product < 0 && return true
            continue
        end

        # Degenerate rays can have a higher-dimensional support kernel. Project
        # the floating ray onto an exact rational combination of that kernel.
        float_kernel = Float64.(exact_kernel)
        target = Float64.(-ray[support] ./ scale)
        coefficients = float_kernel \ target
        all(isfinite, coefficients) || continue

        for tolerance in (1.0e-8, 1.0e-10, 1.0e-6, 1.0e-12)
            exact_coefficients = [
                _float_to_qq(value, tolerance) for value in coefficients
            ]
            lambda_support = exact_kernel * exact_coefficients
            all(value -> value >= 0, lambda_support) || continue

            lhs = transpose(A_support) * lambda_support
            all(iszero, lhs) || continue
            rhs_product = QQ(0)
            for k in eachindex(support)
                rhs_product += b[support[k]] * lambda_support[k]
            end
            rhs_product < 0 && return true
        end
    end

    # Last cheap attempt: component-wise rational reconstruction. Exact
    # verification below makes any unsuccessful approximation harmless.
    for tolerance in (1.0e-8, 1.0e-10, 1.0e-6, 1.0e-12)
        lambda = [
            _float_to_qq(-value / scale, tolerance)
            for value in ray
        ]
        _farkas_certificate_is_valid(lambda, row_sets...) && return true
    end
    return false
end

"""
    soplex_prefilter(tl_hreps, tg_hreps, norm_A, norm_b) -> feasible_pairs

Parallel exact-rational Phase-I filter using the SoPlex C API directly. Each
Julia thread owns one persistent model. For a fixed `tl_hreps[i]`, its base
rows and the normalization rows are installed once; each `tg_hreps[j]` is
added, solved, and removed in place. The original inequalities are feasible
exactly when the rational Phase-I optimum is zero.
"""
function soplex_prefilter(
    tl_hreps :: Vector{Tuple{Matrix{QQFieldElem}, Vector{QQFieldElem}}},
    tg_hreps :: Vector{Tuple{Matrix{QQFieldElem}, Vector{QQFieldElem}}},
    norm_A   :: Matrix{QQFieldElem},
    norm_b   :: Vector{QQFieldElem},
)::Vector{Tuple{Int,Int}}
    isempty(tl_hreps) && return Tuple{Int,Int}[]
    isempty(tg_hreps) && return Tuple{Int,Int}[]

    n_variables = size(norm_A, 2)
    size(norm_A, 1) == length(norm_b) || throw(DimensionMismatch(
        "norm_A has $(size(norm_A, 1)) rows but norm_b has " *
        "length $(length(norm_b))"
    ))
    for (A, b) in Iterators.flatten((tl_hreps, tg_hreps))
        size(A, 2) == n_variables || throw(DimensionMismatch(
            "H-representation has $(size(A, 2)) columns; expected $n_variables"
        ))
        size(A, 1) == length(b) || throw(DimensionMismatch(
            "H-representation has $(size(A, 1)) rows but RHS has " *
            "length $(length(b))"
        ))
    end

    max_base_rows = maximum(size(A, 1) for (A, _) in tl_hreps) +
                    size(norm_A, 1)
    max_dynamic_rows = maximum(size(A, 1) for (A, _) in tg_hreps)
    max_rows = max_base_rows + max_dynamic_rows

    # TG and normalization rows are reused for every outer cone. Convert their
    # QQ coefficients to the C representation only once per experiment.
    tg_rows = [_encode_soplex_rows(A, b) for (A, b) in tg_hreps]
    normalization_rows = _encode_soplex_rows(norm_A, norm_b)

    # With :static scheduling, every thread uses only its own persistent model.
    models = _SoPlexPhase1Model[]
    try
        # `threadid()` can exceed `nthreads()` when Julia has an interactive
        # thread pool, so size this table by the largest possible thread id.
        for _ in 1:Threads.maxthreadid()
            push!(models, _make_soplex_phase1_model(n_variables, max_rows))
        end

        feasible_by_i = [Tuple{Int,Int}[] for _ in eachindex(tl_hreps)]
        @showprogress Threads.@threads :static for i in eachindex(tl_hreps)
            model = models[Threads.threadid()]
            A1, b1 = tl_hreps[i]
            base_rows = _encode_soplex_rows(A1, b1)
            _reset_soplex_phase1!(
                model, base_rows, normalization_rows
            )

            local_pairs = feasible_by_i[i]
            for j in eachindex(tg_rows)
                _test_and_remove_soplex!(model, tg_rows[j]) &&
                    push!(local_pairs, (i, j))
            end
        end
        return reduce(vcat, feasible_by_i; init=Tuple{Int,Int}[])
    finally
        foreach(_destroy_soplex!, models)
    end
end

"""
    hybrid_prefilter(tl_hreps, tg_hreps, norm_A, norm_b) -> feasible_pairs

Parallel exact feasibility filter. HiGHS first solves each floating-point
system. Its result is used directly only after an exact QQ witness check:

  * an optimal primal point must reconstruct to a point satisfying every row;
  * an infeasible dual ray must reconstruct to a valid Farkas certificate.

All other cases are sent to the persistent per-thread exact SoPlex Phase-I
model. Consequently the returned pair list is exact, while most easy systems
avoid an exact SoPlex solve. With `soplex_fallback=false` these pairs are
returned as feasible instead, so that the list contains every feasible pair
and the caller settles the uncertified ones exactly.
"""
function hybrid_prefilter(
    tl_hreps :: Vector{Tuple{Matrix{QQFieldElem}, Vector{QQFieldElem}}},
    tg_hreps :: Vector{Tuple{Matrix{QQFieldElem}, Vector{QQFieldElem}}},
    norm_A   :: Matrix{QQFieldElem},
    norm_b   :: Vector{QQFieldElem};
    soplex_fallback :: Bool = true,
)::Vector{Tuple{Int,Int}}
    isempty(tl_hreps) && return Tuple{Int,Int}[]
    isempty(tg_hreps) && return Tuple{Int,Int}[]

    n_variables = size(norm_A, 2)
    size(norm_A, 1) == length(norm_b) || throw(DimensionMismatch(
        "norm_A has $(size(norm_A, 1)) rows but norm_b has " *
        "length $(length(norm_b))"
    ))
    for (A, b) in Iterators.flatten((tl_hreps, tg_hreps))
        size(A, 2) == n_variables || throw(DimensionMismatch(
            "H-representation has $(size(A, 2)) columns; expected $n_variables"
        ))
        size(A, 1) == length(b) || throw(DimensionMismatch(
            "H-representation has $(size(A, 1)) rows but RHS has " *
            "length $(length(b))"
        ))
    end

    max_base_rows = maximum(size(A, 1) for (A, _) in tl_hreps) +
                    size(norm_A, 1)
    max_dynamic_rows = maximum(size(A, 1) for (A, _) in tg_hreps)
    max_rows = max_base_rows + max_dynamic_rows

    # Reused by every outer cone and every thread.
    tg_rows_exact = [_encode_soplex_rows(A, b) for (A, b) in tg_hreps]
    norm_rows_exact = _encode_soplex_rows(norm_A, norm_b)
    tg_rows_float = [_encode_highs_rows(A, b) for (A, b) in tg_hreps]
    norm_rows_float = _encode_highs_rows(norm_A, norm_b)

    models = _HybridThreadModel[]
    try
        for _ in 1:Threads.maxthreadid()
            push!(models, _make_hybrid_thread_model(n_variables, max_rows))
        end

        feasible_by_i = [Tuple{Int,Int}[] for _ in eachindex(tl_hreps)]

        # Rows: HiGHS optimal, HiGHS infeasible, other HiGHS status,
        # exact primal witness, exact Farkas certificate, SoPlex fallback,
        # feasible SoPlex fallback.
        stats_by_i = zeros(Int, 7, length(tl_hreps))

        @showprogress Threads.@threads :static for i in eachindex(tl_hreps)
            model = models[Threads.threadid()]
            A1, b1 = tl_hreps[i]
            base_exact = _encode_soplex_rows(A1, b1)
            base_float = _encode_highs_rows(A1, b1)
            _reset_hybrid_highs!(model, base_float, norm_rows_float)

            n_base = size(A1, 1) + size(norm_A, 1)
            soplex_base_ready = false
            local_pairs = feasible_by_i[i]

            for j in eachindex(tg_hreps)
                A2, b2 = tg_hreps[j]
                n_dynamic = size(A2, 1)
                _add_highs_rows!(model, tg_rows_float[j])

                pair_is_feasible = false
                needs_soplex = false
                try
                    run_status = HiGHS.Highs_run(model.highs)
                    model_status = HiGHS.Highs_getModelStatus(model.highs)
                    exact_system = ((A1, b1), (norm_A, norm_b), (A2, b2))

                    if run_status == HiGHS.kHighsStatusOk &&
                       model_status == kHighsModelStatusOptimal
                        stats_by_i[1, i] += 1
                        solution_status = HiGHS.Highs_getSolution(
                            model.highs, model.primal,
                            C_NULL, C_NULL, C_NULL,
                        )
                        if solution_status == HiGHS.kHighsStatusOk &&
                           _highs_primal_is_exactly_feasible(
                               model.primal, norm_A, norm_b, exact_system...
                           )
                            stats_by_i[4, i] += 1
                            pair_is_feasible = true
                        else
                            needs_soplex = true
                        end

                    elseif run_status == HiGHS.kHighsStatusOk &&
                           model_status == kHighsModelStatusInfeasible
                        stats_by_i[2, i] += 1
                        has_ray = Ref{HiGHS.HighsInt}(0)
                        ray_status = HiGHS.Highs_getDualRay(
                            model.highs, has_ray, model.dual_ray
                        )
                        ray = view(model.dual_ray, 1:n_base + n_dynamic)
                        if ray_status != HiGHS.kHighsStatusError &&
                           has_ray[] == 1 &&
                           _highs_ray_is_exact_infeasibility_certificate(
                               ray, exact_system...
                           )
                            stats_by_i[5, i] += 1
                            pair_is_feasible = false
                        else
                            needs_soplex = true
                        end
                    else
                        stats_by_i[3, i] += 1
                        needs_soplex = true
                    end

                    if needs_soplex && !soplex_fallback
                        stats_by_i[6, i] += 1
                        pair_is_feasible = true
                    elseif needs_soplex
                        stats_by_i[6, i] += 1
                        if !soplex_base_ready
                            _reset_soplex_phase1!(
                                model.soplex, base_exact, norm_rows_exact
                            )
                            soplex_base_ready = true
                        end
                        pair_is_feasible = _test_and_remove_soplex!(
                            model.soplex, tg_rows_exact[j]
                        )
                        pair_is_feasible && (stats_by_i[7, i] += 1)
                    end

                    pair_is_feasible && push!(local_pairs, (i, j))
                finally
                    _remove_highs_dynamic_rows!(model, n_base, n_dynamic)
                end
            end
        end

        totals = vec(sum(stats_by_i; dims=2))
        @printf(
            "    hybrid routing: HiGHS optimal=%d (%d exact witnesses), infeasible=%d (%d exact rays), other=%d; uncertified=%d (%d feasible by SoPlex)\n",
            totals[1], totals[4], totals[2], totals[5], totals[3],
            totals[6], totals[7],
        )
        return reduce(vcat, feasible_by_i; init=Tuple{Int,Int}[])
    finally
        foreach(_destroy_hybrid_model!, models)
    end
end

# ══════════════════════════════════════════════════════════════════════════════
# Phase 2 — Exact Oscar intersections → graph
# ══════════════════════════════════════════════════════════════════════════════

struct GraphResult
    g            :: SimpleGraph
    list_P       :: Vector{Polyhedron{QQFieldElem}}
    node_coords  :: Vector{Vector{QQFieldElem}}
    edge_sources :: Dict{Tuple{Int,Int}, Tuple{Int,Int}}
    n_isolated   :: Int    # vertices from dim-0 intersections
    n_fp         :: Int    # LP false positives (dim = -1)
    n_skipped    :: Int    # dim ≥ 2 or rays
end

"""
    build_intersection_graph(feasible_pairs, tl_translated, tg_polys, norm_poly)
        -> GraphResult

Exact QQ intersections for every feasible (i,j) pair.
"""
function build_intersection_graph(
    feasible_pairs  :: Vector{Tuple{Int,Int}},
    tl_translated   :: Dict{Int, Polyhedron},
    tg_polys        :: Vector{<:Polyhedron},
    norm_poly       :: Polyhedron,
)::GraphResult

    g             = SimpleGraph()
    node_coords   = Vector{Vector{QQFieldElem}}()
    coord_to_node = Dict{Vector{QQFieldElem}, Int}()
    edge_sources  = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()

    cnt_isolated = Ref(0)
    cnt_fp       = Ref(0)
    cnt_skipped  = Ref(0)

    function get_or_add!(v)
        key = Vector{QQFieldElem}(collect(v))
        get!(coord_to_node, key) do
            add_vertex!(g)
            id = nv(g)
            push!(node_coords, key)
            id
        end
    end
    list_polyhedra = Polyhedron{QQFieldElem}[]

    @showprogress for (i, j) in feasible_pairs
        inter = intersect(intersect(tl_translated[i], tg_polys[j]), norm_poly)

        d     = dim(inter)
        push!(list_polyhedra, inter)

        if d == 0
            verts = collect(Oscar.vertices(inter))
            isempty(verts) ? (cnt_fp[] += 1; continue) : nothing
            get_or_add!(verts[1])
            cnt_isolated[] += 1

        elseif d == 1
            verts = collect(Oscar.vertices(inter))
            if length(verts) < 2
                cnt_skipped[] += 1; continue
            end
            n1 = get_or_add!(verts[1])
            n2 = get_or_add!(verts[2])
            if !has_edge(g, n1, n2)
                add_edge!(g, n1, n2)
                edge_sources[(min(n1,n2), max(n1,n2))] = (i, j)
            end

        elseif d == -1
            cnt_fp[] += 1

        else
            println("Skipping intersection with dim = $d (i=$i, j=$j)")
            cnt_skipped[] += 1
        end
    end
    return GraphResult(g, list_polyhedra, node_coords, edge_sources,
                       cnt_isolated[], cnt_fp[], cnt_skipped[])
end

# ══════════════════════════════════════════════════════════════════════════════
# Single experiment
# ══════════════════════════════════════════════════════════════════════════════

struct ExperimentResult
    t_int          :: Vector{Int}
    S              :: Int
    n_feasible_lp  :: Int
    n_nodes        :: Int
    n_edges        :: Int
    n_components   :: Int
    n_isolated     :: Int
    n_fp           :: Int
    n_skipped      :: Int
end

"`r` as a NamedTuple, which Oscar saves; `load_experiment_result` reads it back."
experiment_record(r::ExperimentResult) =
    NamedTuple{fieldnames(ExperimentResult)}(ntuple(i -> getfield(r, i), fieldcount(ExperimentResult)))

"The `ExperimentResult` saved as an `experiment_record` in `path`."
function load_experiment_result(path::AbstractString)
    record = load(path)
    record isa NamedTuple && keys(record) == fieldnames(ExperimentResult) ||
        error("$path is not a record of ExperimentResult")
    return ExperimentResult(record...)
end

"""
    run_experiment(
        experiment_id, t_int, tl_polys, tg_polys, tl_hreps_base, tg_hreps;
        output_dir,
        feasibility_solver=:highs,
    ) -> ExperimentResult

Run one full experiment for translation vector t_int:
  1. Build shifted H-reps for the linear space (b_shifted = b + A*t).
  2. Feasibility pre-filter using `:highs` (fast, floating point), `:soplex`
     (parallel exact rational arithmetic), `:hybrid` (HiGHS with exact
     witnesses and persistent SoPlex fallback), or `:certified` (the same
     without SoPlex: the pairs HiGHS cannot certify are kept, so no feasible
     pair is lost).
  3. Exact Oscar intersections, which discard the infeasible pairs kept.
Returns summary statistics and saves the graph, the polyhedra and the result.
"""

function run_experiment(
    experiment_id  :: Int,
    t_int          :: Vector{Int},
    tl_polys       :: Vector,
    tg_polys       :: Vector,
    tl_hreps_base  :: Vector,
    tg_hreps       :: Vector;
    output_dir::AbstractString,
    feasibility_solver::Symbol = :highs,
)::ExperimentResult
    feasibility_solver in (:highs, :soplex, :hybrid, :certified) || throw(ArgumentError(
        "feasibility_solver must be :highs, :soplex, :hybrid, or :certified; " *
        "got $feasibility_solver"
    ))

    N = length(t_int)
    S = sum(t_int)

    # Normalization hyperplane ∑xᵢ = S
    norm_poly = polyhedron(
        matrix(QQ, [ones(Int, 1, N); -ones(Int, 1, N)]),
        [QQ(S), QQ(-S)]
    )

    # Phase 1: feasibility pre-filter
    feasible_pairs = if feasibility_solver === :highs
        float_hrep_type = Tuple{Matrix{Float64}, Vector{Float64}}
        all(h -> h isa float_hrep_type, tl_hreps_base) || throw(ArgumentError(
            ":highs requires Float64 H-representations from get_hrep_f64"
        ))
        all(h -> h isa float_hrep_type, tg_hreps) || throw(ArgumentError(
            ":highs requires Float64 H-representations from get_hrep_f64"
        ))

        t_f64 = Float64.(t_int)
        norm_A = vcat(ones(Float64, 1, N), -ones(Float64, 1, N))
        norm_b = Float64[S; -S]
        tl_hreps_shifted = [
            (A, translate_hrep(A, b, t_f64))
            for (A, b) in tl_hreps_base
        ]

        print("  [HiGHS filter]  ")
        lp_prefilter(tl_hreps_shifted, tg_hreps, norm_A, norm_b)
    else
        rational_hrep_type = Tuple{Matrix{QQFieldElem}, Vector{QQFieldElem}}
        all(h -> h isa rational_hrep_type, tl_hreps_base) || throw(ArgumentError(
            ":$feasibility_solver requires exact H-representations from " *
            "get_hrep_qq"
        ))
        all(h -> h isa rational_hrep_type, tg_hreps) || throw(ArgumentError(
            ":$feasibility_solver requires exact H-representations from " *
            "get_hrep_qq"
        ))

        t_qq = QQ.(t_int)
        norm_A = QQ.(vcat(ones(Int, 1, N), -ones(Int, 1, N)))
        norm_b = QQFieldElem[QQ(S), QQ(-S)]
        tl_hreps_shifted = [
            (A, translate_hrep(A, b, t_qq))
            for (A, b) in tl_hreps_base
        ]

        if feasibility_solver === :soplex
            print("  [SoPlex exact filter]  ")
            soplex_prefilter(tl_hreps_shifted, tg_hreps, norm_A, norm_b)
        elseif feasibility_solver === :certified
            print("  [HiGHS, uncertified pairs kept]  ")
            hybrid_prefilter(tl_hreps_shifted, tg_hreps, norm_A, norm_b; soplex_fallback=false)
        else
            print("  [HiGHS + SoPlex hybrid]  ")
            hybrid_prefilter(tl_hreps_shifted, tg_hreps, norm_A, norm_b)
        end
    end
    println("  → $(length(feasible_pairs)) feasible pairs")

    # Translate linear-space cones exactly (QQ), only for feasible pairs.
    pt_poly     = convex_hull(matrix(QQ, reshape(t_int, 1, N)))
    needed_l7   = sort!(unique(first.(feasible_pairs)))
    linear_space_exact = Dict{Int, Polyhedron}(
        i => minkowski_sum(tl_polys[i], pt_poly)
        for i in needed_l7
    )

    # Phase 2: exact intersections
    print("  [Oscar intersect] ")
    result = build_intersection_graph(
        feasible_pairs,
        linear_space_exact,
        tg_polys,
        norm_poly,
    )
    g  = result.g
    cc = length(connected_components(g))

    experiment_result = ExperimentResult(
        t_int, S,
        length(feasible_pairs),
        nv(g), ne(g), cc,
        result.n_isolated, result.n_fp, result.n_skipped,
    )

    # Write the result marker last. Its presence means that both larger output
    # files were completely written and this experiment can be resumed safely.
    try
        # The SimpleGraph as an Oscar graph, which Oscar saves.
        atomic_oscar_save(
            joinpath(output_dir, "graph_$(experiment_id).mrdi"),
            graph_from_edges(Undirected, [[Graphs.src(e), Graphs.dst(e)] for e in Graphs.edges(g)], nv(g)),
        )
        atomic_oscar_save(
            joinpath(output_dir, "list_polyhedra_$(experiment_id).mrdi"),
            result.list_P,
        )
        atomic_oscar_save(
            joinpath(output_dir, "result_$(experiment_id).mrdi"),
            experiment_record(experiment_result),
        )
    finally
        # Drop the largest Julia-side reference collections before returning.
        # The runner performs a full GC after this function's stack is gone.
        empty!(result.list_P)
        empty!(result.node_coords)
        empty!(result.edge_sources)
        empty!(linear_space_exact)
        empty!(feasible_pairs)
    end

    return experiment_result
end

function intersect_polyhedral_complexes_shifted(
    TG::PolyhedralComplex{QQFieldElem},
    TL::PolyhedralComplex{QQFieldElem},
    t_int::Vector{Int},
    )::Vector{Polyhedron{QQFieldElem}}
    t_f64 = Float64.(t_int)

    N = length(t_f64)
    S = sum(t_int)

    # constructing polytope for fixed hyperplane ∑xᵢ = S
    norm_A = vcat( ones(Float64, 1, N), -ones(Float64, 1, N))
    norm_b = Float64[S; -S]
    norm_poly = polyhedron(
        matrix(QQ, [ones(Int, 1, N); -ones(Int, 1, N)]),
        [QQ(S), QQ(-S)]
    )

    tg_polys = collect(maximal_polyhedra(TG))
    tl_polys = collect(maximal_polyhedra(TL))

    println("Extracting $(length(tl_polys)) tropical linear space H-reps (base, un-shifted)...")
    tl_hreps_base = Vector{Tuple{Matrix{Float64}, Vector{Float64}}}(
        undef, length(tl_polys)
    )
    @showprogress for i in eachindex(tl_polys)
        local_poly=tl_polys[i]
        tl_hreps_base[i] = get_hrep_f64(local_poly)
    end

    println("Extracting $(length(tg_polys)) tropical grassmanian H-reps...")
    tg_hreps = [get_hrep_f64(P) for P in tg_polys]
    println("$(length(tl_hreps_base)) × $(length(tg_hreps)) = $(length(tl_hreps_base)*length(tg_hreps)) pairs per experiment\n")


    tl_hreps_shifted = [
        (A, translate_hrep(A, b, t_f64))
        for (A, b) in tl_hreps_base
    ]

    feasible_pairs = lp_prefilter(tl_hreps_shifted, tg_hreps, norm_A, norm_b)

    pt_poly     = convex_hull(matrix(QQ, reshape(t_int, 1, N)))

    needed_tl   = sort!(unique(first.(feasible_pairs)))
    tl_exact    = Dict{Int, Polyhedron}(
        i => minkowski_sum(tl_polys[i], pt_poly)
        for i in needed_tl
    )

    result = build_intersection_graph(feasible_pairs, tl_exact, tg_polys, norm_poly)
    return result.list_P
end
