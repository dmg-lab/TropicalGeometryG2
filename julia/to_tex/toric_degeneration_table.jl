###############################################################################
# The paper table of the Newton--Okounkov polytopes of the toric degenerations
# of the G2-Grassmannian, in characteristic 0 and in characteristic 2.
#
# Run from julia/ with `julia --project=.`. The script reads the stage-6
# censuses of both fields and takes the orbit labels, their order, and the
# matching of the orbits of both fields from the G2-tree atlas, by
# studies/g2/atlas/flag/cell_orbit_labels.jl. It prints the tabular, or
# writes it to --output=FILE; main.tex supplies the table float and its
# caption, and loads booktabs for its rules.
###############################################################################

# The censuses, the atlas labels and their numbers, unless a script that
# includes this file has included them already.
isdefined(@__MODULE__, :toric_table_support_rows) || include(joinpath(
    @__DIR__, "..", "studies", "g2", "atlas", "flag", "cell_orbit_labels.jl",
))

const TORIC_TABLE_CHARACTERISTIC = (qq=0, gf2=2)
# The color of the orbit labels, which are bold old-style numbers, as in the
# table of maximal trees of toric_degeneration_tree_table.jl and \celltype in
# main.tex: black.
const TORIC_TABLE_LABEL_COLOR = "000000"
# The table has the booktabs form of the tables of main.tex: the \arraystretch
# of its rows, and the width that its column gaps stretch to fill. The stretch
# is larger than in the other tables, so that the rows whose polytope takes two
# lines stand apart: it spaces the rows, not the lines within an entry. (1.1
# since 2026-11-14, as main.tex has it.)
const TORIC_TABLE_STRETCH = 1.1
# Entries of the type column shortened in main.tex (2026-11-14), so that they fit
# on one line: the phrase of toric_table_type_tex => the printed entry.
const TORIC_TABLE_SHORTENED = Dict(
    "Wedge at a triangular face of the pyramid over the multiplex \$M^{3,5}\$" =>
        "Wedge at a triangular face of the pyramid over \$M^{3,5}\$",
)
const TORIC_TABLE_WIDTH = "\\linewidth"

# The census marks the Chevalley polytope of G2/P2; the atlas build certifies
# that it is unimodularly equivalent to this lecture hall simplex.
const TORIC_TABLE_CHEVALLEY = "Lecture hall simplex \$P_{(1,3,2,3,1)}\$, Chevalley polytope"

"""The argument of `name(...)` if `description` has this form."""
function toric_table_argument(description, name)
    startswith(description, name * "(") && endswith(description, ")") || return nothing
    return chop(description; head=length(name) + 1, tail=1)
end

"""The base and the wedged face of the argument `base; face=face` of a wedge."""
function toric_table_wedge(argument)
    depth = 0
    for index in eachindex(argument)
        character = argument[index]
        depth += character == '(' ? 1 : character == ')' ? -1 : 0
        rest = SubString(argument, index)
        depth == 0 && startswith(rest, "; face=") || continue
        return (String(SubString(argument, 1, prevind(argument, index))),
                String(chop(rest; head=length("; face="), tail=0)))
    end
    error("wedge without a face: $argument")
end

"""The f-vector of the base of a pyramid with f-vector `f_vector`."""
function toric_table_pyramid_base(f_vector)
    f_vector === nothing && return nothing
    base = Int[]
    for (index, faces) in enumerate(f_vector[1:(end - 1)])
        push!(base, faces - (index == 1 ? 1 : base[end]))
    end
    base[end] + 1 == f_vector[end] || error("not the f-vector of a pyramid: $f_vector")
    return base
end

toric_table_f_vector_tex(values) = "\$(" * join(values, ",") * ")\$"

"""The face at which a wedge is taken, in words."""
function toric_table_face(description)
    description == "point" && return "a vertex"
    description == "segment" && return "an edge"
    description == "3-gon" && return "a triangular face"
    description == "4-gon" && return "a quadrilateral face"
    return "a face that is " * toric_table_phrase(description)
end

"""
The polytope type of a description of `RecursivePolytopeTypes` in words, as a
noun phrase with its article. `f_vector` is the f-vector of the polytope when
it is known; it then names an unrecognized base below pyramids.
"""
function toric_table_phrase(description, f_vector=nothing)
    depth = 0
    while (base = toric_table_argument(description, "Pyramid")) !== nothing
        description = base
        f_vector = toric_table_pyramid_base(f_vector)
        depth += 1
    end
    if depth > 0
        name = depth == 1 ? "pyramid" : depth == 2 ? "double pyramid" :
               depth == 3 ? "triple pyramid" : "$depth-fold pyramid"
        return "the $name over " * toric_table_phrase(description, f_vector)
    end
    wedge = toric_table_argument(description, "Wedge")
    if wedge !== nothing
        base, face = toric_table_wedge(wedge)
        return "the wedge at " * toric_table_face(face) * " of " * toric_table_phrase(base)
    end
    prism = toric_table_argument(description, "Prism")
    prism === nothing || return "the prism over " * toric_table_phrase(prism)
    dual = toric_table_argument(description, "Dual")
    dual === nothing || return "the polar dual of " * toric_table_phrase(dual)
    unknown = match(r"^Unrecognized\((\d+)D; (\d+) vertices, (\d+) facets\)$", description)
    if unknown !== nothing
        f_vector === nothing &&
            return "a $(unknown[1])-polytope with $(unknown[2]) vertices and $(unknown[3]) facets"
        string.((length(f_vector), f_vector[1], f_vector[end])) == Tuple(unknown.captures) ||
            error("f-vector $f_vector does not fit $description")
        # The braces keep the formula on one line.
        return "a $(unknown[1])-polytope with \${f=(" * join(f_vector, ",") * ")}\$"
    end
    description == "3-gon" && return "a triangle"
    description == "4-gon" && return "a quadrilateral"
    polygon = match(r"^(\d+)-gon$", description)
    polygon === nothing || return "a \$$(polygon[1])\$-gon"
    simplex = match(r"^(\d+)-simplex$", description)
    simplex === nothing || return "a \$$(simplex[1])\$-simplex"
    multiplex = match(r"^multiplex M\((\d+), (\d+)\)$", description)
    multiplex === nothing || return "the multiplex \$M^{$(multiplex[1]),$(multiplex[2])}\$"
    description == "4-gonal antiprism" && return "a square antiprism"
    antiprism = match(r"^(\d+)-gonal antiprism$", description)
    antiprism === nothing || return "a \$$(antiprism[1])\$-gonal antiprism"
    name = replace(description, "Δ" => "\$\\Delta\$", "×" => "\$\\times\$", "_" => " ")
    return (lowercase(first(name)) in "aeiou" ? "an " : "a ") * name
end

"""The entry of the type column: the phrase of the whole polytope, capitalized."""
function toric_table_type_tex(description, f_vector, chevalley::Bool)
    simplex = occursin(r"^\d+-simplex$", description)
    chevalley && !simplex && error("the Chevalley polytope is recorded as $description")
    chevalley && return TORIC_TABLE_CHEVALLEY
    simplex && return "Simplex"
    phrase = toric_table_phrase(description, f_vector)
    for article in ("the ", "an ", "a ")
        startswith(phrase, article) && (phrase = chopprefix(phrase, article); break)
    end
    return get(TORIC_TABLE_SHORTENED, uppercasefirst(phrase), uppercasefirst(phrase))
end

"""
The toric support orbits of both fields in atlas order: the rows of
`toric_table_support_rows`, each with the data of its polytope, which must
agree in the fields in which the orbit is toric.
"""
function toric_degeneration_rows()
    census = toric_table_census()
    return map(toric_table_support_rows(census)) do row
        records = unique!([
            (
                f_vector=parse.(Int, split(census[field][orbit]["body_f_vector"], ';')),
                lattice_points=parse(Int, census[field][orbit]["body_lattice_points"]),
                type=census[field][orbit]["body_type"],
                chevalley=census[field][orbit]["chevalley_equivalent"] == "true",
            ) for (field, orbit) in row.orbits
        ])
        length(records) == 1 ||
            error("the polytope of atlas label $(row.label) differs between the fields")
        merge(row, only(records))
    end
end

"""
    toric_degeneration_table_tex()

The input-ready tabular: one row per toric support orbit, in the order of
the G2-tree atlas and numbered by its labels followed by the orbit size in
parentheses, as in the tables of trees, and the characteristics in which the
orbit occurs, or "all" if it occurs in both.
"""
function toric_degeneration_table_tex()
    rows = toric_degeneration_rows()
    io = IOBuffer()
    println(io, "% Generated by julia/to_tex/toric_degeneration_table.jl from")
    for path in (TORIC_TABLE_CENSUS..., TORIC_TABLE_ATLAS)
        println(io, "%   $path")
    end
    for field in TORIC_TABLE_FIELDS
        sizes = [row.orbit_size for row in rows if any(==(field) ∘ first, row.orbits)]
        println(io, "% Characteristic $(TORIC_TABLE_CHARACTERISTIC[field]) ",
                "($(TORIC_TABLE_FIELD_NAMES[field])): $(length(sizes)) orbits, ",
                "$(sum(sizes; init=0)) toric degenerations.")
    end
    # Five columns fill the text width of main.tex with a smaller column gap,
    # which stretches to the width. The type column is ragged: its phrases are
    # short for justified lines. The heading of the column of numbers is centred.
    println(io, "\\begingroup\\setlength{\\tabcolsep}{2pt}%")
    println(io, "\\renewcommand{\\arraystretch}{$TORIC_TABLE_STRETCH}%")
    println(io, "\\begin{tabular*}{$TORIC_TABLE_WIDTH}{@{\\extracolsep{\\fill}}",
            "l c >{\\raggedright\\arraybackslash}p{9.7cm} r c@{}}")
    println(io, "  \\toprule")
    println(io, "  Type & \$f\$-vector & Newton--Okounkov polytope & ",
            "\\multicolumn{1}{c}{\$N\$} & Char \\tabularnewline")
    println(io, "  \\midrule")
    for row in rows
        fields = [field for (field, _) in row.orbits]
        characteristics = issetequal(fields, TORIC_TABLE_FIELDS) ? "all" :
            "\$" * join((TORIC_TABLE_CHARACTERISTIC[field] for field in fields), ",") * "\$"
        source = join(("$(TORIC_TABLE_FIELD_NAMES[field]) orbit $orbit"
                       for (field, orbit) in row.orbits), ", ")
        println(io, "  ", join((
            "\\textcolor[HTML]{$TORIC_TABLE_LABEL_COLOR}{\\textbf{\\oldstylenums{$(toric_table_label_number(row.label))}}($(row.orbit_size))}",
            toric_table_f_vector_tex(row.f_vector),
            toric_table_type_tex(row.type, row.f_vector, row.chevalley),
            row.lattice_points,
            characteristics,
        ), " & "), " \\tabularnewline % $source")
    end
    println(io, "  \\bottomrule")
    println(io, "\\end{tabular*}%")
    println(io, "\\endgroup")
    return String(take!(io))
end

"""Write the input-ready tabular to `path`."""
function write_toric_degeneration_table_tex(path::AbstractString)
    isempty(path) && throw(ArgumentError("output path must not be empty"))
    tex = toric_degeneration_table_tex()
    mkpath(dirname(abspath(path)))
    write(path, tex)
    return path
end

function toric_degeneration_table_main(args=ARGS)
    if "--help" in args
        println("Usage: julia --project=. to_tex/toric_degeneration_table.jl " *
                "[--output=FILE]")
        return nothing
    end
    for arg in args
        startswith(arg, "--output=") ||
            throw(ArgumentError("unknown option $arg; use --help"))
    end
    outputs = [arg for arg in args if startswith(arg, "--output=")]
    length(outputs) <= 1 || throw(ArgumentError("pass only one --output"))
    if isempty(outputs)
        print(toric_degeneration_table_tex())
    else
        path = split(only(outputs), '='; limit=2)[2]
        write_toric_degeneration_table_tex(path)
        println("Wrote ", path)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    toric_degeneration_table_main()
end
