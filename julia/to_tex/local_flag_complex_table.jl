###############################################################################
# The paper table of the local flag complexes of the tropical G2 flag variety,
# in characteristic 0 and in characteristic 2.
#
# Above every maximal cell S of the tropical G2-Grassmannian stands the
# subcomplex of the tropical G2 flag variety made of the maximal flag cells
# that project onto S, with all their faces. Its support is the preimage of S.
# The script computes this complex for one cell per orbit of maximal cells,
# in every field in which the orbit occurs, with its f-vector modulo lineality
# and its combinatorial type, and groups the orbits by isomorphism type.
#
# Run from the repository root with julia --project=julia. The script reads
# the normalized coarsest flag fan, the projection census and the stage-6
# Grassmannian census of each field, and takes the orbit numbers from the
# G2-tree atlas. It prints the tabularx, or writes it to --output=FILE;
# main.tex supplies the table float and its caption, and loads the packages
# tabularx and booktabs.
###############################################################################

using Oscar
using TropicalCensus

include(joinpath(@__DIR__, "toric_degeneration_table.jl"))

const LOCAL_FLAG_DIRECTORIES = (
    qq=(flag="flag", grassmannian="grassmannian"),
    gf2=(flag="flag_gf2", grassmannian="grassmannian_gf2"),
)

# The space between two types in the list of the types of an f-vector, where
# the list may break.
const LOCAL_FLAG_TYPE_SEP = "\\quad"

# The \arraystretch of the rows, in the booktabs form of the tables of main.tex.
const LOCAL_FLAG_TABLE_STRETCH = 0.9

# The words of the maximal-cell types of the census, singular and plural.
const LOCAL_FLAG_CELL_WORDS = Dict(
    "3-simplex" => ("tetrahedron", "tetrahedra"),
    "Pyramid(4-gon)" => ("square pyramid", "square pyramids"),
    "Δ1 × Δ2" => ("triangular prism", "triangular prisms"),
)

local_flag_data_file(field, target, name) =
    TropicalCensus.data_path("g2", LOCAL_FLAG_DIRECTORIES[field][target], name)

"""
The flag fan of `field` modulo lineality: its rays, its lineality, its
maximal cells, and the flag cells over every maximal cell of the
Grassmannian, from the certified projection census `07_projection_cells.csv`.
The cells of both fans are indexed as `fan_cells` of the normalized coarsest
fans, as in the projection stage.
"""
function local_flag_field_data(field::Symbol)
    fan = load(local_flag_data_file(field, :flag, "03_normalized_coarsest_fan.mrdi"))
    data = rays_modulo_lineality(fan)
    cells = TropicalCensus.fan_cells(fan)
    over = Dict{Int,Vector{Int}}()
    for row in TropicalCensus.csv_rows(local_flag_data_file(field, :flag, "07_projection_cells.csv"))
        push!(get!(over, parse(Int, row["grassmannian_cell"]), Int[]), parse(Int, row["flag_cell"]))
    end
    sum(length, values(over)) == length(cells) || error("the projection census misses flag cells over $field")
    census = TropicalCensus.csv_rows(local_flag_data_file(field, :grassmannian, "06_canonical_interior_degenerations.csv"))
    orbits = Dict(parse(Int, row["orbit"]) => parse.(Int, split(row["cells"], ';')) for row in census)
    symmetry = TropicalCensus.load_symmetry_artifacts(local_flag_data_file(field, :grassmannian, "05_symmetry"))
    issetequal(Set.(values(orbits)), Set.(symmetry.cell_orbits)) ||
        error("the census and the symmetry orbits of the Grassmannian over $field disagree")
    return (
        field, cells, over, orbits,
        rays=[QQ.(collect(ray)) for ray in data.rays_modulo_lineality],
        lineality=matrix(QQ, data.lineality_basis),
    )
end

"""
The local flag complex over the Grassmannian maximal cell `cell`: its faces
by dimension, each a sorted list of flag rays, its f-vector, the types of its
maximal cells, and its vertex-face incidence. The faces of every flag cell
come from the facets of its own cone, so no other cell of the fan is used.
"""
function local_flag_complex(data, cell::Int)
    flag_cells = data.over[cell]
    faces = [Set{Vector{Int}}() for _ in 0:3]
    types = Dict{String,Int}()
    for flag_cell in flag_cells
        levels, incidence = TropicalCensus.cell_faces(data.cells[flag_cell], data.rays, data.lineality)
        length(levels) == 4 || error("flag cell $flag_cell over $(data.field) is not 3-dimensional")
        for (level, list) in enumerate(levels)
            union!(faces[level], list)
        end
        type = describe_type(determine_polytope(incidence, 3))
        types[type] = get(types, type, 0) + 1
    end
    faces = [sort!(collect(level)) for level in faces]
    vertices = [only(face) for face in faces[1]]
    position = Dict(ray => index for (index, ray) in enumerate(vertices))
    # A face is the cone over its rays, so the vertex-face incidence of the
    # faces of positive dimension determines the face lattice.
    columns = reduce(vcat, faces[2:end])
    incidence = IncidenceMatrix([[position[ray] for ray in face] for face in columns])
    return (; cell, flag_cells, faces, f_vector=length.(faces), types, incidence)
end

"""Whether two local flag complexes are combinatorially isomorphic."""
local_flag_isomorphic(a, b) =
    a.f_vector == b.f_vector && a.types == b.types &&
    Polymake.graph.isomorphic(a.incidence, b.incidence)::Bool

"""
The maximal-cell orbits of both fields in atlas order: the `number` of the
tables of trees and the census `orbits` of the fields in which they occur.
"""
function local_flag_orbit_rows()
    rows = []
    for row in toric_table_csv_rows(TORIC_TABLE_ATLAS)
        fields = split(row["fields"], '+')
        all(field -> field in ("qq", "gf2"), fields) || continue
        number = toric_table_label_number(row["label"])
        any(other -> other.number == number, rows) && continue
        orbits = [field => parse(Int, row["$(field)_orbit"]) for field in TORIC_TABLE_FIELDS
                  if string(field) in fields]
        push!(rows, (; number, orbits))
    end
    sort!(rows; by=row -> row.number)
    [row.number for row in rows] == collect(1:length(rows)) ||
        error("the atlas does not number the maximal-cell orbits 1, 2, ...")
    return rows
end

"""
The isomorphism types of the local flag complexes. Each type has the complex
of its first orbit and its `members`: the complexes of all orbits of this
type, each with its orbit `number` of the tables of trees, its `field`, its
census `orbit` and its `orbit_size`. An orbit of both fields is computed in
both, over one representative cell.
"""
function local_flag_types(data=NamedTuple{TORIC_TABLE_FIELDS}(map(local_flag_field_data, TORIC_TABLE_FIELDS)),
                          rows=local_flag_orbit_rows())
    for field in TORIC_TABLE_FIELDS
        listed = sort!([orbit for row in rows for (other, orbit) in row.orbits if other == field])
        listed == sort!(collect(keys(data[field].orbits))) ||
            error("$TORIC_TABLE_ATLAS does not list the orbits of the $(TORIC_TABLE_FIELD_NAMES[field]) census once each")
    end
    types = []
    for row in rows, (field, orbit) in row.orbits
        members = data[field].orbits[orbit]
        complex = local_flag_complex(data[field], first(members))
        complex = merge(complex, (; number=row.number, field, orbit, orbit_size=length(members)))
        index = findfirst(type -> local_flag_isomorphic(type.complex, complex), types)
        if index === nothing
            push!(types, (; complex, members=[complex]))
        else
            push!(types[index].members, complex)
        end
    end
    return types
end

"""
The types in table order, by f-vector and then by their least orbit number,
each with its `letter`, a fraktur capital for each f-vector, and its `label`:
the letter and, when several types share the f-vector, an index.
"""
function local_flag_labelled_types(types=local_flag_types())
    order = sort(types; by=type -> (type.complex.f_vector, minimum(m.number for m in type.members)))
    vectors = unique(type.complex.f_vector for type in order)
    length(vectors) <= 26 || error("more f-vectors than letters")
    return map(order) do type
        group = findfirst(==(type.complex.f_vector), vectors)
        letter = "\\mathfrak " * string(Char(Int('A') + group - 1))
        same = [other for other in order if other.complex.f_vector == type.complex.f_vector]
        label = length(same) == 1 ? letter : "$(letter)_{$(findfirst(other -> other === type, same))}"
        merge(type, (; letter, label))
    end
end

"""
The maximal cells of a complex in words, in the order of the census types,
one line for each type.
"""
function local_flag_cells_tex(counts)
    issubset(keys(counts), keys(LOCAL_FLAG_CELL_WORDS)) || error("unknown maximal-cell type in $(keys(counts))")
    lines = ["$(counts[type]) $(LOCAL_FLAG_CELL_WORDS[type][counts[type] == 1 ? 1 : 2])"
             for type in ("3-simplex", "Pyramid(4-gon)", "Δ1 × Δ2") if haskey(counts, type)]
    length(lines) == 1 && return only(lines)
    return "\\makecell[tl]{" * join(lines, ",\\\\ ") * "}"
end

"""
The orbit numbers of a type. An orbit whose complex has this type in only
some of the fields in which it occurs carries their characteristics.
"""
function local_flag_orbits_tex(type, rows)
    occurs = Dict(row.number => [field for (field, _) in row.orbits] for row in rows)
    numbers = sort!(unique(m.number for m in type.members))
    entries = map(numbers) do number
        fields = [m.field for m in type.members if m.number == number]
        entry = "\\textcolor{g2CellLabel}{\\bfseries\\boldmath\\oldstylenums{$number}}"
        issetequal(fields, occurs[number]) && return entry
        entry * "\\,(char.~" * join((TORIC_TABLE_CHARACTERISTIC[field] for field in fields), ",") * ")"
    end
    return join(entries, ", ")
end

"""
The types `group`, which share an f-vector: a list of each label with its
orbits, which breaks only between two types.
"""
local_flag_subtypes_tex(group, rows) =
    join(("\\mbox{\$$(type.label)\$: $(local_flag_orbits_tex(type, rows))}" for type in group),
         LOCAL_FLAG_TYPE_SEP * " ")

"""
    local_flag_complex_table_tex()

The input-ready tabularx of the text width: one row per f-vector of local
flag complexes, with the f-vector, the maximal cells of the complexes and the
isomorphism types with this f-vector, each with the orbits of maximal cells of
the Grassmannian, numbered as in the tables of trees, above which it stands.
"""
function local_flag_complex_table_tex()
    data = NamedTuple{TORIC_TABLE_FIELDS}(map(local_flag_field_data, TORIC_TABLE_FIELDS))
    rows = local_flag_orbit_rows()
    types = local_flag_labelled_types(local_flag_types(data, rows))
    for field in TORIC_TABLE_FIELDS
        # Every flag cell lies over exactly one Grassmannian cell, so the orbits
        # account for all flag cells.
        total = sum(m.orbit_size * length(m.flag_cells) for type in types for m in type.members if m.field == field)
        total == length(data[field].cells) ||
            error("the local flag complexes over $(TORIC_TABLE_FIELD_NAMES[field]) have $total maximal cells, not $(length(data[field].cells))")
    end
    io = IOBuffer()
    println(io, "% Generated by julia/to_tex/local_flag_complex_table.jl from")
    for field in TORIC_TABLE_FIELDS, (target, name) in ((:flag, "03_normalized_coarsest_fan.mrdi"),
            (:flag, "07_projection_cells.csv"), (:grassmannian, "06_canonical_interior_degenerations.csv"))
        println(io, "%   ", relpath(local_flag_data_file(field, target, name), TORIC_TABLE_ROOT))
    end
    println(io, "%   $TORIC_TABLE_ATLAS")
    println(io, "% $(length(types)) isomorphism types of local flag complexes over the $(length(rows)) orbits ",
            "of maximal cells; f-vectors modulo lineality.")
    println(io, "% Each type: the census orbit, the representative Grassmannian cell and the number of flag cells over it.")
    for type in types
        println(io, "%   ", type.label, ": ", join(("$(m.number) = $(TORIC_TABLE_FIELD_NAMES[m.field]) orbit $(m.orbit) " *
                "(cell $(m.cell), $(length(m.flag_cells)) flag cells)" for m in type.members), "; "))
    end
    println(io, "\\begingroup")
    println(io, "\\definecolor{g2CellLabel}{HTML}{$TORIC_TABLE_LABEL_COLOR}")
    println(io, "\\renewcommand{\\arraystretch}{$LOCAL_FLAG_TABLE_STRETCH}")
    println(io, "\\begin{tabularx}{\\linewidth}{@{}l l >{\\raggedright\\arraybackslash}X@{}}")
    println(io, "  \\toprule")
    println(io, "  \$f\$-vector & Maximal cells & Subtypes and orbits \\tabularnewline")
    # The only rule inside the table: under the heading.
    println(io, "  \\midrule")
    for letter in unique(type.letter for type in types)
        group = [type for type in types if type.letter == letter]
        # In the data the f-vector fixes the maximal cells; the first type
        # stands for its group.
        all(type -> type.complex.types == first(group).complex.types, group) ||
            error("the types $letter share an f-vector but not their maximal cells")
        println(io, "  ", join((
            toric_table_f_vector_tex(first(group).complex.f_vector),
            local_flag_cells_tex(first(group).complex.types),
            local_flag_subtypes_tex(group, rows),
        ), " & "), " \\tabularnewline")
    end
    println(io, "  \\bottomrule")
    println(io, "\\end{tabularx}%")
    println(io, "\\endgroup")
    return String(take!(io))
end

function local_flag_complex_table_main(args=ARGS)
    if "--help" in args
        println("Usage: julia --project=julia julia/to_tex/local_flag_complex_table.jl [--output=FILE]")
        return nothing
    end
    for arg in args
        startswith(arg, "--output=") || throw(ArgumentError("unknown option $arg; use --help"))
    end
    outputs = [arg for arg in args if startswith(arg, "--output=")]
    length(outputs) <= 1 || throw(ArgumentError("pass only one --output"))
    tex = local_flag_complex_table_tex()
    if isempty(outputs)
        print(tex)
    else
        path = split(only(outputs), '='; limit=2)[2]
        mkpath(dirname(abspath(path)))
        write(path, tex)
        println("Wrote ", path)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    local_flag_complex_table_main()
end
