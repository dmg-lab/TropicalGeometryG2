# The second coefficient field of the local G2 flag cell atlas. Rays
# R_1, ..., R_51 are those of the QQ census; the rays that occur only over
# GF(2) continue the numbering. The support cells of both fields are written
# in these atlas rays, so that a cell present in both fields has one name.

using TOML

const GF2_GRASS_SYMMETRY = TropicalCensus.data_path("g2", "grassmannian_gf2", "05_symmetry")
const GF2_FLAG_FAN = TropicalCensus.data_path("g2", "flag_gf2", "02_coarsest_fan.mrdi")
const GF2_FLAG_COARSEN_PROVENANCE =
    TropicalCensus.data_path("g2", "flag_gf2", "provenance", "coarsen.toml")
const QQ_FLAG_FVECTOR = [1, 141, 1352, 3584, 2714]

field_tex(field::Symbol) = field == :qq ? "\\mathbb{Q}" : "\\mathbb{F}_2"

"""
The Grassmannian support of the GF(2) census in atlas rays: the ray matrix
extending `balanced_rays`, the atlas index of every census ray, the ray
orbits extending `qq_ray_orbits`, and one spec per census cell orbit. A spec
of an orbit present over QQ takes the representative and member order of its
QQ spec, recorded as `qq_orbit`; both fields then draw the same cell.
"""
function gf2_support_data(balanced_rays, qq_specs, qq_ray_orbits)
    census = census_atlas_rays(:gf2, balanced_rays)
    rays, ray_index = census.rays, census.ray_index
    first_new_ray = nrows(balanced_rays) + 1

    symmetry = TropicalCensus.load_symmetry_artifacts(GF2_GRASS_SYMMETRY)
    known = Dict(Set(orbit) => number for (number, orbit) in enumerate(qq_ray_orbits))
    ray_orbits = [collect(orbit) for orbit in qq_ray_orbits]
    for orbit in symmetry.ray_orbits
        atlas_orbit = sort!(ray_index[orbit])
        haskey(known, Set(atlas_orbit)) && continue
        all(>=(first_new_ray), atlas_orbit) ||
            error("a GF(2) ray orbit mixes common and GF(2)-only rays")
        push!(ray_orbits, atlas_orbit)
    end
    sort!(reduce(vcat, ray_orbits)) == collect(1:nrows(rays)) ||
        error("the ray orbits of both fields do not partition the atlas rays")

    cells = census_atlas_cells(:gf2, ray_index)
    qq_orbit_of = Dict{Tuple{Vararg{Int}},Int}()
    for (orbit_number, spec) in enumerate(qq_specs), member in spec.members
        qq_orbit_of[Tuple(sort!(collect(member)))] = orbit_number
    end
    specs = Any[]
    for (census_orbit, orbit) in enumerate(symmetry.cell_orbits)
        members = cells[orbit]
        qq_orbits = unique(get(qq_orbit_of, member, 0) for member in members)
        length(qq_orbits) == 1 || error("GF(2) support orbit $census_orbit is only partly present over QQ")
        qq_orbit = only(qq_orbits)
        if qq_orbit > 0
            qq_spec = qq_specs[qq_orbit]
            Set(Tuple(sort!(collect(member))) for member in qq_spec.members) == Set(members) ||
                error("GF(2) support orbit $census_orbit is not the whole QQ orbit $qq_orbit")
            representative, members = qq_spec.representative, qq_spec.members
        else
            representative = first(members)
            members = Tuple(members)
        end
        kind = length(representative) == 3 ? :triangle :
               length(representative) == 4 ? :square :
               error("GF(2) support cell $representative is neither a triangle nor a square")
        push!(specs, (
            kind=kind,
            representative=representative,
            members=members,
            census_orbit=census_orbit,
            qq_orbit=qq_orbit > 0 ? qq_orbit : nothing,
        ))
    end
    return (rays=rays, ray_index=ray_index, ray_orbits=ray_orbits, specs=specs)
end

"The f-vector modulo lineality of the coarsest GF(2) flag fan, from its census."
function gf2_flag_f_vector()
    metadata = TOML.parsefile(GF2_FLAG_COARSEN_PROVENANCE)["metadata"]
    return Int.(metadata["f_vector_modulo_lineality_including_zero"])
end
