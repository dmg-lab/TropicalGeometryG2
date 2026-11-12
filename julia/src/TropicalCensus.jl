module TropicalCensus

using Dates
using Distributed
using LinearAlgebra
using Oscar
using Printf
using Random
using SHA
using TOML

include("paths.jl")
include("progress.jl")
include("trees.jl")
include("symmetry.jl")
include("coarsening.jl")
include("orbit_intersection.jl")
using .OrbitIntersection: expand_restricted_prevariety, canonical_all_ones_ray
include("parallel_tropicalization.jl")
include("degenerations.jl")
include("artifacts.jl")
include("projections.jl")
# The Chevalley polytopes of Spacek--Wang live with their study.
include(joinpath("..", "studies", "g2", "degenerations_grassmannian", "Spacek_Wang",
                 "chevalley_polytopes.jl"))
include("census.jl")
include("recursive_polytope_types.jl")

using .ChevalleyPolytopes
using .RecursivePolytopeTypes

include("summary.jl")

export TropicalProblem,
       TropicalizationConfig,
       run_tropicalization,
       parse_gfan_fan,
       gfan_prevariety,
       expand_restricted_prevariety,
       canonical_all_ones_ray,
       CoarseningStrategy,
       SupportLinkCoarsening,
       EqualSpanCoarsening,
       coarsen_fan,
       coarsen_prevariety_fan,
       prevariety_dimension_skeleton,
       coarsen_fan_file,
       open_equal_span_walls,
       open_wall_certificate,
       TreeMetric,
       TreeGauge,
       normalize_ray,
       normalize_fan,
       trees_from_fan,
       tree_from_min_pluecker,
       tree_distance,
       tree_edge_lengths,
       pendant_lengths,
       internal_splits,
       star_pluecker,
       primitive_integer,
       canonical_cell_weights,
       normalized_cell_trees,
       taxon_symmetry,
       all_permutations,
       action_orbits,
       coordinate_permutation,
       tropical_coordinate_permutation,
       symmetry_coordinate_permutation,
       signed_permutation_group,
       signed_coordinate_images,
       monomial_support_automorphisms,
       quotient_ray_signature,
       preserves_lineality,
       DegenerationConfig,
       canonical_initial_ideal,
       classify_degenerations,
       affine_unimodular_equivalent,
       verify_canonical_weight_action,
       degeneration_artifact_paths,
       certify_ray_projection,
       verify_ray_projection_certificate,
       write_ray_projection_certificate,
       project_path,
       data_path,
       study_path,
       work_path,
       atlas_style_path,
       CensusStage,
       tropicalize,
       coarsen,
       normalize,
       trees,
       symmetry,
       degenerations,
       atlas,
       projections,
       summarize,
       zenodo,
       CensusOptions,
       StudyTarget,
       StudyCensus,
       parse_census_options,
       census_configuration,
       run_census,
       write_provenance,
       stage_provenance_path,
       coarsest_report_path,
       summary_path,
       census_summary,
       write_census_summary,
       ChevalleyData,
       chevalley_data,
       chevalley_polytope,
       parabolic_indices,
       minimal_coset_representative,
       paper_word,
       all_paper_words,
       heap_relations,
       degree_of_embedding,
       normalised_volume_of,
       has_idp,
       verify,
       is_adapted_basis,
       reduced_subwords,
       decomposition_polytope,
       check_decomposition,
       save_polymake,
       save_chevalley_polytope,
       polytope_name,
       weight_multiplicities,
       determine_polytope,
       describe_type

end
