const PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))
const REPOSITORY_ROOT = normpath(joinpath(PACKAGE_ROOT, ".."))

"Return a path relative to the Julia project, independently of the caller's cwd."
project_path(parts...) = normpath(joinpath(PACKAGE_ROOT, parts...))

"Return a path below the single canonical data directory."
data_path(parts...) = project_path("data", parts...)

"Return a path below a named study directory."
study_path(study::Union{Symbol, AbstractString}, parts...) =
    project_path("studies", lowercase(string(study)), parts...)

"Return a path below the ignored resumable-work directory."
work_path(parts...) = project_path("work", parts...)

"Shared TeX style used by every local-cell and tree atlas."
atlas_style_path() = project_path("assets", "tree_atlas_style.tex")
