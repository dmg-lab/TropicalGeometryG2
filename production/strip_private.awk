###############################################################################
# Print a file without its private blocks: the lines from a line that opens a
# private block to the line that closes it, both included. The markers depend
# on the extension of the file, each on a line of its own, blanks around it
# allowed:
#
#   <!-- private -->  <!-- /private -->   .md, .html
#   # private         # /private          .jl, .py, .sh, .slurm, .toml
#   % private         % /private          .tex
#
# Blocks do not nest. Exit status: 0 when the file mentions its opening marker
# (its public version is printed); 1 when it does not, or has another
# extension (the file is public as it is; what is printed is to be ignored);
# 2 for a nested, unopened or unclosed block (a message on stderr). Used by
# the release scripts, make_zenodo_tarball.sh among them. From julia/:
#
#   awk -f ../production/strip_private.awk FILE > PUBLIC
###############################################################################

BEGIN {
    file = ARGV[1]
    if (file ~ /\.(md|html)$/)                  { opener = "<!-- private -->"; closer = "<!-- /private -->" }
    else if (file ~ /\.(jl|py|sh|slurm|toml)$/) { opener = "# private";        closer = "# /private" }
    else if (file ~ /\.tex$/)                   { opener = "% private";        closer = "% /private" }
    else { status = 1; exit 1 }
}

function fail(message) {
    print file ": " message > "/dev/stderr"
    status = 2
    exit 2
}

index($0, opener) { mentioned = 1 }
{ line = $0; sub(/^[ \t]+/, "", line); sub(/[ \t\r]+$/, "", line) }
line == opener { if (inside) fail("nested private block at line " NR); inside = 1; next }
line == closer { if (!inside) fail("unopened private block at line " NR); inside = 0; next }
!inside { print }

END {
    if (status) exit status
    if (inside) fail("unclosed private block")
    exit (mentioned ? 0 : 1)
}
