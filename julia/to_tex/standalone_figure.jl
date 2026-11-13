###############################################################################
# A figure of main.tex as a standalone document: the preamble of main.tex in
# the class standalone, and the figure without the spaces that the line ends
# of its file leave around it. The PDF of the document is the figure at its
# size in main.tex with a border, which \figurepdf in main.tex trims.
#
#   julia to_tex/standalone_figure.jl FIGURE.tex [--output=FILE]
#
# writes the document of a hand-drawn figure; the to_tex scripts write theirs
# with --standalone, through g2_standalone_figure. A document compiles from
# the root of the repository, where the figures that it inputs are.
###############################################################################

const G2_STANDALONE_MAIN_TEX = normpath(joinpath(@__DIR__, "..", "..", "main.tex"))
# The border of a figure, the trim of \figurepdf in main.tex.
const G2_STANDALONE_BORDER_BP = 72

# The lines around the fragment in a document, which g2_standalone_fragment
# finds; what follows the second may change.
const G2_STANDALONE_BEFORE = "\\setbox\\gtwostandalonefigure\\hbox{%\n"
const G2_STANDALONE_END = "}% The end of the figure.\n"
const G2_STANDALONE_AFTER = G2_STANDALONE_END * raw"""
% The one box of the figure, without the spaces around it: e-TeX's
% \lastnodetype is 11 for glue and -1 for an empty list.
\setbox\gtwostandalonefigure\hbox{\unhbox\gtwostandalonefigure
  \loop\ifnum\lastnodetype=11 \unskip\repeat
  \global\setbox1\lastbox
  \loop\ifnum\lastnodetype=11 \unskip\repeat
  \xdef\gtwostandaloneleft{\the\lastnodetype}}%
\ifvoid1 \errmessage{The figure does not end with a box}\fi
\ifnum\gtwostandaloneleft=-1 \else \errmessage{The figure is more than one box between spaces}\fi
\box1
\end{document}
"""

"""
    g2_standalone_figure(fragment; main_tex=G2_STANDALONE_MAIN_TEX)

The standalone document of the figure `fragment`: the preamble of `main_tex`
in the class `standalone`, and the fragment in a box. The fragment draws one
box between spaces, which a paragraph of main.tex drops; the page is that box
with a border of $(G2_STANDALONE_BORDER_BP)bp, which \\figurepdf in main.tex trims
but does not clip, so that ink up to the border beyond the box, such as half a
line on its edge, shows as in the picture. Anything else beside the box is an
error.
"""
function g2_standalone_figure(fragment::AbstractString;
                              main_tex::AbstractString=G2_STANDALONE_MAIN_TEX)
    pieces = split(read(main_tex, String), "\\begin{document}"; limit=2)
    length(pieces) == 2 || error("cannot locate the preamble of $main_tex")
    preamble = pieces[1]
    class = match(r"^\\documentclass(?:\[([^\]]*)\])?\{([^}]*)\}"m, preamble)
    class === nothing && error("cannot locate the document class of $main_tex")
    options = class[1] === nothing || isempty(class[1]) ? "" : class[1] * ","
    preamble = preamble[1:prevind(preamble, class.offset)] *
               "\\documentclass[$(options)class=$(class[2]),border=$(G2_STANDALONE_BORDER_BP)bp]{standalone}" *
               preamble[class.offset+ncodeunits(class.match):end]
    return "% A figure of main.tex at its size there, " *
           "written by julia/to_tex/standalone_figure.jl.\n" * preamble * raw"""
% The page geometry of the class standalone, which the preamble of main.tex changes.
\topmargin=0pt \headheight=0pt \headsep=0pt \footskip=0pt
\oddsidemargin=0pt \evensidemargin=0pt \topskip=0pt
\newbox\gtwostandalonefigure
\begin{document}
""" * G2_STANDALONE_BEFORE * fragment * (endswith(fragment, '\n') ? "" : "\n") * G2_STANDALONE_AFTER
end

"""
    g2_standalone_fragment(document)

The fragment of a document of `g2_standalone_figure`, which that function
wraps again in the preamble of main.tex as it is now; `nothing` for another
document.
"""
function g2_standalone_fragment(document::AbstractString)
    before = findfirst(G2_STANDALONE_BEFORE, document)
    after = findlast("\n" * G2_STANDALONE_END, document)
    (before === nothing || after === nothing || last(before) >= first(after)) && return nothing
    return document[nextind(document, last(before)):first(after)]
end

function g2_standalone_figure_main(args=ARGS)
    usage = "Usage: julia to_tex/standalone_figure.jl FIGURE.tex [--output=FILE]"
    if "--help" in args
        println(usage)
        return nothing
    end
    figures = filter(arg -> !startswith(arg, "--"), args)
    outputs = [arg[length("--output=")+1:end] for arg in args if startswith(arg, "--output=")]
    length(figures) == 1 && length(outputs) <= 1 &&
        all(arg -> !startswith(arg, "--") || startswith(arg, "--output="), args) ||
        error(usage)
    tex = g2_standalone_figure(read(only(figures), String))
    isempty(outputs) ? print(tex) : write(only(outputs), tex)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    g2_standalone_figure_main()
end
