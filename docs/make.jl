using Documenter
using DocumenterCitations
using MonteCarloScattering

bib = CitationBibliography(joinpath(@__DIR__, "refs.bib"))

makedocs(;
    sitename = "MonteCarloScattering",
    plugins = [bib],
)
