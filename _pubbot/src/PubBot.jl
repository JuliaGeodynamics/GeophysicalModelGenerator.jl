"""
Keeps the Publications page of the documentation (`docs/src/man/publications.md`) up to date
with papers that use GeophysicalModelGenerator.jl.

- `weekly`: scans OpenAlex for recent works that name or cite the package, or are by its
  developers, and edits the Publications page.
- `add_dois`: formats and inserts publications from a list of DOIs.
"""
module PubBot

using Dates, Downloads, JSON, TOML, Unicode
using Poppler_jll: pdftotext

const BOTDIR = dirname(@__DIR__)
const PAGE = joinpath(dirname(BOTDIR), "docs", "src", "man", "publications.md")
const SEEN = joinpath(BOTDIR, "seen_dois.txt")

include("util.jl")
include("openalex.jl")
include("fulltext.jl")
include("format.jl")
include("page.jl")
include("main.jl")

end
