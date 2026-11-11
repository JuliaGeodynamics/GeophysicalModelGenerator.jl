#=
Full text comes from open-access PDFs. Copernicus and arXiv serve PDFs to scripts; Wiley/AGU,
OUP and Elsevier do not, so for those OpenAlex's own full-text index is the only evidence.
=#

const CACHE = joinpath(BOTDIR, ".cache")

"Copernicus PDFs live next to the article landing page, named after the DOI suffix."
function copernicus_pdf(doi::AbstractString)
    status, _, landing = http_get(doi_url(doi); method = "HEAD")
    status == 200 || return nothing
    return rstrip(landing, '/') * "/" * split(doi, '/'; limit = 2)[2] * ".pdf"
end

function pdf_to_text(bytes::Vector{UInt8})
    return mktemp() do path, io
        write(io, bytes)
        close(io)
        try
            read(`$(pdftotext()) -q -enc UTF-8 $path -`, String)
        catch
            ""
        end
    end
end

# Publishers that refuse PDF downloads from scripts (403), not worth a request.
const BLOCKED_HOSTS = r"^https?://([^/]+\.)?(wiley\.com|oup\.com|sciencedirect\.com|tandfonline\.com)/"

"Text of the paper's open-access PDF, or `nothing`. Cached in `_pubbot/.cache`."
function fulltext(w::Work)
    path = joinpath(CACHE, replace(w.doi, r"[^A-Za-z0-9._-]" => "_") * ".txt")
    if isfile(path)
        text = read(path, String)
        return isempty(text) ? nothing : text
    end
    text = ""
    for url in (startswith(w.doi, "10.5194/") ? [nothing; w.pdf_urls] : w.pdf_urls)
        url === nothing && (url = copernicus_pdf(w.doi); url === nothing && continue)
        occursin(BLOCKED_HOSTS, url) && continue
        status, body, _ = http_get(url; timeout = 45, retries = 1)
        status == 200 && length(body) > 4 && body[1:4] == b"%PDF" || continue
        # Undo line-break hyphenation so the package name matches.
        text = replace(pdf_to_text(body), r"(\w)-\s*\n\s*(\w)" => s"\1\2", r"\s+" => " ")
        isempty(text) || break
    end
    isempty(text) && return nothing
    mkpath(CACHE)
    write(path, text)
    return text
end

Base.@kwdef mutable struct Assessment
    work::Work
    evidence::Vector{String} = String[]
    pdf::Bool = false                 # full text checked
    bucket::Symbol = :none            # :strong (✅), :review (🟡) or :none (❌)
    entry::Union{String, Nothing} = nothing
    year::Union{Int, Nothing} = nothing
    note::String = ""
end

"""
    assess(w, cfg, refs, text, mentions)

✅ `:strong` when `w` names the package (title, abstract or PDF `text`) or is a package paper;
🟡 `:review` when only the OpenAlex full-text index (`mentions`) finds the name or it only cites a
package paper (`refs`); ❌ `:none` otherwise. The index stems words, so its hits need
confirmation.
"""
function assess(w::Work, cfg::Config, refs::Set{String}, text::Union{String, Nothing}, mentions::Set{String})
    a = Assessment(; work = w, pdf = text !== nothing)
    w.doi in cfg.dois && push!(a.evidence, "package paper")
    occursin(cfg.pattern, w.title * " " * w.abstract) && push!(a.evidence, "title/abstract")
    text !== nothing && occursin(cfg.pattern, text) && push!(a.evidence, "full text")
    named = !isempty(a.evidence)
    w.id in mentions && push!(a.evidence, "OpenAlex full-text index")
    isdisjoint(w.refs, refs) || push!(a.evidence, "cites package paper")
    a.bucket = named ? :strong : isempty(a.evidence) ? :none : :review
    return a
end

"""
Demote proposed preprints whose published version is already listed, evaluated before or
proposed in this scan; titles can differ too much for `select_candidates` to merge them.
"""
function drop_published_preprints!(
        assessments::Vector{Assessment}, cfg::Config;
        on_page::Set{String}, seen::Set{String}
    )
    known = on_page ∪ seen ∪ Set(a.work.doi for a in assessments if a.bucket != :none)
    for a in assessments
        a.bucket != :none && is_preprint(a.work, cfg) || continue
        published = crossref_published_version(a.work.doi)
        (published === nothing || !(published in known)) && continue
        a.bucket = :none
        a.note = "preprint of $published"
    end
    return assessments
end
