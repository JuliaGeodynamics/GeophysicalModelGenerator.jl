const OPENALEX = "https://api.openalex.org"
const WORK_FIELDS = "id,doi,title,publication_date,type,authorships,referenced_works,abstract_inverted_index,locations"

Base.@kwdef mutable struct Work
    id::String                                # OpenAlex ID, e.g. "W4417170000"
    doi::String                               # normalized
    title::String
    date::Date
    type::String
    authors::Vector{String}                   # display names
    refs::Set{String} = Set{String}()         # OpenAlex IDs of referenced works
    abstract::String = ""
    pdf_urls::Vector{String} = String[]
    other_dois::Vector{String} = String[]     # other versions merged into this work
end

short_id(url::AbstractString) = String(last(split(url, '/')))

function openalex(path::AbstractString, params)
    params = collect(Pair{String, String}, params)
    key = get(ENV, "OPENALEX_API_KEY", "")
    isempty(key) || push!(params, "api_key" => key)
    return get_json("$OPENALEX/$path?" * query(params))
end

"""
All works matching an OpenAlex `filter`, following cursor pagination. Pages of 50: larger
pages with abstracts take long enough to stall.
"""
function openalex_works(filter::AbstractString; select::AbstractString = WORK_FIELDS)
    results = []
    cursor = "*"
    while cursor isa AbstractString
        r = openalex("works", ["filter" => filter, "select" => select, "per-page" => "50", "cursor" => cursor])
        isempty(r["results"]) && break
        append!(results, r["results"])
        cursor = get(r["meta"], "next_cursor", nothing)
    end
    return results
end

function abstract_text(inverted)
    inverted isa AbstractDict || return ""
    positions = Tuple{Int, String}[]
    for (word, idxs) in inverted, i in idxs
        push!(positions, (i, word))
    end
    return join(last.(sort!(positions)), " ")
end

function parse_work(x)
    pdfs = String[]
    for loc in something(get(x, "locations", nothing), [])
        url = get(loc, "pdf_url", nothing)
        url isa AbstractString && push!(pdfs, url)
    end
    return Work(;
        id = short_id(x["id"]), doi = normalize_doi(something(x["doi"], "")),
        title = clean_text(something(x["title"], "")), date = Date(x["publication_date"]), type = x["type"],
        authors = [string(something(get(a["author"], "display_name", nothing), "")) for a in x["authorships"]],
        refs = Set(short_id.(something(get(x, "referenced_works", nothing), []))),
        abstract = abstract_text(get(x, "abstract_inverted_index", nothing)), pdf_urls = unique!(pdfs)
    )
end

"OpenAlex IDs of the package papers."
function package_ids(cfg::Config)
    isempty(cfg.dois) && return Set{String}()
    ids = Set(short_id(x["id"]) for x in openalex_works("doi:" * join(cfg.dois, "|"); select = "id"))
    isempty(ids) && error("OpenAlex knows none of the package papers $(cfg.dois)")
    return ids
end

"""
Recent works by the developers, citing a package paper (`refs`), or whose full text (OpenAlex
index) names the package. Returns the works and the IDs of those the index finds the name in.
"""
function fetch_works(cfg::Config, refs::Set{String}; since::Date)
    recent = "from_publication_date:$since,type:$(join(cfg.types, "|"))"
    works, ids = Work[], Set{String}()
    function add!(f)
        found = parse_work.(openalex_works("$f,$recent"))
        for w in found
            w.id in ids || (push!(ids, w.id); push!(works, w))
        end
        return found
    end
    add!("authorships.author.orcid:" * join(cfg.developers, "|"))
    isempty(refs) || add!("cites:" * join(refs, "|"))
    mentions = Set(w.id for w in add!("fulltext.search:\"$(cfg.match)\""))
    return works, mentions
end

family_names(w::Work) = Set(last(words(name)) for name in w.authors if !isempty(words(name)))

"Two records of the same paper, e.g. a preprint and its published version."
same_paper(a::Work, b::Work) =
    similar_titles(a.title, b.title) && !isempty(family_names(a) ∩ family_names(b))

is_preprint(w::Work, cfg::Config) = w.type == "preprint" || preprint_server(w.doi, cfg) !== nothing

"""
Drop repository copies and non-papers, merge versions of the same paper (the published
version wins), and skip works already listed or evaluated before.
"""
function select_candidates(works::Vector{Work}, cfg::Config; on_page::Set{String}, seen::Set{String})
    keep = filter(works) do w
        !isempty(w.doi) && !any(startswith(w.doi, p) for p in cfg.ignore_doi_prefixes) &&
            !any(occursin(r, w.title) for r in cfg.ignore_titles)
    end
    unique!(w -> w.doi, keep)
    sort!(keep; by = w -> (is_preprint(w, cfg), w.date))
    merged = Work[]
    for w in keep
        k = findfirst(m -> same_paper(m, w), merged)
        k === nothing ? push!(merged, w) : push!(merged[k].other_dois, w.doi)
    end
    # Versions of listed preprints are handled by `published_versions`, not added again.
    return filter(merged) do w
        !(w.doi in on_page || w.doi in seen || any(d in on_page for d in w.other_dois))
    end
end

"DOI of the published version Crossref records for a preprint, or `nothing`."
function crossref_published_version(doi::AbstractString)
    params = ["filter" => "relation.type:has-preprint,relation.object:$doi", "select" => "DOI", "rows" => "5"]
    items = try
        get_json("https://api.crossref.org/works?" * query(params))["message"]["items"]
    catch err
        @warn "Crossref lookup failed" doi exception = err
        return nothing
    end
    return isempty(items) ? nothing : normalize_doi(items[1]["DOI"])
end

"""
Preprints listed on the page that now have a published version: `old DOI => new DOI`.
Uses Crossref relations (Copernicus), else a title match among the scanned works (arXiv, …).
Swaps whose target was already proposed once (in `seen`) are not proposed again.
"""
function published_versions(on_page::Set{String}, seen::Set{String}, cfg::Config, works::Vector{Work})
    preprints = sort!([d for d in on_page if preprint_server(d, cfg) !== nothing])
    isempty(preprints) && return Pair{String, String}[]
    titles = Dict{String, String}()
    for chunk in Iterators.partition(preprints, 50)
        for x in openalex_works("doi:" * join(chunk, "|"); select = "doi,title")
            titles[normalize_doi(x["doi"])] = clean_text(something(x["title"], ""))
        end
    end
    swaps = Pair{String, String}[]
    for d in preprints
        new = crossref_published_version(d)
        if new === nothing && haskey(titles, d)
            k = findfirst(works) do w
                !isempty(w.doi) && !is_preprint(w, cfg) &&
                    !any(startswith(w.doi, p) for p in cfg.ignore_doi_prefixes) &&
                    similar_titles(w.title, titles[d])
            end
            k === nothing || (new = works[k].doi)
        end
        new === nothing || new in on_page || new in seen || push!(swaps, d => new)
    end
    return swaps
end
