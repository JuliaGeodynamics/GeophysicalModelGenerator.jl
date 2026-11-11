# Crossref's recommended DOI pattern, without parentheses so markdown links parse cleanly.
const DOI_RE = r"10\.\d{4,9}/[-._;/:A-Za-z0-9]+"

"Bare lowercase DOI, without resolver prefix or trailing punctuation."
function normalize_doi(s::AbstractString)
    d = replace(strip(s), r"^(https?://)?(dx\.)?doi\.org/"i => "", r"^doi:\s*"i => "")
    return lowercase(rstrip(d, ['.', ',', ';', ':']))
end

"All DOIs found in `text`, normalized, in order of appearance."
extract_dois(text::AbstractString) = unique!([normalize_doi(m.match) for m in eachmatch(DOI_RE, text)])

doi_url(doi::AbstractString) = "https://doi.org/" * doi

"Collapse whitespace and drop HTML/JATS tags from metadata strings."
function clean_text(s::AbstractString)
    s = replace(s, r"<[^>]+>" => "")
    s = replace(s, "&amp;" => "&", "&lt;" => "<", "&gt;" => ">", "&quot;" => "\"", "‐" => "-", "‑" => "-")
    return strip(replace(s, r"\s+" => " "))
end

"First string of a metadata field that may be a string, a list or missing."
function first_string(x)
    x isa AbstractString && return String(x)
    x isa AbstractVector && !isempty(x) && x[1] isa AbstractString && return String(x[1])
    return ""
end

"Accent-free lowercase words, for fuzzy comparisons."
function words(s::AbstractString)
    t = Unicode.normalize(String(s); stripmark = true, casefold = true)
    return [m.match for m in eachmatch(r"[a-z0-9]+", t)]
end

"Title words, with British spelling mapped to American (\"modelling\" → \"modeling\")."
function title_words(s::AbstractString)
    american(w) = replace(w, r"ll(ing|ed)$" => s"l\1", r"is(e|ed|es|ing|ation|ations)$" => s"iz\1", r"our$" => "or")
    return Set(american.(words(s)))
end

"""
Titles of two versions of a paper: nearly the same words, or one title extending the other
(the published version often gains a suffix like "(version 1.0.0)"). Titles whose numbers or
roman numerals differ ("Part I" / "Part II") never match.
"""
function similar_titles(a::AbstractString, b::AbstractString)
    wa, wb = title_words(a), title_words(b)
    (isempty(wa) || isempty(wb)) && return false
    numbers(ws) = filter(w -> occursin(r"^(\d+|i|ii|iii|iv)$", w), ws)
    na, nb = numbers(wa), numbers(wb)
    isempty(na) || isempty(nb) || na == nb || return false
    common, shorter = length(wa ∩ wb), min(length(wa), length(wb))
    return common >= 0.8 * length(wa ∪ wb) || shorter >= 6 && common >= 0.9 * shorter
end

"Append `key=value` to the GitHub Actions step outputs, when running in Actions."
function set_output(key, value)
    path = get(ENV, "GITHUB_OUTPUT", "")
    isempty(path) && return
    return open(io -> println(io, "$key=$value"), path, "a")
end

struct Config
    match::String
    pattern::Regex
    dois::Vector{String}
    developers::Vector{String}
    lookback_months::Int
    types::Vector{String}
    ignore_doi_prefixes::Vector{String}
    ignore_titles::Vector{Regex}
    preprint_servers::Vector{Pair{String, String}}
    author_fixes::Vector{Pair{String, String}}
end

function load_config()
    c = TOML.parsefile(joinpath(BOTDIR, "config.toml"))
    servers = sort!([String(k) => String(v) for (k, v) in c["preprint_servers"]]; by = s -> -length(first(s)))  # longest prefix wins
    return Config(
        c["match"], Regex("(?<![A-Za-z0-9])\\Q" * c["match"] * "\\E(?![A-Za-z0-9])"),
        normalize_doi.(c["dois"]), String.(c["developers"]), c["lookback_months"], String.(c["types"]),
        String.(c["ignore_doi_prefixes"]), [Regex(r, "i") for r in c["ignore_titles"]], servers,
        [String(k) => String(v) for (k, v) in c["author_fixes"]]
    )
end

"Name of the preprint server hosting `doi`, or `nothing` for other DOIs."
function preprint_server(doi::AbstractString, cfg::Config)
    k = findfirst(s -> startswith(doi, first(s)), cfg.preprint_servers)
    return k === nothing ? nothing : last(cfg.preprint_servers[k])
end

# HTTP

const USER_AGENT = "GeophysicalModelGenerator-pubbot (+https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl)"
const DOWNLOADER = Ref{Union{Downloads.Downloader, Nothing}}(nothing)

# Curl gives up after 20 s without data by default; OpenAlex can pause longer on large pages.
function downloader()
    if DOWNLOADER[] === nothing
        d = Downloads.Downloader()
        d.easy_hook = (easy, _) -> Downloads.Curl.setopt(easy, Downloads.Curl.CURLOPT_LOW_SPEED_TIME, 60)
        DOWNLOADER[] = d
    end
    return DOWNLOADER[]
end

"""
    http_get(url; headers, method, timeout, retries) -> (status, body, final_url)

Follows redirects and retries on rate limiting, server errors and network failures.
Returns status 0 when the request could not be made at all.
"""
function http_get(
        url::AbstractString; headers = Pair{String, String}[], method = "GET",
        timeout = 60, retries = 3
    )
    hdrs = ["User-Agent" => USER_AGENT; headers]
    for attempt in 1:retries
        io = IOBuffer()
        resp = try
            Downloads.request(url; method, output = io, headers = hdrs, timeout, throw = false, downloader = downloader())
        catch err
            err
        end
        if resp isa Downloads.Response
            retry = resp.status == 429 || resp.status >= 500
            retry && attempt < retries && (sleep(2^attempt); continue)
            return (resp.status, take!(io), resp.url)
        end
        attempt < retries && (sleep(2^attempt); continue)
        @warn "Request failed" url exception = resp
    end
    return (0, UInt8[], String(url))
end

function get_json(url::AbstractString)
    status, body, _ = http_get(url)
    status == 200 || error("GET $url returned status $status")
    return JSON.parse(String(body))
end

function urlencode(s::AbstractString)
    io = IOBuffer()
    for b in codeunits(s)
        if UInt8('A') <= b <= UInt8('Z') || UInt8('a') <= b <= UInt8('z') ||
                UInt8('0') <= b <= UInt8('9') || b in (UInt8('-'), UInt8('_'), UInt8('.'), UInt8('~'))
            write(io, b)
        else
            print(io, '%', uppercase(string(b; base = 16, pad = 2)))
        end
    end
    return String(take!(io))
end

query(params) = join((k * "=" * urlencode(string(v)) for (k, v) in params), "&")

"CSL-JSON metadata for a DOI (Crossref, DataCite, …), or `nothing` if it does not resolve or lists no authors."
function get_csl(doi::AbstractString)
    status, body, _ = http_get(doi_url(doi); headers = ["Accept" => "application/vnd.citationstyles.csl+json"])
    status == 200 || return nothing
    csl = try
        JSON.parse(String(body))
    catch
        return nothing
    end
    csl isa AbstractDict && !isempty(something(get(csl, "author", nothing), [])) || return nothing
    return csl
end
