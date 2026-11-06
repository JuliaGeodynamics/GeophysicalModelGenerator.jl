# Entries of publications.md are `- ` lines under `### YYYY` headings, newest year first.

struct Page
    path::String
    lines::Vector{String}
end

load_page(path::AbstractString = PAGE) = Page(path, readlines(path))
save_page(p::Page) = write(p.path, rstrip(join(p.lines, "\n")) * "\n")

page_dois(p::Page) = Set(extract_dois(join(p.lines, "\n")))

function year_heading(line::AbstractString)
    m = match(r"^###\s+(\d{4})\s*$", line)
    return m === nothing ? nothing : parse(Int, m[1])
end

is_entry(line::AbstractString) = startswith(line, "- ")

"Insert `entry` at the top of the `### year` section, creating the section if needed."
function insert_entry!(p::Page, year::Integer, entry::AbstractString)
    r = eachindex(p.lines)
    for i in r
        y = year_heading(p.lines[i])
        y === nothing && continue
        if y == year
            k = i + 1
            k <= last(r) && isempty(strip(p.lines[k])) && (k += 1)
            block = k <= last(r) && is_entry(p.lines[k]) ? [entry] : [entry, ""]
            splice!(p.lines, k:(k - 1), block)
            return p
        elseif y < year
            splice!(p.lines, i:(i - 1), ["### $year", "", entry, ""])
            return p
        end
    end
    # Oldest year so far: append at the end of the list.
    k = last(r) + 1
    block = ["### $year", "", entry, ""]
    k > 1 && !isempty(strip(p.lines[k - 1])) && pushfirst!(block, "")
    splice!(p.lines, k:(k - 1), block)
    return p
end

"Remove the entry citing `doi`, and its year heading if the section becomes empty."
function remove_entry!(p::Page, doi::AbstractString)
    i = findfirst(l -> is_entry(l) && doi in extract_dois(l), p.lines)
    i === nothing && return false
    deleteat!(p.lines, i)
    remove_empty_sections!(p)
    return true
end

function remove_empty_sections!(p::Page)
    r = eachindex(p.lines)
    heads = [i for i in r if year_heading(p.lines[i]) !== nothing]
    for (n, h) in Iterators.reverse(collect(enumerate(heads)))
        stop = n < length(heads) ? heads[n + 1] - 1 : last(r)
        any(is_entry, p.lines[(h + 1):stop]) && continue
        # Drop the heading and the blank lines that follow it.
        e = h + 1
        while e <= stop && isempty(strip(p.lines[e]))
            e += 1
        end
        deleteat!(p.lines, h:(e - 1))
    end
    return p
end
