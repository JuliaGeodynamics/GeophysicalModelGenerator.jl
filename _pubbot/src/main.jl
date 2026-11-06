function load_seen(path::AbstractString = SEEN)
    isfile(path) || return Set{String}()
    return Set(normalize_doi(l) for l in eachline(path) if !isempty(strip(l)) && !startswith(l, "#"))
end

function save_seen(path::AbstractString, seen::Set{String})
    return open(path, "w") do io
        println(io, "# DOIs already evaluated by the publication bot (_pubbot); they are not proposed again.")
        println(io, "# Delete a line to have that DOI evaluated again by the next scan.")
        foreach(d -> println(io, d), sort!(collect(seen)))
    end
end

"Fill in the formatted entry; demote to review when metadata cannot be fetched."
function add_entry!(a::Assessment, cfg::Config)
    csl = get_csl(a.work.doi)
    if csl === nothing
        a.note = "DOI metadata unavailable, add manually"
        a.bucket == :strong && (a.bucket = :review)
        return a
    end
    a.entry = entry_line(csl, cfg)
    a.year = something(csl_year(csl), year(a.work.date))
    return a
end

"""
    weekly(; since, dry_run, report_path, page_path, seen_path)

Scan for new publications, edit `publications.md` and `seen_dois.txt` (unless `dry_run`),
and write the pull request description to `report_path` (stdout if `nothing`).
"""
function weekly(;
        since::Union{Date, Nothing} = nothing, dry_run::Bool = false,
        report_path::Union{AbstractString, Nothing} = nothing,
        page_path::AbstractString = PAGE, seen_path::AbstractString = SEEN
    )
    cfg = load_config()
    since = something(since, today() - Month(cfg.lookback_months))
    page = load_page(page_path)
    on_page = page_dois(page)
    seen = load_seen(seen_path)

    refs = package_ids(cfg)
    works, mentions = fetch_works(cfg, refs; since)
    candidates = select_candidates(works, cfg; on_page, seen)
    @info "Fetched works" since records = length(works) candidates = length(candidates)

    texts = asyncmap(fulltext, candidates; ntasks = 6)
    assessments = [assess(w, cfg, refs, t, mentions) for (w, t) in zip(candidates, texts)]
    drop_published_preprints!(assessments, cfg; on_page, seen)
    sort!(assessments; by = a -> a.work.date, rev = true)
    foreach(a -> a.bucket == :none || add_entry!(a, cfg), assessments)

    swaps = NamedTuple[]
    for (old, new) in published_versions(on_page, seen, cfg, works)
        csl = get_csl(new)
        csl === nothing && continue
        title = clean_text(first_string(get(csl, "title", "")))
        push!(swaps, (; old, new, title, entry = entry_line(csl, cfg), year = something(csl_year(csl), year(today()))))
    end

    # Oldest first, so the newest paper ends up at the top of its year.
    for a in Iterators.reverse(assessments)
        a.bucket == :strong && insert_entry!(page, a.year, a.entry)
    end
    for s in swaps
        remove_entry!(page, s.old)
        insert_entry!(page, s.year, s.entry)
    end

    changed = !isempty(swaps) || any(a -> a.bucket != :none, assessments)
    if changed && !dry_run
        save_page(page)
        for a in assessments
            push!(seen, a.work.doi)
            union!(seen, a.work.other_dois)
        end
        foreach(s -> push!(seen, s.new), swaps)
        save_seen(seen_path, seen)
    end
    report = weekly_report(assessments, swaps; since, nworks = length(works))
    report_path === nothing ? print(report) : write(report_path, report)
    set_output("changed", changed && !dry_run)
    return (; assessments, swaps, changed)
end

"""
    add_dois(dois; issue, report_path, page_path, seen_path)

Format and insert the publications with the given DOIs.
"""
function add_dois(
        dois::AbstractVector{<:AbstractString}; issue = nothing,
        report_path::Union{AbstractString, Nothing} = nothing,
        page_path::AbstractString = PAGE, seen_path::AbstractString = SEEN
    )
    cfg = load_config()
    page = load_page(page_path)
    on_page = page_dois(page)
    seen = load_seen(seen_path)
    added, skipped = Tuple{String, String}[], Tuple{String, String}[]
    for doi in unique(normalize_doi.(dois))
        if doi in on_page
            push!(skipped, (doi, "already listed"))
            continue
        end
        csl = get_csl(doi)
        if csl === nothing
            push!(skipped, (doi, "no DOI metadata with authors"))
            continue
        end
        entry = entry_line(csl, cfg)
        insert_entry!(page, something(csl_year(csl), year(today())), entry)
        push!(on_page, doi)
        push!(seen, doi)
        push!(added, (doi, entry))
    end
    changed = !isempty(added)
    if changed
        save_page(page)
        save_seen(seen_path, seen)
    end
    report = add_report(added, skipped; issue)
    report_path === nothing ? print(report) : write(report_path, report)
    set_output("changed", changed)
    return (; added, skipped, changed)
end

# Pull request descriptions

link(doi) = "[$doi]($(doi_url(doi)))"
escape_md(s) = replace(s, r"([*_`\[\]<>])" => s"\\\1")

function summary_line(a::Assessment)
    evidence = isempty(a.evidence) ? "none" : join(a.evidence, "; ")
    return "**$(escape_md(a.work.title))** · $(link(a.work.doi))  \n  evidence: $evidence · " *
        (a.pdf ? "PDF checked" : "no PDF")
end

details(io, summary, body) = println(io, "<details><summary>$summary</summary>\n\n", body, "\n</details>\n")

# shortcut: no truncation to GitHub's 65536-character PR body limit; add one if a scan's ❌ list ever hits it.
function weekly_report(assessments::Vector{Assessment}, swaps; since::Date, nworks::Int)
    strong = filter(a -> a.bucket == :strong, assessments)
    review = filter(a -> a.bucket == :review, assessments)
    none = filter(a -> a.bucket == :none, assessments)
    io = IOBuffer()
    println(io, "## Publication scan ($(today()))\n")
    println(
        io, "Works published since $since by the developers, citing the package paper or naming the package: ",
        "$nworks records, $(length(assessments)) new papers after removing duplicates and already evaluated DOIs.\n"
    )

    println(io, "### ✅ Added to publications.md ($(length(strong)))\n")
    isempty(strong) || println(io, "To drop an entry, delete its line in `publications.md` before merging; it will not be proposed again.\n")
    foreach(a -> println(io, "- ", summary_line(a)), strong)
    println(io)

    if !isempty(swaps)
        println(io, "### 🔁 Preprints replaced by their published version ($(length(swaps)))\n")
        foreach(s -> println(io, "- $(link(s.old)) → $(link(s.new)): $(escape_md(s.title))"), swaps)
        println(io)
    end

    println(io, "### 🟡 Needs review ($(length(review)))\n")
    if !isempty(review)
        println(io, "To add one, paste its line into `publications.md` on this branch.\n")
        body = join(
            ("- " * summary_line(a) * "\n" * (a.entry === nothing ? "  _$(a.note)_" : "  ```md\n  $(a.entry)\n  ```") for a in review),
            "\n"
        )
        details(io, "Show $(length(review)) papers", body)
    end

    println(io, "### ❌ No signal ($(length(none)))\n")
    if !isempty(none)
        body = join(("- " * summary_line(a) * (isempty(a.note) ? "" : " · $(a.note)") for a in none), "\n")
        details(io, "Show $(length(none)) papers", body)
    end
    print(io, "---\nMerging records every DOI above in `_pubbot/seen_dois.txt`. Generated by `_pubbot/weekly.jl`.\n")
    return String(take!(io))
end

function add_report(added, skipped; issue = nothing)
    io = IOBuffer()
    println(io, "## Add publications\n")
    for (doi, entry) in added
        println(io, "- ✅ $(link(doi))\n  ```md\n  $entry\n  ```")
    end
    for (doi, reason) in skipped
        println(io, "- ⏭️ $(link(doi)): $reason")
    end
    isempty(added) && isempty(skipped) && println(io, "No DOI found in the request.")
    issue === nothing || println(io, "\nCloses #$issue")
    return String(take!(io))
end
