# Add publications by DOI.
#
#   julia --project=_pubbot _pubbot/add_doi.jl [--report FILE] DOI...
#
# Without DOI arguments, DOIs are read from the ISSUE_BODY environment variable
# (the "Add a publication" issue form); ISSUE_NUMBER adds "Closes #N" to the report.

using PubBot

report_path = nothing
dois = String[]
let args = copy(ARGS)
    while !isempty(args)
        arg = popfirst!(args)
        arg == "--report" ? (global report_path = popfirst!(args)) : append!(dois, PubBot.extract_dois(arg))
    end
end
if isempty(dois)
    body = get(ENV, "ISSUE_BODY", "")
    # Only the "DOI(s)" field of the form, not DOIs mentioned in the note.
    section = match(r"###\s*DOI\(s\)\s*\n(.*?)(?=\n###\s|\z)"s, body)
    append!(dois, PubBot.extract_dois(section === nothing ? body : section[1]))
end
issue = get(ENV, "ISSUE_NUMBER", "")

PubBot.add_dois(dois; report_path, issue = isempty(issue) ? nothing : issue)
