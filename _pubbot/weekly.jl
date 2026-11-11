# Weekly publication scan.
#
#   julia --project=_pubbot _pubbot/weekly.jl [--dry-run] [--since YYYY-MM-DD] [--report FILE]
#
# --dry-run  only print the report, leave publications.md and seen_dois.txt untouched

using PubBot, Dates

function parse_args(args)
    opts = Dict{Symbol, Any}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--dry-run"
            opts[:dry_run] = true
        elseif arg == "--since"
            opts[:since] = Date(args[i += 1])
        elseif arg == "--report"
            opts[:report_path] = args[i += 1]
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end
    return opts
end

PubBot.weekly(; parse_args(ARGS)...)
