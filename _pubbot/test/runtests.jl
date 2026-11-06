using PubBot, Test, Dates
using PubBot: normalize_doi, extract_dois, similar_titles, initials, entry_line, load_config,
    Page, page_dois, insert_entry!, remove_entry!, Work, assess, select_candidates

const CFG = load_config()

@testset "DOIs" begin
    @test normalize_doi("https://doi.org/10.1029/2025JB031240") == "10.1029/2025jb031240"
    @test normalize_doi("doi: 10.5194/gmd-15-5757-2022.") == "10.5194/gmd-15-5757-2022"
    @test extract_dois("see 10.1234/abc?x=\"1\" and <script>") == ["10.1234/abc"]
end

@testset "Titles" begin
    @test similar_titles("Thermal modelling of magma reservoirs", "Thermal modeling of magma reservoirs")
    @test similar_titles(
        "Thermal evolution of crustal magma reservoirs and their zircon record",
        "Thermal evolution of crustal magma reservoirs and their zircon record with MTK (version 1.0.0)"
    )
    @test !similar_titles("Magma reservoirs", "Magma reservoirs in the upper crust of Chile")  # too short to extend
    @test !similar_titles("Sill emplacement, Part I: theory", "Sill emplacement, Part II: theory")
end

@testset "Formatting" begin
    @test initials("Boris J. P.") == "B. J. P."
    @test initials("Jean-Pierre") == "J.-P."
    @test PubBot.clean_text("Thermo‐Mechanical <i>Controls</i> &amp; 10") == "Thermo-Mechanical Controls & 10"
    authors = [
        Dict("family" => "Kaus", "given" => "Boris J. P."),
        Dict("family" => "Monserrat", "non-dropping-particle" => "de", "given" => "Albert"),
        Dict("family" => "Aellig", "given" => "Pascal"),
    ]
    article = Dict(
        "DOI" => "10.1130/GES02563.1", "author" => authors, "title" => "Zircon age spectra",
        "container-title" => "Geosphere", "volume" => "19", "issue" => "4", "page" => "1006-1031",
        "issued" => Dict("date-parts" => [[2023, 6]])
    )
    @test entry_line(article, CFG) ==
        "- Kaus, B. J. P., de Montserrat, A., and Aellig, P. (2023). **Zircon age spectra.** Geosphere, 19(4), 1006–1031. [https://doi.org/10.1130/GES02563.1](https://doi.org/10.1130/GES02563.1)"
    preprint = Dict(
        "DOI" => "10.5194/egusphere-2025-1", "author" => authors[3:3], "title" => "Sills?",
        "issued" => Dict("date-parts" => [[2025]])
    )
    @test entry_line(preprint, CFG) ==
        "- Aellig, P. (2025). **Sills?** EGUsphere (preprint). [https://doi.org/10.5194/egusphere-2025-1](https://doi.org/10.5194/egusphere-2025-1)"
end

const PAGE_TEXT = """
# Publications

Intro.

### 2026

- A (2026). [https://doi.org/10.1234/a](https://doi.org/10.1234/a)

### 2024

- B (2024). [https://doi.org/10.1234/b](https://doi.org/10.1234/b)
"""

page() = Page("", split(chomp(PAGE_TEXT), "\n"))
text(p) = rstrip(join(p.lines, "\n")) * "\n"  # as `save_page` writes it

@testset "Page editing" begin
    p = page()
    @test page_dois(p) == Set(["10.1234/a", "10.1234/b"])
    insert_entry!(p, 2026, "- new2026")
    @test occursin("### 2026\n\n- new2026\n- A (2026)", text(p))
    insert_entry!(p, 2025, "- new2025")
    @test occursin("(https://doi.org/10.1234/a)\n\n### 2025\n\n- new2025\n\n### 2024", text(p))
    insert_entry!(p, 2027, "- new2027")
    @test occursin("Intro.\n\n### 2027\n\n- new2027\n\n### 2026", text(p))
    insert_entry!(p, 2020, "- new2020")
    @test endswith(text(p), "(https://doi.org/10.1234/b)\n\n### 2020\n\n- new2020\n")

    p = page()
    @test remove_entry!(p, "10.1234/b")
    @test text(p) == replace(PAGE_TEXT, "\n### 2024\n\n- B (2024). [https://doi.org/10.1234/b](https://doi.org/10.1234/b)\n" => "")
    @test !remove_entry!(p, "10.1234/x")

    # The real page: every entry has a DOI, and the introduction names none.
    real = PubBot.load_page()
    @test length(page_dois(real)) == count(PubBot.is_entry, real.lines)
end

work(; title = "A paper", doi = "10.1234/x", refs = Set{String}(), abstract = "", type = "article", authors = ["Jane Doe"]) =
    Work(; id = "W1", doi, title, date = Date(2026, 1, 1), type, authors, refs, abstract)

@testset "Classification" begin
    refs = Set(["W98"])
    none = Set{String}()
    bucket(w, text = nothing, mentions = none) = assess(w, CFG, refs, text, mentions).bucket
    # Named in the PDF or the abstract: listed.
    @test bucket(work(), "Models were run with GeophysicalModelGenerator.jl.") == :strong
    @test bucket(work(; abstract = "using GeophysicalModelGenerator")) == :strong
    @test bucket(work(), "Models were run with GeophysicalModelGenerator.jl.", Set(["W1"])) == :strong
    @test bucket(work(; doi = "10.21105/joss.06763")) == :strong
    # Found only by OpenAlex's stemmed full-text index: a human decides.
    @test bucket(work(), nothing, Set(["W1"])) == :review
    @test bucket(work(), "A geophysical model generator.", Set(["W1"])) == :review
    # Citing the package paper alone: a human decides.
    @test bucket(work(; refs)) == :review
    @test bucket(work(; refs), "Rheology as in the cited package paper.") == :review
    # Neither named nor cited, or not a whole word.
    @test bucket(work(), "A Julia GPU solver") == :none
    @test bucket(work(), "GeophysicalModelGeneratorX") == :none
    @test bucket(work(), "A geophysical model generator.") == :none
end

@testset "Candidate selection" begin
    pre = work(; doi = "10.5194/egusphere-2025-1", title = "Thermal modelling of magma reservoirs", type = "preprint")
    pub = work(; doi = "10.5194/gmd-19-1-2026", title = "Thermal modeling of magma reservoirs")
    zen = work(; doi = "10.5281/zenodo.1", title = "Thermal modeling of magma reservoirs")
    thanks = work(; doi = "10.1029/thanks", title = "Thank You to Our 2025 Peer Reviewers")
    c = select_candidates([pre, pub, zen, thanks], CFG; on_page = Set{String}(), seen = Set{String}())
    @test length(c) == 1 && c[1].doi == pub.doi && c[1].other_dois == [pre.doi]
    @test isempty(select_candidates([pre, pub], CFG; on_page = Set([pre.doi]), seen = Set{String}()))
    @test isempty(select_candidates([pre, pub], CFG; on_page = Set{String}(), seen = Set([pub.doi])))
end
