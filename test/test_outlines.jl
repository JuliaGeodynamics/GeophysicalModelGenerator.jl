using Test, GeophysicalModelGenerator

mktempdir() do dir
    file(name, content) = (path = joinpath(dir, name); write(path, content); path)

    # Map view of the reference setup of Schuler et al. (2025, doi:10.1029/2024GC011996), drawn in Inkscape for geomIO,
    # from Zenodo 10.5281/zenodo.13987488 (CC-BY 4.0), with its embedded background image removed. Its reference path
    # maps the drawing to 0-6000 x 0-3500 km.
    svg = "test_files/MapviewSetup.svg"
    outlines = read_outlines(svg; tolerance = 2.0)
    @test length(outlines) == 35 && count(o -> o.closed, outlines) == 25
    @test outlines[1].name == "contmargin" && outlines[end].name == "weak_Hellenides"
    @test all(o -> all(-100 .< o.x .< 6100) && all(-100 .< o.y .< 3600), outlines)
    @test sum(o -> length(o.x), read_outlines(svg; tolerance = 20.0)) < sum(o -> length(o.x), outlines)

    @test_throws "must be positive" read_outlines(svg; tolerance = -1)
    @test_throws "does not exist" read_outlines(joinpath(dir, "missing.svg"))
    numeric = [(; name = "007", closed = true, x = [0.0, 1.0, 1.0], y = [0.0, 0.0, 1.0])]
    write_outlines(joinpath(dir, "num.csv"), numeric)
    @test read_outlines(joinpath(dir, "num.csv")) == numeric
    overshoot = file("overshoot.svg", """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100"><path id="c" d="M 10,50 C 90,50 -40,50 50,50"/></svg>""")
    @test length(read_outlines(overshoot; bounds = (0, 100, 0, 100), tolerance = 1.0)[1].x) > 2

    # plates in drawing order, so that later outlines are painted over earlier ones
    Grid = CartData(xyz_grid(0:50.0:6000, 0:50.0:3500, -300:25.0:0))
    Phases, Temp = zeros(Int, size(Grid)), zeros(size(Grid))
    for (i, o) in enumerate(outlines)
        o.closed && add_plate!(Phases, Temp, Grid; xlim = o.x, ylim = o.y, zlim = (-100.0, 0.0), phase = ConstantPhase(i))
    end
    phase_at(x, y) = outlines[Phases[findfirst(≈(x), Grid.x.val[:, 1, 1]), findfirst(≈(y), Grid.y.val[1, :, 1]), end]].name
    @test phase_at(2500, 300) == "Africa" && phase_at(2800, 1600) == "Adria" && phase_at(5700, 2850) == "EuropeAsia"

    # a plain SVG, as written by any editor (e.g. Affinity): shapes named by their own or their group's id or data-name,
    # with transforms, mapped by a named reference line or by `bounds` on the viewBox; `defs` are not outlines
    f = file(
        "plain.svg", """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 50">
        <defs><clipPath id="clip"><rect width="100" height="50"/></clipPath></defs>
        <line id="ref" x1="0" y1="50" x2="100" y2="0"/>
        <g id="square" transform="matrix(1,0,0,1,10,0)"><path d="M0,50l10,0 0-10h-10z"/></g>
        <path data-name="trench" d="M50,25h10v-5"/>
        <g transform="rotate(90 50 25)"><rect id="box" x="40" y="20" width="20" height="10"/></g>
        <circle id="disc" cx="50" cy="25" r="10"/>
        <path id="smooth" d="M0,25 C10,20 20,20 30,25 S50,30 60,25"/>
        <path id="quad" d="M0,25 Q25,0 50,25 T100,25"/>
        <path id="redo" d="M0,0 L10,0 L10,10 Z L20,20"/>
        <use href="#square" x="20"/>
        </svg>"""
    )
    o = read_outlines(f; reference = "ref" => ((0, 0), (1000, 500)))
    @test [x.name for x in o] == ["square", "trench", "box", "disc", "smooth", "quad", "redo", "redo", "square"]
    sq, tr, box, disc, smooth, quad, redo1, redo2, clone = o
    @test sum(extrema(smooth.y)) ≈ 500 && sum(extrema(quad.y)) ≈ 500   # S and T mirror the previous control point
    @test redo1.closed && length(redo1.x) == 3 && !redo2.closed && redo2.x ≈ [0, 200] && redo2.y ≈ [500, 300]   # drawing on after Z
    @test clone.x ≈ sq.x .+ 200 && clone.y ≈ sq.y   # <use> of a group
    @test_throws "not in the file" read_outlines(file("use.svg", """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1"><use href="#gone"/></svg>"""); bounds = (0, 1, 0, 1))
    @test sq.closed && sq.x ≈ [100, 200, 200, 100] && sq.y ≈ [0, 0, 100, 100]
    @test !tr.closed && tr.x ≈ [500, 600, 600] && tr.y ≈ [250, 250, 300]
    @test box.closed && sort(box.x) ≈ [450, 450, 550, 550] && sort(box.y) ≈ [150, 150, 350, 350]
    @test disc.closed && all(abs.(hypot.(disc.x .- 500, disc.y .- 250) .- 100) .<= 1)   # the arcs of the circle
    @test all(((a, b),) -> a.x ≈ b.x && a.y ≈ b.y, zip(read_outlines(f; bounds = (0, 1000, 0, 500))[2:end], o))
    @test_throws "no shape with a CoordRef" read_outlines(f)
    @test_throws "either `bounds` or `reference`" read_outlines(f; bounds = (0, 1, 0, 1), reference = "ref" => ((0, 0), (1, 1)))
    @test_throws "unsupported SVG path command" read_outlines(file("bad.svg", """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1"><path d="M 0,0 B 1,1"/></svg>"""); bounds = (0, 1, 0, 1))
    @test_throws "unsupported outline format" read_outlines(joinpath(dir, "x.blend"))
    @test_throws "flags of an SVG arc" read_outlines(file("arc.svg", """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1"><path d="M 0,0 a 5 5 0 0110 0"/></svg>"""); bounds = (0, 1, 0, 1))
    @test_throws "symbols with a viewBox are not supported" read_outlines(
        file("symbol.svg", """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1"><symbol id="s" viewBox="0 0 2 2"><rect width="1" height="1"/></symbol><use href="#s"/></svg>"""); bounds = (0, 1, 0, 1)
    )
    @test_throws "Z takes no numbers" read_outlines(file("z.svg", """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1"><path d="M0,0 L1,1 Z 5,5"/></svg>"""); bounds = (0, 1, 0, 1))

    # no consecutive duplicate points
    dup = read_outlines(
        file(
            "dup.svg", """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">
            <rect id="pill" width="20" height="40" rx="10"/><path id="repeat" d="M0,0 L10,0 L10,0 L20,5"/></svg>"""
        ); bounds = (0, 100, 0, 100)
    )
    @test all(o -> all(i -> (o.x[i], o.y[i]) != (o.x[i + 1], o.y[i + 1]), 1:(length(o.x) - 1)), dup)
    @test dup[2].x ≈ [0, 10, 20]

    # Inkscape writes an id on every element: a group's label names the shapes in it, unless they have their own label
    inkscape = file(
        "inkscape.svg", """<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" viewBox="0 0 10 10">
        <g inkscape:label="Adria" id="layer1"><path id="path123" d="M0,0 L1,1"/><path id="path124" inkscape:label="trench" d="M0,0 L1,1"/></g>
        <path id="path125" d="M0,0 L1,1"/></svg>"""
    )
    @test [o.name for o in read_outlines(inkscape; bounds = (0, 10, 0, 10))] == ["Adria", "trench", "path125"]

    # names: a group's id does not override a shape's own id, unlike labels and Inkscape layers; a <use> keeps the name
    # of its target unless it has a label; hidden elements are skipped
    svgfile(name, body; viewbox = "0 0 100 100") = file(name, """<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" viewBox="$viewbox">$body</svg>""")
    names = svgfile(
        "names.svg", """<g id="g12"><path id="trench" d="M0,0 L1,1"/><path d="M0,0 L1,1"/></g>
        <g inkscape:groupmode="layer" id="layer7"><path id="path1" d="M0,0 L1,1"/></g>
        <g inkscape:groupmode="layer" id="layer8" inkscape:label="Africa"><path id="path2" d="M0,0 L1,1"/></g>
        <use id="use12" href="#trench"/><use id="use13" inkscape:label="trench2" href="#trench"/>
        <g inkscape:groupmode="layer" inkscape:label="hidden" style="fill:none; display : none"><path id="h1" d="M0,0 L1,1"/></g>
        <path id="h2" display="none" d="M0,0 L1,1"/>"""
    )
    @test [o.name for o in read_outlines(names; bounds = (0, 100, 0, 100))] == ["trench", "g12", "layer7", "Africa", "trench", "trench2"]

    # an open path that returns to its start is closed
    ret = only(read_outlines(svgfile("return.svg", """<path id="r" d="M0,0 L10,0 L10,10 L0,0"/>"""); bounds = (0, 100, 0, 100)))
    @test ret.closed && ret.x ≈ [0, 10, 10] && ret.y ≈ [100, 100, 90]

    # transforms of nested groups compose: translate(10) after scale(2,3), and skewX(45) and scale(2) on a shape
    tf = read_outlines(
        svgfile(
            "transforms.svg", """<g transform="translate(10)"><g transform="scale(2,3)"><path id="nested" d="M1,1 L2,1"/></g></g>
            <path id="skew" transform="skewX(45)" d="M0,10 L10,10"/><path id="scaled" transform="scale(2)" d="M1,1 L2,2"/>"""
        ); bounds = (0, 100, 0, 100)
    )
    @test tf[1].x ≈ [12, 14] && tf[1].y ≈ [97, 97]
    @test tf[2].x ≈ [10, 20] && tf[2].y ≈ [90, 90]
    @test tf[3].x ≈ [2, 4] && tf[3].y ≈ [98, 96]

    # the large arc of an ellipse with radii 20 and 10, rotated by 30°, centered at (50, 50), from angle 0 to 90° the
    # long way round, i.e. through angles 0 to -270°
    c, s = cosd(30), sind(30)
    p0, p1 = (50 + 20c, 50 + 20s), (50 - 10s, 50 + 10c)
    arc = only(read_outlines(svgfile("ellipse.svg", """<path id="arc" d="M$(p0[1]),$(p0[2]) A 20 10 30 1 0 $(p1[1]),$(p1[2])"/>"""); bounds = (0, 100, 0, 100), tolerance = 0.01))
    u, v = (c .* (arc.x .- 50) .+ s .* (50 .- arc.y)) ./ 20, (-s .* (arc.x .- 50) .+ c .* (50 .- arc.y)) ./ 10
    θ = atand.(v, u)
    @test all(abs.(hypot.(u, v) .- 1) .< 1.0e-3) && !arc.closed
    @test [arc.x[1], 100 - arc.y[1]] ≈ collect(p0) && [arc.x[end], 100 - arc.y[end]] ≈ collect(p1)
    @test !any(5 .< θ .< 85) && any(θ .< -175) && any(-95 .< θ .< -85) && length(θ) > 10

    # relative commands equal their absolute counterparts
    rel = read_outlines(
        svgfile(
            "relative.svg", """<path id="abs" d="M10,10 C10,20 20,20 20,10 S30,0 30,10 M10,50 Q20,30 30,50 T50,50 M10,80 A10,5 30 0 1 30,80"/>
            <path id="rel" d="M10,10 c0,10 10,10 10,0 s10,-10 10,0 M10,50 q10,-20 20,0 t20,0 M10,80 a10,5 30 0 1 20,0"/>"""
        ); bounds = (0, 100, 0, 100)
    )
    @test length(rel) == 6 && all(i -> rel[i].x ≈ rel[i + 3].x && rel[i].y ≈ rel[i + 3].y, 1:3)

    # malformed or empty drawings
    bad(body; kw...) = read_outlines(svgfile("bad.svg", body); bounds = (0, 100, 0, 100), kw...)
    @test_throws "refers to itself" bad("""<g id="a"><use href="#a"/></g>""")
    @test_throws "no outlines" bad("""<g id="a" display="none"><path d="M0,0 L1,1"/></g>""")
    @test_throws "root element" read_outlines(file("html.svg", """<html><path d="M0,0 L1,1"/></html>"""); bounds = (0, 1, 0, 1))
    @test_throws "must start with M or m" bad("""<path id="p" d="L0,0 L1,1"/>""")
    @test_throws "path command L of \"p\" needs 2 numbers" bad("""<path id="p" d="M0,0 L1"/>""")
    @test_throws "path command L of \"p\" needs 2 numbers" bad("""<path id="p" d="M0,0 L1 Z"/>""")
    @test_throws "\"#\", which is neither a command nor a number" bad("""<path id="p" d="M0,0 L1,#2"/>""")
    @test_throws "attribute width of <rect> must be 1 number(s) without units" bad("""<rect width="10mm" height="1"/>""")
    @test_throws "viewBox" read_outlines(svgfile("vb.svg", """<path d="M0,0 L1,1"/>"""; viewbox = "0 0 100mm 100"); bounds = (0, 1, 0, 1))
    @test_throws "malformed SVG transform" bad("""<path transform="translate(1,2" d="M0,0 L1,1"/>""")
    @test_throws "arguments of rotate" bad("""<path transform="rotate(1,2)" d="M0,0 L1,1"/>""")
    @test_throws "nonzero extent" bad("""<path d="M0,0 L1,1"/>"""; bounds = (0, 0, 0, 100))
    coordref(ref) = read_outlines(svgfile("coordref.svg", """<path CoordRef="$ref" d="M0,100 L100,0"/><path d="M0,0 L1,1"/>"""))
    @test only(coordref("0,0 10,20")).x ≈ [0, 0.1]
    hiddenref = read_outlines(svgfile("hiddenref.svg", """<g style="display:none"><path CoordRef="0,0 10,20" d="M0,100 L100,0"/></g><path id="v" d="M0,0 L1,1"/>"""))
    @test only(hiddenref).name == "v" && only(hiddenref).x ≈ [0, 0.1]
    @test_throws "CoordRef" coordref("0,0 10")
    @test_throws "nonzero extent" coordref("0,0 0,10")

    # CSV round trip, with comments
    outlines = [
        (; name = "plate A", closed = true, x = [0.0, 10, 10], y = [0.0, 0, 5]), (; name = "trench", closed = false, x = [1.0, 2], y = [3.0, 4]),
        (; name = "trench", closed = false, x = [5.0, 6], y = [7.0, 8]),
    ]
    f = joinpath(dir, "outlines.csv")
    write_outlines(f, outlines; comment = "test outlines\nunits: km")
    @test startswith(read(f, String), "# test outlines\n# units: km\nid,name,closed,x,y\n")
    @test read_outlines(f) == outlines
    @test_throws "must not contain commas" write_outlines(f, [(; name = "a,b", closed = false, x = [0.0], y = [0.0])])
    @test_throws "columns of" read_outlines(file("bad.csv", "name,x,y\na,1,2\n"))
    @test_throws "must be 0 or 1, got \"true\"" read_outlines(file("closed.csv", "id,name,closed,x,y\n1,a,true,0,0\n1,a,true,1,1\n"))
    @test_throws "at least 2" write_outlines(f, [(; name = "a", closed = false, x = [0.0], y = [0.0])])
    @test_throws "quotes" write_outlines(f, [(; name = "a\"b", closed = false, x = [0.0, 1.0], y = [0.0, 1.0])])
end
