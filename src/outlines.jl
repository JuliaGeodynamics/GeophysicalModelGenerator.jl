# Reading and writing named map-view outlines (polygons and polylines), e.g. plate boundaries and trenches

using LightXML

export read_outlines, write_outlines

const Outline = @NamedTuple{name::String, closed::Bool, x::Vector{Float64}, y::Vector{Float64}}

"""
    outlines = read_outlines(filename; tolerance = 1.0, bounds = nothing, reference = nothing)

Reads named map-view outlines, such as plate boundaries, weak zones and trenches, from a file. Returns a vector of
named tuples `(; name, closed, x, y)` in the order of the file, which is the drawing order of an SVG file, so that
later outlines can be painted over earlier ones. A closed outline is a polygon, e.g. for
`add_plate!(Phase, Temp, Grid; xlim = o.x, ylim = o.y, ...)`; only an open outline is valid as `Trace` of `add_fault!`
or of a `Trench`, which also reject traces that turn by more than 90° between consecutive segments.

The format follows from the extension of `filename`:
- `.csv`: one row per vertex with the columns `id,name,closed,x,y`, as written by [`write_outlines`](@ref); the rows of
  an outline share its `id`, and `closed` is 0 or 1. Lines starting with `#` are comments.
- `.svg`: every visible shape (`path`, `rect`, `circle`, `ellipse`, `line`, `polyline`, `polygon`) outside `defs`,
  clip paths, masks, patterns, symbols and markers is an outline, so that drawings from any vector graphics program
  (Inkscape, Affinity, Illustrator, ...) can be read. Elements hidden by `display: none`, such as hidden Inkscape
  layers, give no outlines, though they may hold the reference shape, and shapes exported as images, e.g.
  because they are semi-transparent, are not read. An outline is named by its attribute `label` (geomIO, or Inkscape's
  `inkscape:label`) or `data-name` (Illustrator), else by the label of its closest labeled group or Inkscape layer (or
  the layer's `id`), else by its own `id`, else by the `id` of its closest group (the layer name in Affinity). A
  `<use>` with a label names its copy. A path with several subpaths, e.g. a polygon with holes, gives one outline of
  the same name per subpath. A subpath that ends at its start is closed. Curves and arcs are flattened into polylines
  that deviate less than `tolerance` from them, in model units. Drawing coordinates are mapped to model coordinates by one of
  - `reference = name => ((X0, Y0), (X1, Y1))`: the shape `name`, typically a straight line, starts at the model
    coordinates `(X0, Y0)` and ends at `(X1, Y1)`;
  - `bounds = (xmin, xmax, ymin, ymax)`: the model coordinates of the `viewBox` of the drawing (its page);
  - a shape with an attribute `CoordRef = "X0,Y0 X1,Y1"`, as in geomIO, used if neither keyword is given.
  The reference shape is not returned as an outline and must span a nonzero extent in x and y. Lengths with units
  (e.g. `width="10mm"` of a `rect`), arc flags written without a separator before the next number
  (`a 5 5 0 0110 0`) and `<use>` of a `<symbol>` with a `viewBox` are not supported and throw an error.

Example
========
```julia-repl
julia> outlines = read_outlines("MapviewSetup.svg"; tolerance = 4.0);
julia> adria = only(filter(o -> o.name == "Adria", outlines));
julia> add_plate!(Phases, Temp, Grid; xlim = adria.x, ylim = adria.y, zlim = (-110, 0), phase = LithosphericPhases(Layers = [30, 80], Phases = [10, 3, 1]))
```
"""
function read_outlines(filename::AbstractString; tolerance = 1.0, bounds = nothing, reference = nothing)
    tolerance > 0 || throw(ArgumentError("`tolerance` must be positive, got $tolerance"))
    ext = lowercase(splitext(filename)[2])
    !isfile(filename) && ext in (".csv", ".svg") && throw(ArgumentError("outline file $filename does not exist"))
    ext == ".csv" && return read_outlines_csv(filename)
    ext == ".svg" && return read_outlines_svg(filename, tolerance, bounds, reference)
    throw(ArgumentError("unsupported outline format \"$ext\" of $filename; supported are .csv and .svg"))
end

"""
    write_outlines(filename, outlines; comment = nothing)

Writes `outlines`, a vector of named tuples `(; name, closed, x, y)` as returned by [`read_outlines`](@ref), to the
CSV file `filename`, with one row per vertex and the columns `id,name,closed,x,y`. `comment` is written as `#` lines
at the top, e.g. to record the source and units of the outlines.
"""
function write_outlines(filename::AbstractString, outlines; comment = nothing)
    for o in outlines
        occursin(r"[,\"\r\n#]", o.name) && throw(ArgumentError("outline names must not contain commas, quotes, newlines or #, got $(repr(o.name))"))
        length(o.x) == length(o.y) || throw(DimensionMismatch("outline \"$(o.name)\" has $(length(o.x)) x and $(length(o.y)) y coordinates"))
        length(o.x) >= 2 || throw(ArgumentError("outline \"$(o.name)\" has $(length(o.x)) points; an outline needs at least 2"))
    end
    open(filename, "w") do io
        isnothing(comment) || foreach(l -> println(io, "# ", l), split(comment, '\n'))
        println(io, "id,name,closed,x,y")
        for (id, o) in enumerate(outlines), (x, y) in zip(o.x, o.y)
            println(io, id, ',', o.name, ',', Int(o.closed), ',', x, ',', y)
        end
    end
    return nothing
end

function read_outlines_csv(filename)
    data, header = readdlm(filename, ',', String; header = true, comments = true, comment_char = '#')
    strip.(vec(header)) == ["id", "name", "closed", "x", "y"] ||
        throw(ArgumentError("the columns of $filename must be id,name,closed,x,y, got $(join(header, ','))"))
    outlines = Outline[]
    for i in axes(data, 1)
        id, closed = strip(data[i, 1]), strip(data[i, 3])
        closed in ("0", "1") || throw(ArgumentError("the column closed of $filename must be 0 or 1, got \"$closed\" in the row of outline $id"))
        if i == firstindex(data, 1) || id != strip(data[i - 1, 1])
            push!(outlines, (; name = data[i, 2], closed = closed == "1", x = Float64[], y = Float64[]))
        end
        push!(outlines[end].x, parse(Float64, data[i, 4]))
        push!(outlines[end].y, parse(Float64, data[i, 5]))
    end
    return outlines
end

# ---- SVG ---------------------------------------------------------------------------------------------------------------

# 2D affine map (a, b, c, d, e, f): (x, y) -> (a x + c y + e, b x + d y + f), as SVG's matrix(a, b, c, d, e, f)
const Affine = NTuple{6, Float64}
const IDENTITY = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
svg_apply(m::Affine, p) = (m[1] * p[1] + m[3] * p[2] + m[5], m[2] * p[1] + m[4] * p[2] + m[6])
svg_linear(m::Affine, v) = (m[1] * v[1] + m[3] * v[2], m[2] * v[1] + m[4] * v[2])
svg_compose(m::Affine, n::Affine) = (svg_linear(m, (n[1], n[2]))..., svg_linear(m, (n[3], n[4]))..., svg_apply(m, (n[5], n[6]))...)
svg_translation(x, y) = (1.0, 0.0, 0.0, 1.0, x, y)

# the numbers in `s`, separated by whitespace or commas, whose count must be one of `counts`; `what` names `s` in errors
function svg_floats(s, what, counts)
    v = tryparse.(Float64, split(strip(s), r"[\s,]+"; keepempty = false))
    all(!isnothing, v) && length(v) in counts ||
        throw(ArgumentError("$what must be $(join(counts, " or ")) number(s) without units, got \"$s\""))
    return Float64.(v)
end

# the attribute `a` of the element `e` as a number, `default` if it is missing
svg_number(e, a, default = "0") = only(svg_floats(something(attribute(e, a), default), "the attribute $a of <$(name(e))>", 1))

function svg_transform(s)
    isnothing(s) && return IDENTITY
    re = r"(\w+)\s*\(([^)]*)\)"
    occursin(r"^[\s,]*$", replace(s, re => "")) || throw(ArgumentError("malformed SVG transform \"$s\""))
    m = IDENTITY
    for t in eachmatch(re, s)
        f = t[1]
        counts = f == "matrix" ? 6 : f in ("translate", "scale") ? (1, 2) : f == "rotate" ? (1, 3) :
            f in ("skewX", "skewY") ? 1 : throw(ArgumentError("unsupported SVG transform \"$f\""))
        v = svg_floats(t[2], "the arguments of $f in the SVG transform \"$s\"", counts)
        n = if f == "matrix"
            Affine(v)
        elseif f == "translate"
            svg_translation(v[1], get(v, 2, 0.0))
        elseif f == "scale"
            (v[1], 0.0, 0.0, get(v, 2, v[1]), 0.0, 0.0)
        elseif f == "rotate"   # by v[1] degrees, about (v[2], v[3]) if given
            r = (cosd(v[1]), sind(v[1]), -sind(v[1]), cosd(v[1]), 0.0, 0.0)
            length(v) == 3 ? svg_compose(svg_compose(svg_translation(v[2], v[3]), r), svg_translation(-v[2], -v[3])) : r
        elseif f == "skewX"
            (1.0, 0.0, tand(v[1]), 1.0, 0.0, 0.0)
        else   # skewY
            (1.0, tand(v[1]), 0.0, 1.0, 0.0, 0.0)
        end
        m = svg_compose(m, n)
    end
    return m
end

const SVG_SHAPES = ("path", "rect", "circle", "ellipse", "line", "polyline", "polygon")
const SVG_NOT_RENDERED = ("defs", "clipPath", "mask", "pattern", "symbol", "marker")

function svg_ids!(ids, e)
    id = attribute(e, "id")
    isnothing(id) || (ids[id] = e)
    foreach(c -> svg_ids!(ids, c), child_elements(e))
    return ids
end

# All shapes with their transform from their coordinates to those of the drawing, their name, and whether they or an
# ancestor are hidden by `display: none`, in document order. A shape is named by its own label, else by `label` (that of its closest labeled group or Inkscape layer, or
# the layer's id), else by its own id, else by `id` (that of its closest group with an id); Inkscape writes an id on
# every element. `ids` maps ids to elements, for the copies drawn by `use`, and `used` holds the targets of the `use`
# elements being expanded.
function svg_shapes!(shapes, e, m, label, id, ids, used, hidden)
    name(e) in SVG_NOT_RENDERED && return shapes
    hidden = hidden || attribute(e, "display") == "none" || occursin(r"display\s*:\s*none", something(attribute(e, "style"), ""))
    m = svg_compose(m, svg_transform(attribute(e, "transform")))
    ownlabel = something(attribute(e, "label"), attribute(e, "data-name"), Some(nothing))
    label = something(ownlabel, attribute(e, "groupmode") == "layer" ? attribute(e, "id") : nothing, Some(label))
    name(e) in SVG_SHAPES && push!(shapes, (; element = e, transform = m, name = something(label, attribute(e, "id"), id, name(e)), hidden))
    id = something(attribute(e, "id"), Some(id))
    if name(e) == "use"   # a copy of the element it refers to, e.g. an Inkscape clone, placed like a group
        href = something(attribute(e, "href"), "")   # also matches xlink:href
        target = get(ids, lstrip(href, '#'), nothing)
        isnothing(target) && throw(ArgumentError("a <use> element refers to \"$href\", which is not in the file"))
        target in used && throw(ArgumentError("the <use> of \"$href\" refers to itself, directly or through other <use> elements"))
        m = svg_compose(m, svg_translation(svg_number(e, "x"), svg_number(e, "y")))
        name(target) == "symbol" && !isnothing(attribute(target, "viewBox")) &&
            throw(ArgumentError("a <use> refers to the <symbol> \"$href\", which has a viewBox; symbols with a viewBox are not supported"))
        n = length(shapes)
        foreach(c -> svg_shapes!(shapes, c, m, label, id, ids, (used..., target), hidden), name(target) == "symbol" ? child_elements(target) : (target,))
        isnothing(ownlabel) || (shapes[(n + 1):end] = [(; s..., name = ownlabel) for s in shapes[(n + 1):end]])
        return shapes
    end
    foreach(c -> svg_shapes!(shapes, c, m, label, id, ids, used, hidden), child_elements(e))
    return shapes
end

# path data of a shape
function svg_pathdata(e)
    n = name(e)
    n == "path" && return something(attribute(e, "d"), "")
    if n == "rect"
        x, y, w, h = svg_number.(Ref(e), ("x", "y", "width", "height"))
        # a missing corner radius equals the other one
        rx, ry = svg_number(e, "rx", something(attribute(e, "ry"), "0")), svg_number(e, "ry", something(attribute(e, "rx"), "0"))
        rx, ry = min(rx, w / 2), min(ry, h / 2)
        a = "A $rx,$ry 0 0 1"
        return "M $(x + rx),$y H $(x + w - rx) $a $(x + w),$(y + ry) V $(y + h - ry) $a $(x + w - rx),$(y + h) " *
            "H $(x + rx) $a $x,$(y + h - ry) V $(y + ry) $a $(x + rx),$y Z"
    elseif n == "circle" || n == "ellipse"
        cx, cy = svg_number(e, "cx"), svg_number(e, "cy")
        rx, ry = n == "circle" ? (svg_number(e, "r"), svg_number(e, "r")) : (svg_number(e, "rx"), svg_number(e, "ry"))
        return "M $(cx - rx),$cy A $rx,$ry 0 1 0 $(cx + rx),$cy A $rx,$ry 0 1 0 $(cx - rx),$cy Z"
    elseif n == "line"
        return "M $(svg_number(e, "x1")),$(svg_number(e, "y1")) L $(svg_number(e, "x2")),$(svg_number(e, "y2"))"
    else   # polyline, polygon
        return "M " * something(attribute(e, "points"), "") * (n == "polygon" ? " Z" : "")
    end
end

# Cubic Bézier segments of the elliptical arc from p0 to p with radii rx, ry and x axis rotated by φ degrees, following
# the endpoint to center conversion of the SVG implementation notes
function arc_beziers(p0, rx, ry, φ, large, sweep, p)
    p0 == p && return NTuple{4, NTuple{2, Float64}}[]
    (rx == 0 || ry == 0) && return [(p0, p0, p, p)]
    rx, ry, c, s = abs(rx), abs(ry), cosd(φ), sind(φ)
    dx, dy = (p0 .- p) ./ 2
    x1, y1 = c * dx + s * dy, -s * dx + c * dy
    Λ = x1^2 / rx^2 + y1^2 / ry^2
    Λ > 1 && (rx *= sqrt(Λ); ry *= sqrt(Λ))   # radii too small to reach p are scaled up
    k = sqrt(max(0.0, (rx^2 * ry^2 - rx^2 * y1^2 - ry^2 * x1^2) / (rx^2 * y1^2 + ry^2 * x1^2)))
    large == sweep && (k = -k)
    cx1, cy1 = k * rx * y1 / ry, -k * ry * x1 / rx
    center = (c * cx1 - s * cy1, s * cx1 + c * cy1) .+ (p0 .+ p) ./ 2
    u, v = ((x1 - cx1) / rx, (y1 - cy1) / ry), ((-x1 - cx1) / rx, (-y1 - cy1) / ry)
    θ0, Δθ = atan(u[2], u[1]), atan(u[1] * v[2] - u[2] * v[1], u[1] * v[1] + u[2] * v[2])
    !sweep && Δθ > 0 && (Δθ -= 2π)
    sweep && Δθ < 0 && (Δθ += 2π)
    point(θ) = center .+ (rx * c * cos(θ) - ry * s * sin(θ), rx * s * cos(θ) + ry * c * sin(θ))
    tangent(θ) = (-rx * c * sin(θ) - ry * s * cos(θ), -rx * s * sin(θ) + ry * c * cos(θ))
    n = ceil(Int, abs(Δθ) / (π / 2))
    δ = Δθ / n
    α = 4 / 3 * tan(δ / 4)
    return [
        (
            i == 1 ? p0 : point(θ0 + (i - 1) * δ), point(θ0 + (i - 1) * δ) .+ α .* tangent(θ0 + (i - 1) * δ),
            point(θ0 + i * δ) .- α .* tangent(θ0 + i * δ), i == n ? p : point(θ0 + i * δ),
        ) for i in 1:n
    ]
end

# Subpaths of the path data of `shape`: their segments, each a cubic Bézier curve given by 4 points (lines and quadratic
# curves are converted), and whether they are closed
function svg_subpaths(shape)
    d = svg_pathdata(shape.element)
    re = r"[A-Za-z]|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
    stray = filter(c -> !isspace(c) && c != ',', replace(d, re => " "))
    isempty(stray) || throw(ArgumentError("the path data of \"$(shape.name)\" contains \"$(first(stray, 20))\", which is neither a command nor a number"))
    tokens = [m.match for m in eachmatch(re, d)]
    isempty(tokens) || tokens[1] in ("M", "m") ||
        throw(ArgumentError("the path data of \"$(shape.name)\" must start with M or m, got \"$(tokens[1])\""))
    Segment = NTuple{4, NTuple{2, Float64}}
    subpaths = @NamedTuple{segments::Vector{Segment}, closed::Bool}[]
    pos, start = (0.0, 0.0), (0.0, 0.0)
    ctrl = nothing   # (kind, point): the last control point of a preceding cubic or quadratic segment
    i, cmd = 1, ' '
    function take(n)
        args = tokens[i:min(end, i + n - 1)]
        length(args) == n && !any(t -> isletter(t[1]), args) ||
            throw(ArgumentError("the path command $cmd of \"$(shape.name)\" needs $n numbers, got \"$(join(args, ' '))\""))
        i += n
        return parse.(Float64, args)
    end
    line(p, q) = (p, p, q, q)
    while i <= length(tokens)
        if isletter(tokens[i][1])
            cmd = tokens[i][1]
            i += 1
        end
        rel = islowercase(cmd)
        off = rel ? pos : (0.0, 0.0)
        C = uppercase(cmd)
        next = nothing
        if C != 'M' && C != 'Z' && subpaths[end].closed   # drawing on after Z starts a new subpath
            push!(subpaths, (; segments = Segment[], closed = false))
        end
        if C == 'M'
            x, y = take(2)
            pos = start = (off[1] + x, off[2] + y)
            push!(subpaths, (; segments = Segment[], closed = false))
            cmd = rel ? 'l' : 'L'      # further pairs are line-tos
        elseif C == 'L' || C == 'H' || C == 'V'
            p = C == 'L' ? (off[1], off[2]) .+ Tuple(take(2)) : C == 'H' ? ((rel ? pos[1] : 0.0) + take(1)[1], pos[2]) : (pos[1], (rel ? pos[2] : 0.0) + take(1)[1])
            push!(subpaths[end].segments, line(pos, p))
            pos = p
        elseif C == 'C' || C == 'S'
            v = take(C == 'C' ? 6 : 4)
            c1 = C == 'C' ? (off[1] + v[1], off[2] + v[2]) : !isnothing(ctrl) && ctrl[1] == :cubic ? 2 .* pos .- ctrl[2] : pos
            c2, p = C == 'C' ? ((off[1] + v[3], off[2] + v[4]), (off[1] + v[5], off[2] + v[6])) : ((off[1] + v[1], off[2] + v[2]), (off[1] + v[3], off[2] + v[4]))
            push!(subpaths[end].segments, (pos, c1, c2, p))
            pos, next = p, (:cubic, c2)
        elseif C == 'Q' || C == 'T'
            v = take(C == 'Q' ? 4 : 2)
            c = C == 'Q' ? (off[1] + v[1], off[2] + v[2]) : !isnothing(ctrl) && ctrl[1] == :quadratic ? 2 .* pos .- ctrl[2] : pos
            p = (off[1] + v[end - 1], off[2] + v[end])
            push!(subpaths[end].segments, (pos, pos .+ 2 / 3 .* (c .- pos), p .+ 2 / 3 .* (c .- p), p))   # degree elevation
            pos, next = p, (:quadratic, c)
        elseif C == 'A'
            flags = (get(tokens, i + 3, "0"), get(tokens, i + 4, "0"))
            all(in(("0", "1")), flags) ||
                throw(ArgumentError("the flags of an SVG arc must each be 0 or 1, separated from the next number (as in \"a 5 5 0 0 1 10 0\"); got $(join(flags, ' '))"))
            rx, ry, φ, large, sweep, x, y = take(7)
            p = (off[1] + x, off[2] + y)
            append!(subpaths[end].segments, arc_beziers(pos, rx, ry, φ, large != 0, sweep != 0, p))
            pos = p
        elseif C == 'Z'
            i <= length(tokens) && !isletter(tokens[i][1]) &&
                throw(ArgumentError("the SVG path command Z takes no numbers, got \"$(tokens[i])\" after it"))
            pos == start || push!(subpaths[end].segments, line(pos, start))
            subpaths[end] = (; subpaths[end].segments, closed = true)
            pos = start
        else
            throw(ArgumentError("unsupported SVG path command \"$cmd\""))
        end
        ctrl = next
    end
    return filter(sp -> !isempty(sp.segments), subpaths)
end

# Points of a cubic Bézier segment after the first, except one equal to the last point of `pts`, subdivided until its
# control points are within `tol` of its chord
function svg_flatten!(pts, seg, tol)
    p0, pend = first(seg), last(seg)
    dx, dy = pend .- p0
    L = hypot(dx, dy)
    along(q) = L == 0 ? 0.0 : (dx * (q[1] - p0[1]) + dy * (q[2] - p0[2])) / L^2   # projection onto the chord, 0 to 1 between its ends
    dist(q) = L == 0 ? hypot((q .- p0)...) : abs(dx * (q[2] - p0[2]) - dy * (q[1] - p0[1])) / L
    if all(q -> dist(q) <= tol && -1.0e-9 <= along(q) <= 1 + 1.0e-9, seg[2:(end - 1)])
        pend == pts[end] || push!(pts, pend)
    else   # de Casteljau at t = 1/2
        levels = [collect(seg)]
        while length(levels[end]) > 1
            l = levels[end]
            push!(levels, [(l[k] .+ l[k + 1]) ./ 2 for k in 1:(length(l) - 1)])
        end
        svg_flatten!(pts, Tuple(first.(levels)), tol)
        svg_flatten!(pts, Tuple(reverse(last.(levels))), tol)
    end
    return pts
end

# first and last point of a shape, in drawing coordinates
function svg_endpoints(shape)
    sps = svg_subpaths(shape)
    isempty(sps) && throw(ArgumentError("the reference shape \"$(shape.name)\" is empty"))
    return svg_apply(shape.transform, first(first(sps).segments[1])), svg_apply(shape.transform, last(last(sps).segments[end]))
end

# map from drawing coordinates to model coordinates that takes the points p and q of the drawing to P and Q
function svg_reference_map(p, q, P, Q)
    p[1] != q[1] && p[2] != q[2] && P[1] != Q[1] && P[2] != Q[2] ||
        throw(ArgumentError("the reference must span a nonzero extent in x and y, but it maps the drawing points $p and $q to the model points $P and $Q"))
    return pt -> (P[1] + (pt[1] - p[1]) / (q[1] - p[1]) * (Q[1] - P[1]), P[2] + (pt[2] - p[2]) / (q[2] - p[2]) * (Q[2] - P[2]))
end

function read_outlines_svg(filename, tolerance, bounds, reference)
    isnothing(bounds) || isnothing(reference) || throw(ArgumentError("give either `bounds` or `reference`, not both"))
    doc = parse_file(filename)
    try
        svg = root(doc)
        name(svg) == "svg" || throw(ArgumentError("the root element of $filename is <$(name(svg))>, not <svg>"))
        shapes = @NamedTuple{element::XMLElement, transform::Affine, name::String, hidden::Bool}[]
        ids = svg_ids!(Dict{String, XMLElement}(), svg)
        foreach(c -> svg_shapes!(shapes, c, IDENTITY, nothing, nothing, ids, (), false), child_elements(svg))
        coordref = filter(s -> !isnothing(attribute(s.element, "CoordRef")), shapes)
        notoutline = coordref
        tomodel = if !isnothing(reference)
            refname, (P, Q) = reference
            refs = filter(s -> s.name == refname, shapes)              # a reference may lie in a hidden layer
            length(refs) == 1 || throw(ArgumentError("$filename has $(length(refs)) shapes named \"$refname\"; the reference must be exactly one"))
            notoutline = [notoutline; refs]
            svg_reference_map(svg_endpoints(only(refs))..., P, Q)
        elseif !isnothing(bounds)
            vb = attribute(svg, "viewBox")
            isnothing(vb) && throw(ArgumentError("$filename has no viewBox to map to `bounds`"))
            x0, y0, w, h = svg_floats(vb, "the viewBox of $filename", 4)
            xmin, xmax, ymin, ymax = bounds
            svg_reference_map((x0, y0 + h), (x0 + w, y0), (xmin, ymin), (xmax, ymax))   # SVG's y points down
        elseif !isempty(coordref)
            length(coordref) == 1 || throw(ArgumentError("$filename has $(length(coordref)) reference shapes with a CoordRef attribute; expected at most one"))
            v = svg_floats(attribute(only(coordref).element, "CoordRef"), "the CoordRef of \"$(only(coordref).name)\"", 4)
            svg_reference_map(svg_endpoints(only(coordref))..., (v[1], v[2]), (v[3], v[4]))
        else
            throw(ArgumentError("$filename has no shape with a CoordRef attribute; give `reference = name => ((X0, Y0), (X1, Y1))` or `bounds = (xmin, xmax, ymin, ymax)`"))
        end
        outlines = Outline[]
        for s in shapes
            (s.hidden || any(r -> r.element == s.element, notoutline)) && continue
            for sp in svg_subpaths(s)
                segs = [tomodel.(svg_apply.(Ref(s.transform), seg)) for seg in sp.segments]
                pts = [first(first(segs))]
                foreach(seg -> svg_flatten!(pts, seg, tolerance), segs)
                closed = sp.closed || length(pts) > 2 && pts[end] == pts[1]
                closed && pts[end] == pts[1] && pop!(pts)
                push!(outlines, (; name = s.name, closed, x = first.(pts), y = last.(pts)))
            end
        end
        isempty(outlines) && throw(ArgumentError("$filename contains no outlines, i.e. visible shapes other than the reference"))
        return outlines
    finally
        free(doc)
    end
end
