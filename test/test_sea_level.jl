using Test

x = SeaLevel(:Spratt_800ka)

@test length(x)     === length(x.elevation)
@test size(x)       === size(x.elevation)
@test eachindex(x)  === eachindex(x.elevation)
@test axes(x)       === axes(x.elevation)
@test curve_name(x) === :Spratt_800ka
@test x[1]          == x.elevation[1]
@test x[1,1]        == (x.elevation[1], x.age[1])
@test x[1,799]      == (x.elevation[1], x.age[799])
@test x[(1,)]       == (x.elevation[1], x.age[1])
@test x[(799,)]     == (x.elevation[799], x.age[799])

x_rev = SeaLevel(:Spratt_800ka; flip_elevation = true, flip_age = true)
@test x.elevation == reverse(x_rev.elevation)
@test x.age       == reverse(x_rev.age)