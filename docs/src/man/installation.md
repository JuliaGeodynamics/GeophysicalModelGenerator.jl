# Installation instructions

GeophysicalModelGenerator.jl is written in the [julia](https://julialang.org) programming language, which is an extremely powerful, modern, scientific computing language. Julia works on all major operating systems, is free, fast, and has a very active user basis (with many useful packages). In case you haven't heard about julia yet, you are not alone. Yet, perhaps a look at [this](https://www.nature.com/articles/d41586-019-02310-3) or [this](https://thenextweb.com/news/watch-out-python-julia-programming-coding-language-coming-for-crown-syndication) article, which explains nicely why it has an enormous potential for computational geosciences as well.

### 1. Install julia
In order to use then package you obviously need to install julia. We recommend downloading and installing binaries from the [julia](https://julialang.org) webpage.


### 2. Install Visual Studio Code
The julia files itself are text files (just like matlab scripts). You may want to edit or modify them at some stage, for which you can use any text editor for that. We prefer to use the freely available [Visual Studio Code](https://code.visualstudio.com) as it has a build-in terminal and is the comes with the (official) julia debugger (install the Julia extension for that).

### 3. Getting started with julia
You start julia on the command line with:
```
kausb$ julia
```
This will start the command-line interface of julia:
```julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.0 (2021-03-24)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> 
```

From the julia prompt, you start the package manager by typing `]`:
```julia
(@v1.6) pkg> 
```
And you return to the command line with a backspace.

Also useful is that julia has a build-in terminal, which you can reach by typing `;` on the command line:
```julia
julia>;
shell> 
```
In the shell, you can use the normal commands like listing the content of a directory, or the current path:
```julia
shell> ls
LICENSE         Manifest.toml   Project.toml    README.md       docs            src             test            tutorial
shell> pwd
/Users/kausb/.julia/dev/GeophysicalModelGenerator
```
As before, return to the main command line (called `REPL`) with a backspace.

If you want to see help information for any julia function, type `?` followed by the command. 
An example for `tan` is:
```julia
help?> tan
search: tan tanh tand atan atanh atand instances transpose transcode contains UnitRange ReentrantLock StepRange StepRangeLen trailing_ones trailing_zeros

  tan(x)

  Compute tangent of x, where x is in radians.

  ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

  tan(A::AbstractMatrix)

  Compute the matrix tangent of a square matrix A.

  If A is symmetric or Hermitian, its eigendecomposition (eigen) is used to compute the tangent. Otherwise, the tangent is determined by calling exp.

  Examples
  ≡≡≡≡≡≡≡≡≡≡

  julia> tan(fill(1.0, (2,2)))
  2×2 Matrix{Float64}:
   -1.09252  -1.09252
   -1.09252  -1.09252
``` 

If you are in a directory that has a julia file (which have the extension `*.jl`), you can open that file with Visual Studio Code:
```julia
shell> code runtests.jl
```
Execute the file with:
```julia
julia> include("runtests")
```
Note that you do not include the `*.jl` extension.


### 4. Install GeophysicalModelGenerator.jl
In order to install GeophysicalModelGenerator.jl, start julia and go to the package manager:
```julia
julia> ]
(@v1.6) pkg> add https://github.com/JuliaGeodynamics/GeophysicalModelGenerator
```
This will automatically install various other packages it relies on (using the correct version).

If you want, you can test if it works on your machine by running the test suite in the package manager:
```julia
julia> ]
(@v1.6) pkg> test GeophysicalModelGenerator
```
Note that we run these tests automatically on Windows, Linux and Mac every time we add a new feature to GeophysicalModelGenerator (using different julia versions). This Continuous Integration (CI) ensures that new features do not break others in the package. The results can be seen [here](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl/actions).

The installation of `GMG` only needs to be done once, and will precompile the package and all other dependencies.

If you, at a later stage, want to upgrade to the latest version of `GMG`, you can type:
```julia
julia> ]
(@v1.6) pkg> update GeophysicalModelGenerator
```

You can load GeophysicalModelGenerator, for example to create cross-sections, with:
```julia
julia> using GeophysicalModelGenerator
```

### 5. Other useful packages
As you will work your way through the tutorials you will see that we often use external packages, for example to load ascii data files into julia. You will find detailed instructions in the respective tutorials. 

If you already want to install some of those, here our favorites. Install them through the package manager:

- [CSV](https://github.com/JuliaData/CSV.jl): Read comma-separated data files into julia.  
- [Plots](https://github.com/JuliaPlots/Plots.jl): Create all kinds of plots in julia (quite an extensive package, but very useful to have). 
- [JLD2](https://github.com/JuliaIO/JLD2.jl): This allows saving julia objects (such as a tomographic model) to a binary file and load it again at a later stage.
- [Geodesy](https://github.com/JuliaGeo/Geodesy.jl): Convert UTM coordinates to latitude/longitude/altitude.
- [NetCDF](https://github.com/JuliaGeo/NetCDF.jl): Read NetCDF files.
- [GMT](https://github.com/GenericMappingTools/GMT.jl): A julia interface to the Generic Mapping Tools (GMT), which is a highly popular package to create (geophysical) maps. Note that installing `GMT.jl` is more complicated than installing the other packages listed above, as you first need to have a working version of `GMT` on your machine (it is not yet installed automatically). Installation instructions for Windows/Linux are on their webpage. On a mac, we made the best experiences by downloading the binaries from their webpage and not using a package manager to install GMT.




