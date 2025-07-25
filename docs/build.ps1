julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path = @__DIR__)); Pkg.instantiate(); Pkg.build()'
julia --project=docs/ docs/make.jl
