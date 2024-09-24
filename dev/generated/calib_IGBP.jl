using PML, Ipaper

using Plots
x = range(0, stop=6π, length=1000)
y1 = sin.(x)
y2 = cos.(x)
plot(x, [y1, y2])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
