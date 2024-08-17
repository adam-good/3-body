include("3-body.jl")
using .ThreeBody
using Distributions
using GLMakie

const NBODIES = 3

low_bound = -100.0; high_bound = 100.0
mass = 3e10
prob = Body3(
    Body(rand(Uniform(low_bound, high_bound), 3)..., mass),
    Body(rand(Uniform(low_bound, high_bound), 3)..., mass),
    Body(rand(Uniform(low_bound, high_bound), 3)..., mass),
)

set_theme!(theme_black())
points = [
    Observable(Point3f[]),
    Observable(Point3f[]),
    Observable(Point3f[])
]
colors = Observable(Int[])

fig = Figure()
ax = Axis3(
    fig[1,1],
    viewmode=:fit, protrusions = (0,0,0,0),
    limits = (low_bound, high_bound, low_bound, high_bound, low_bound, high_bound)
)
lines = [
    lines!(ax, points[1], color=colors, colormap=:Wistia, transparency=true),
    lines!(ax, points[2], color=colors, colormap=:brg, transparency=true),
    lines!(ax, points[3], color=colors, colormap=:cool, transparency=true)
]

max_plot_len = 5e5
nframes = 240
iter_per_frame=2e4
record(fig, "3-body.mp4", 1:nframes) do frame
    for i = 1:iter_per_frame
        step!(prob)
        
        for (p,b) in zip(points, [prob.b1, prob.b2, prob.b3])
            push!(p[], Point3f(unpack(b.r)))
        end
        push!(colors[], frame)
        if i % iter_per_frame == 0
            @info """Update $(frame * iter_per_frame + i):
                $(prob.b1.ṙ)
                $(prob.b2.ṙ)
                $(prob.b3.ṙ)
                """
        end
    end

    l = length(points[1][])
    if l > max_plot_len
        diff::Int = l - max_plot_len
        points[1][] = points[1][][diff:l]
        points[2][] = points[2][][diff:l]
        points[3][] = points[3][][diff:l]
        colors[] = colors[][diff:l]
    end

    ax.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / nframes)
    notify(colors)
    for p in points
        notify(p)
    end
    for l in lines
        l.colorrange = (0,frame)
    end
end