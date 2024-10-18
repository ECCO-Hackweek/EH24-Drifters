
## Julia setup

if !isdefined(Main,:IndividualDisplacements)
    #using Pkg; Pkg.activate("scripts")
    using IndividualDisplacements, CairoMakie, Climatology, NetCDF
    p0=joinpath(dirname(pathof(IndividualDisplacements)),"..","examples")
    f0=joinpath(p0,"worldwide","OCCA_FlowFields.jl")
    include(f0);
end

# Ocean Circulation setup

P,D=OCCA_FlowFields.setup(backward_in_time=true, nmax=5) #parameters, diagnostics
G=D.Γ #grid parameters
rec=OCCA_FlowFields.custom🔴 #recorder
proc=OCCA_FlowFields.custom🔧 #processor
step_forward! =∫! #integrate over time

# Initial Conditions (positions)

nf=1000; lo=(-160.0,-150.0); la=(30.0,40.0); level=2.5;
df=OCCA_FlowFields.initial_positions(G, nf, lo, la, level)

# Individuals setup

I=Individuals(P,df.x,df.y,df.z,df.f,(🔴=rec,🔧=proc, 𝐷=D))

# 10-day integration

T=(0.0,10*86400.0)
step_forward!(I,T)

# visulation

"""
    myplot(I::Individuals)

Plot the initial and final positions as scatter plot in `lon,lat` or `x,y` plane.
"""
function myplot(I::Individuals)
	🔴_by_t = IndividualDisplacements.DataFrames.groupby(I.🔴, :t)
	set_theme!(theme_black())
	fig=Figure(size = (600, 400))
    a = Axis(fig[1, 1],xlabel="longitude",ylabel="latitude")		
    scatter!(a,🔴_by_t[1].lon,🔴_by_t[1].lat,color=:green2,markersize=4,label="initial positions")
    scatter!(a,🔴_by_t[end].lon,🔴_by_t[end].lat,color=:red,markersize=4,label="final positions")
    axislegend(a)
    return fig
end

myplot(I)
