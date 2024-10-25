#Let's start a temporary environment for this notebook, and add julia packages that we will use
if !isdefined(Main,:IndividualDisplacements)
    using Pkg; Pkg.activate(temp=true)
    Pkg.add.(["IndividualDisplacements", "GLMakie", "Climatology", "NetCDF", "MeshArrays", "GeoJSON", "DataDeps"])
 
    using IndividualDisplacements, GLMakie, Climatology, NetCDF, MeshArrays, GeoJSON, DataDeps
	p0=joinpath(dirname(pathof(IndividualDisplacements)),"..","examples")
	f0=joinpath(p0,"worldwide","OCCA_FlowFields.jl")
	include(f0);
end
γ=GridSpec("PeriodicChannel",MeshArrays.GRID_LL360) # using MeshArrays to get Grid information
Γ=GridLoad(γ;option="full")

# Ocean Circulation setup

P,D=OCCA_FlowFields.setup(backward_in_time=false) #backward_in_time = change true/false depending on if you want to fun forward or backward, namx=Inf by default, nmax=5 for small simulations
G=D.Γ #grid parameters
rec=OCCA_FlowFields.custom🔴 #recorder
proc=OCCA_FlowFields.custom🔧 #processor
step_forward! =∫! #integrate over time

# Southern Atlantic
nf=10*1000; lo=(-40.0,30.0); la=(-6.0,-6.0); level=2.5
#basins=demo.ocean_basins()
#AtlExt=demo.extended_basin(basins,:Atl)
#IndExt=demo.extended_basin(basins,:Ind)
#mymask=AtlExt .* (Γ.XC .> -40.0) .* (Γ.XC .<=30.0) .* (Γ.YC .< -6.0) .* ( Γ.hFacC[:,1] .> 0.0 )

df=OCCA_FlowFields.initial_positions(G, nf, lo, la, level)
I=Individuals(P,df.x,df.y,df.z,df.f,(🔴=rec,🔧=proc, 𝐷=D))
#T=(0.0,10*10*86400.0) #use this for only intial and final positon of particles
#step_forward!(I,T)

for i = 1:50 #use this to get multiple positons of particles
    T=(0.0,i*5*86400.0)
    step_forward!(I,T)
end

# Now we have two ways to add land background. Please try both ways and choose the one you like.
# I recommend the second way.
lon180(x)=Float64(x>180.0 ? x-360.0 : x)  # define a function to convert 0-360 to -180-180 
lon_p = (-70,40) # set longitude range you wanna plot, note indices are in the range of -180-180
lat_p = (-30,30) # set latitude range you wanna plot

# visualisation 1

fil=demo.download_polygons("countries.geojson") # using MeshArrays to get countries' map
pol=MeshArrays.read_polygons(fil)

"""
    myplot1(I::Individuals)

Plot the initial and final positions as scatter plot in `lon,lat` or `x,y` plane.
"""
function myplot1(I::Individuals)
	🔴_by_t = IndividualDisplacements.DataFrames.groupby(I.🔴, :t)
	set_theme!(theme_black())
    fig=Figure(size = (600, 400))
    a = Axis(fig[1, 1],xlabel="longitude",ylabel="latitude")		
    scatter!(a,lon180.(🔴_by_t[1].lon),🔴_by_t[1].lat,color=:green2,markersize=4,label="initial positions") # use lon180 to convert longitude range
    scatter!(a,lon180.(🔴_by_t[end].lon),🔴_by_t[end].lat,color=:red,markersize=4,label="final positions") 
    [lines!(a,l1,color = :white, linewidth = 0.5) for l1 in pol] # countries are -180-180

    xlims!(a,lon_p)
    ylims!(a,lat_p)
    axislegend(a)
    return fig
end

fig = myplot1(I)
file = joinpath(tempdir(),"mask_test1.png")
save(file,fig)

# visualisation 2

lndid = findall(Γ.hFacC[1,1].==0); # find indices of all positions of land

"""
    myplot2(I::Individuals)

Plot the initial and final positions as scatter plot in `lon,lat` or `x,y` plane.
"""
function myplot2(I::Individuals)
	🔴_by_t = IndividualDisplacements.DataFrames.groupby(I.🔴, :t)
	set_theme!(theme_black())
    fig=Figure(size = (600, 400))
    a = Axis(fig[1, 1],xlabel="longitude",ylabel="latitude")	
    
    scatter!(a,lon180.(🔴_by_t[1].lon),🔴_by_t[1].lat,color=:green2,markersize=4,label="initial positions") # use lon180 to convert longitude range
    sca = scatter!(a,lon180.(🔴_by_t[end].lon),🔴_by_t[end].lat,color=🔴_by_t[end].z,markersize=4, colormap = :plasma)
    Colorbar(fig[1, 2], sca; label = "Depth (m)")
#   scatter!(a,lon180.(🔴_by_t[end].lon),🔴_by_t[end].lat,color=:red,markersize=4,label="final positions") 
    scatter!(a,lon180.(Γ.XC[1,1][lndid]),Γ.YC[1,1][lndid],color=:white,markersize=10) # scatter the land background

    xlims!(a,lon_p)
    ylims!(a,lat_p)
    axislegend(a)
    return fig
end

fig2=myplot2(I)

# plot region we focus on using heatmap
function plot_masks(Γ)
    msks1=Γ.hFacC[:,1]*(Γ.XC.>0.0)*(Γ.XC.<40.0)*(Γ.YC.>-30.0)*(Γ.YC.<30.0)
    msks2=Γ.hFacC[:,1]*(Γ.XC.>-70.0+360)*(Γ.XC.<0.0+360)*(Γ.YC.>-30.0)*(Γ.YC.<30.0)
    msks = msks1+msks2
    cmap = [:white, :transparent]
    #println(lndid)

    fig=Figure(size = (600, 400))
    a = Axis(fig[1, 1],xlabel="longitude",ylabel="latitude")

    XC=circshift(Γ.XC[1,1],(-180,0)); XC[XC.>180].-=360
    YC=circshift(Γ.YC[1,1],(-180,0))
    M=circshift(msks[1,1],(-180,0))
    heatmap!(XC[:,1],YC[1,:],M)

    xlims!(a,lon_p)
    ylims!(a,lat_p)

    fig
end

fig_masks=plot_masks(Γ)

# visualisation 3

"""
    myplot3(I::Individuals)

Plot all positions over time as scatter plot in `lon,lat` or `x,y` plane incrementally in movie.
"""
function myplot3(I::Individuals)
🔴_by_t = IndividualDisplacements.DataFrames.groupby(I.🔴, :t)

set_theme!(theme_black())

time=Observable(10)

lon = @lift(vcat([lon180.(🔴_by_t[$time-t+1].lon) for t in 1:3]...))#changing to convert lon range in a way that works with observable
lat = @lift(vcat([🔴_by_t[$time-t+1].lat for t in 1:3]...)) #can change 3 for more particles
z = @lift(vcat([🔴_by_t[$time-t+1].z for t in 1:3]...))

f=Figure(size = (600, 400))
a = Axis(f[1, 1],xlabel="longitude",ylabel="latitude")	

scatter!(a,lon,lat,color=:green2,markersize=2) # use lon180 to convert longitude range
sca = scatter!(a,lon,lat,color=z,markersize=2, colormap = :plasma)
Colorbar(fig[1, 2], sca; label = "Depth (m)")
scatter!(a,lon180.(🔴_by_t[1].lon),🔴_by_t[1].lat,color=:green2,markersize=4,label="initial positions") #scatter initial position
scatter!(a,lon180.(Γ.XC[1,1][lndid]),Γ.YC[1,1][lndid],color=:white,markersize=10) # scatter the land background

xlims!(a,lon_p)
ylims!(a,lat_p)
axislegend(a)

framerate = 20
timestamps = 3:51 #change for number of timesteps

file= joinpath(tempdir(),"time_animation2.mp4")
record(f, file, timestamps;
        framerate = framerate) do t
    time[] = t
end

file
end

file3=myplot3(I::Individuals)

