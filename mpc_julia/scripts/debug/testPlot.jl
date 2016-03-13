# test plotting with pyplot events
push!(LOAD_PATH,pwd())
using CarMpcUtils
using Polynomials, Interpolations
using PyPlot
using PyCall
@pyimport matplotlib.patches as patches
@pyimport PlotUtils as plt_utils

println(pwd())
Map = TrackMap(pwd()*"/maps/RFS_2Lanes_Speed_0128.mat")

fig = figure()
ax1 = fig["add_subplot"](111)
ax1["set_title"]("custom picker for line data", picker=true)
line = ax1["plot"](Map.nodes[1,:]',Map.nodes[2,:]', "o", picker=5)
fig["canvas"]["mpl_connect"]("pick_event", plt_utils.onpick1)
show()
