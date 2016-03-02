# Test script for controller.jl
push!(LOAD_PATH,pwd())
using CarMpcUtils
using CarMpc

### Simulation parameters
simLength = 200

# Map
map = SlamMap("maps/shared_room_map.mat")

# # initialize vehicle
# vehicle = Vehicle([122.724, 44.3005, 2.615, 5.0], model=kinematicBicycleDiscrete)   # RFS 0128
rob = Robot([1, 1, 0, 0.5], model=unicycleDiscrete)   # RFS 0128, map

# # tuning parameters
# tuning = Tuning(dt = 0.2, dtMPC = 0.2, N = 20,
#                 Q = [0.5, 0.5, 10.0, 0.0], R = [20.0, 2.0],
#                 P = [1000.0, 20.0], vRef = 10.0, dSafe = 5.0,
#                 eYRef = 0.0, TTC = 3.0, eYSafe = 0.5)
tuning = Tuning(dt = 1.0, dtMPC = 1.0, N = 2)

# map
# Map = TrackMap("maps/RFS_2Lanes_Speed_0128.mat")

# ### Initialize MPC problem and solve dummy problem
# mpc = initializeMPCProblem(vehicle, tuning)
mpc = initializeMPCProblem(rob, map, tuning)
nz, nu, N = mpc.nz, mpc.nu, tuning.N


# ################
# ##### Main #####
# ################

# ### MPC model parameters updated by subscribers
z0 = vehicle.z
u0 = zeros(nu)

USim = zeros(nu,simLength)
ZSim = [z0 zeros(nz,simLength)]

### Main loop
for t=1:simLength
  ### Reference generation
  # URef, ZRef = generateReference(z0, Map, tuning, mpc)

  ### decide if search terminates


  ### update prob map using observation
  # get obs from ROS publisher
  obs = [100;100]
  updateMap(rob,obj,map)

  ### Update and solve MPC problem
  ### missing part: need to incorporate map into MPC
  # ZOpt, UOpt, solveTime = updateSolveMpcProblem(mpc, ZRef[:,2:N+1], URef, z0, u0)
  ZOpt, UOpt, solveTime = updateSolveMpcProblem(mpc, z0, u0)

  ### Update current input
  u0[1:nu] = UOpt[:,1]

  ### Update ego state
  # updateEgoState!(vehicle,u0,tuning.dtMPC)
  updateEgoState!(rob,u0,tuning.dtMPC)

  ### Variables for logging
  # ZSim[:,t+1] = vehicle.z
  ZSim[:,t+1] = rob.z
  USim[:,t] = u0
end

# using PyPlot
# figure()
# plot(Map.nodes[1,:]',Map.nodes[2,:]')
# plot(ZSim[1,:]',ZSim[2,:]',marker="o")
# show()