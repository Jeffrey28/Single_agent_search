# Test script for controller.jl
push!(LOAD_PATH,pwd())
using CarMpcUtils
using CarMpc

### Simulation parameters
simLength = 200

# initialize vehicle
vehicle = Vehicle([122.724, 44.3005, 2.615, 5.0], model=kinematicBicycleDiscrete)   # RFS 0128

# tuning parameters
tuning = Tuning(dt = 0.2, dtMPC = 0.2, N = 20,
                Q = [0.5, 0.5, 10.0, 0.0], R = [20.0, 2.0],
                P = [1000.0, 20.0], vRef = 10.0, dSafe = 5.0,
                eYRef = 0.0, TTC = 3.0, eYSafe = 0.5)

# map
Map = TrackMap("maps/RFS_2Lanes_Speed_0128.mat")

### Initialize MPC problem and solve dummy problem
mpc = initializeMPCProblem(vehicle, tuning)
nz, nu, N = mpc.nz, mpc.nu, tuning.N

################
##### Main #####
################

### MPC model parameters updated by subscribers
z0 = vehicle.z
u0 = zeros(nu)

USim = zeros(nu,simLength)
ZSim = [z0 zeros(nz,simLength)]

### Main loop
for t=1:simLength
  ### Reference generation
  URef, ZRef = generateReference(z0, Map, tuning, mpc)

  ### Update and solve MPC problem
  ZOpt, UOpt, solveTime = updateSolveMpcProblem(mpc, ZRef[:,2:N+1], URef, z0, u0)

  ### Update current input
  u0[1:nu] = UOpt[:,1]

  ### Update ego state
  updateEgoState!(vehicle,u0,tuning.dtMPC)

  ### Variables for logging
  ZSim[:,t+1] = vehicle.z
  USim[:,t] = u0
end

using PyPlot
figure()
plot(Map.nodes[1,:]',Map.nodes[2,:]')
plot(ZSim[1,:]',ZSim[2,:]',marker="o")
show()
