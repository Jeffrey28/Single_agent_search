#!/usr/bin/env julia

using RobotOS
@rosimport std_msgs.msg: Float64MultiArray, Float64Msg
rostypegen()
using std_msgs.msg

using CarMpcUtils
using CarMpc

### MPC parameters
# initialize vehicle
vehicle = Vehicle([0.0, 0.0, 0.0, 0.0])

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

### Subscriber callbacks
callback(msg::Float64MultiArray) = begin
    init_state = msg.data;
end

sub_state_callback(sub_msg::Float64MultiArray, z0::Array{Float64,1}) = begin
  z0[1:nz] = sub_msg.data
end

################
##### Main #####
################

### MPC model parameters updated by subscribers
z0 = vehicle.z
u0 = zeros(nu)

### ROS subroutines
init_node("mpc_julia")
pub_steer = Publisher{Float64MultiArray}("mpc_out_steer", queue_size=1)
pub_acc = Publisher{Float64MultiArray}("mpc_out_acc", queue_size=1)
sub = Subscriber{Float64MultiArray}("mpc_solver_feed", sub_state_callback, (z0,), queue_size=1)

### Main loop
while !is_shutdown()
  ### Reference generation
  URef, ZRef = generateReference(z0, Map, tuning, mpc)

  ### Update and solve MPC problem
  ZOpt, UOpt, solveTime = updateSolveMpcProblem(mpc, ZRef[:,2:N+1], URef, z0, u0)

  ### Publish outputs
  steer_msg = Float64Msg()
  steer_msg.data = UOpt[1,1]
  publish(pub_steer, steer_msg)

  acc_msg = Float64Msg()
  acc_msg.data = UOpt[2,1]
  publish(pub_acc, acc_msg)

  ### Update current input
  u0[1:nu] = UOpt[:,1]
end
