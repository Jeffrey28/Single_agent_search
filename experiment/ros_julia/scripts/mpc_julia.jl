#!/usr/bin/env julia

using RobotOS
# @rosimport std_msgs.msg: StringMultiArray
@rosimport std_msgs.msg: Float64MultiArray
# gentypes()
rostypegen()
using std_msgs.msg

using JuMP
using Ipopt


N = 20

lead_position = 0;
lead_velocity = 0;

duration = 0.01;

acceleration_ref = 0;
acceleration_min = -4;
acceleration_max = 1.5;
velocity_max = 6;
v_ref = 5 * ones(N) #make this parameter

u_ref = acceleration_ref * ones(N);

mpc = Model(solver=IpoptSolver(print_level=0))
@defVar(mpc, acceleration_min <= u[1:N] <= acceleration_max)
@defVar(mpc, p[1:N])
@defVar(mpc, -velocity_max <= v[1:N] <= velocity_max)


p0 = [0.0]
v0 = [0.0]

# println("@@@@@@@@@@@@ v0: %f =============", v0);


@setObjective(mpc, Min, sum{(u[i] - u_ref[i])^(2) + 100*(v[i] - v_ref[i])^(2), i=1:N})

@addNLConstraint(mpc, p[1] == p0[1] + duration * v0[1])
@addNLConstraint(mpc, v[1] == v0[1] + duration * u[1])

for i = 1:N-1
    @addConstraint(mpc, p[i+1] == p[i] + duration * v[i])
    @addConstraint(mpc, v[i+1] == v[i] + duration * u[i+1])
    @addConstraint(mpc, v[i] <= velocity_max)
    #@addConstraint(mpc, p[i] <= lead_position + lead_velocity * i * duration)
end


callback(msg::Float64MultiArray, pub_obj::Publisher{Float64MultiArray}) = begin
    # println("Got Data!")
    # println(msg.data)

    # N = 20

    init_state = msg.data;

    # # lead_position = init_state[1];
    # # lead_velocity = init_state[2];
    # # p0[1] = init_state[3]
    # # v0[1] = init_state[4]

    v0[1] = msg.data[4]

    println("@@@@@@@@@@@@ v0: %f =============", v0[1])

    # lead_position = 0;
    # lead_velocity = 0;

    # duration = 0.01;

    # acceleration_ref = 0;
    # acceleration_min = -4;
    # acceleration_max = 1.5;
    # velocity_max = 15;
    # v_ref = 5 * ones(N); #make this parameter

    # u_ref = acceleration_ref * ones(N);

    #u_ref = linspace(10, 20, N);

    # @setObjective(mpc, Min, sum{(u[i] - u_ref[i])^(2) + 100*(v[i] - v_ref[i])^(2), i=1:N})

    # @addConstraint(mpc, p[1] == p0 + duration * v0)
    # @addConstraint(mpc, v[1] == v0 + duration * u[1])

    # for i = 1:N-1
    #     @addConstraint(mpc, p[i+1] == p[i] + duration * v[i])
    #     @addConstraint(mpc, v[i+1] == v[i] + duration * u[i+1])
    #     #@addConstraint(mpc, p[i] <= lead_position + lead_velocity * i * duration)
    # end

    #@addConstraint(mpc, p[N] <= lead_position + lead_velocity * N * duration)

    println("Finished MPC")

    #Make this in call-back
    solve(mpc)
    res_msg = Float64MultiArray()
    res_msg.data = [getValue(u[1])]
    publish(pub_obj, res_msg)
end

println("Starting MPC Solver!!!")

init_node("mpc_julia")
pub = Publisher{Float64MultiArray}("mpc_result", queue_size=10)
sub = Subscriber{Float64MultiArray}("mpc_solver_feed", callback, (pub,), queue_size=10)

spin()
