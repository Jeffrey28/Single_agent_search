using JuMP
using Ipopt


# callback() = begin
    N = 20

#    lead_position = init_state[1]
#    lead_velocity = init_state[2]
    p0 = 0.3
    v0 = 0

    lead_position = 0
    lead_velocity = 0

    duration = 1.5

    acceleration_ref = 0;
    acceleration_min = -4;
    acceleration_max = 1.5;
    velocity_max = 15;
    v_ref = 1.4 * ones(N); #make this parameter

    u_ref = acceleration_ref * ones(N);
    #u_ref = linspace(10, 20, N);

    mpc = Model(solver=IpoptSolver())
    @defVar(mpc, acceleration_min <= u[1:N] <= acceleration_max)
    @defVar(mpc, p[1:N])
    # @defVar(mpc, v[1:N])
    @defVar(mpc, -velocity_max <= v[1:N] <= velocity_max)
    @setObjective(mpc, Min, sum{(u[i] - u_ref[i])^(2) + (v[i] - v_ref[i])^(2), i=1:N})

    @addConstraint(mpc, p[1] == p0 + duration * v0)
    @addConstraint(mpc, v[1] == v0 + duration * u[1])

    for i = 1:N-1
        @addConstraint(mpc, p[i+1] == p[i] + duration * v[i])
        @addConstraint(mpc, v[i+1] == v[i] + duration * u[i+1])
        #@addConstraint(mpc, p[i] <= lead_position + lead_velocity * i * duration)
    end

    #@addConstraint(mpc, p[N] <= lead_position + lead_velocity * N * duration)

    println("Finished MPC")
    println(getValue(u))
# callback()
