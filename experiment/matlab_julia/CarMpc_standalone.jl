# this file tests the simple case of objective function

using JuMP
using Ipopt
# using Polynomials, Interpolations

import JuMP.JuMPArray, JuMP.Variable

N = 30 # prediction horizon
# robot 
init_pose = [30;30;0]

# map
x_len = 30 # x length
y_len = 30 # y length
xp = 1:(x_len+1) # points in x direction, i.e. horizontal distance between adjacent points are 0.2
yp = 1:(y_len+1) # points in y direction
ter_pos = [5;5]
ter_r = 3

zLB = [xp[1],yp[1]]
zUB = [xp[end],yp[end]]
uLB = [0,-pi/4]
uUB = [1.5,pi/4]
dt = 1
# formulate MPC problem
model = Model(solver=IpoptSolver(max_iter=1000, max_cpu_time=100.0, print_timing_statistics="yes")) #max_cpu_time=dtMPC, hessian_approximation ="limited-memory"

  # ### Dimensions
  # nz = 3
  # nu = 2

  # ### Dummies
  # # (Parameters for optimization problem that can be modified online)
  # z0 = init_pose
  # u0 = zeros(nu)

  # ### Variables
  # @defVar(model, Z[1:nz,1:N])
  # # @defVar(model, uLB[i] <= U[i=1:nu,1:N] <= uUB[i])
  # @defVar(model, U[1:nu,1:N])

  # ### Objective function

  # # @defNLExpr(obj[j=1:N], abs((Z[1,j]-ter_pos[1])^2+(Z[2,j]-ter_pos[2])^2-ter_r^2))
  # # @setNLObjective(model, Min, sum{obj,i=1:N})

  # # @setNLObjective(model, Min, sum{abs((Z[1,j]-ter_pos[1])^2+(Z[2,j]-ter_pos[2])^2-ter_r^2),j=1:N})
  # @setNLObjective(model, Min, abs((Z[1,N]-ter_pos[1])^2+(Z[2,N]-ter_pos[2])^2-ter_r^2))

  ### Constraints
  # for k=1:N
  #   # Dynamics
  #   if k==1
  #     @addNLConstraints(model, begin
  #                         Z[1,k] == z0[1] + dt*U[1,k]*cos(z0[3]+dt*U[2,k])
  #                         Z[2,k] == z0[2] + dt*U[1,k]*sin(z0[3]+dt*U[2,k])
  #                         Z[3,k] == z0[3] + dt*U[2,k]                          
  #                       end)
  #   else
  #     @addNLConstraints(model, begin
  #                         Z[1,k] == Z[1,k-1] + dt*U[1,k]*cos(Z[3,k-1]+dt*U[2,k])
  #                         Z[2,k] == Z[2,k-1] + dt*U[1,k]*sin(Z[3,k-1]+dt*U[2,k])
  #                         Z[3,k] == Z[3,k-1] + dt*U[2,k]  
  #                       end)
  #   end
  #   # state constraints
  #   @addConstraint(model, zLB[1] <=Z[1,k] <= zUB[1])
  #   @addConstraint(model, zLB[2] <=Z[2,k] <= zUB[2])
  #   @addConstraint(model, uLB[1] <=U[1,k] <= uUB[1])
  #   @addConstraint(model, uLB[2] <=U[2,k] <= uUB[2])
  # end
  # @addNLConstraint(model,(Z[1,N]-ter_pos[1])^2+(Z[2,N]-ter_pos[2])^2 <= ter_r^2)


  maxV = 1.5;

  ### Dimensions
  nz = 2
  nu = 2

  ### Dummies
  # (Parameters for optimization problem that can be modified online)
  z0 = init_pose[1:2]
  u0 = zeros(nu)

  ### Variables
  @defVar(model, Z[1:nz,1:N])
  @defVar(model, U[1:nu,1:N])
  @defVar(model, dummy_t[1:N]>=0)

  ### initial solution
  ux_init = (ter_pos[1]-z0[1])/N;
  uy_init = (ter_pos[2]-z0[2])/N;

  if (ux_init^2+uy_init^2) <= maxV^2
    U_init = [ux_init*ones(1,N);uy_init*ones(1,N)]
  else
    angle = atan2(ter_pos[2]-z0[2],ter_pos[1]-z0[1])
    U_init = [maxV*cos(angle)*ones(1,N);maxV*sin(angle)*ones(1,N)]
  end

  Z_init = zeros(2,N)
  for k=1:N
    if k == 1
      Z_init[1,k] = z0[1] + dt*U_init[1,k]
      Z_init[2,k] = z0[2] + dt*U_init[2,k]
      # println(Z_init)
    else
      Z_init[1,k] = Z_init[1,k-1] + dt*U_init[1,k]
      Z_init[2,k] = Z_init[2,k-1] + dt*U_init[2,k]
      # println(Z_init)
    end
  end

  # for k = 1:N
  #   setValue(U[k],U_init[k])
  #   setValue(Z[k],Z_init[k])
  # end

  ### Objective function

  # @defNLExpr(obj[j=1:N], abs((Z[1,j]-ter_pos[1])^2+(Z[2,j]-ter_pos[2])^2-ter_r^2))
  # @setNLObjective(model, Min, sum{obj,i=1:N})

  # @setNLObjective(model, Min, sum{abs((Z[1,j]-ter_pos[1])^2+(Z[2,j]-ter_pos[2])^2-ter_r^2),j=1:N})
  @setNLObjective(model, Min, sum{dummy_t[i],i=1:N}) # (Z[1,N]-ter_pos[1])^2+(Z[2,N]-ter_pos[2])^2-ter_r^2

  ### Constraints
  for k=1:N
    # Dynamics
    if k==1
      @addNLConstraints(model, begin
                          Z[1,k] == z0[1] + dt*U[1,k]
                          Z[2,k] == z0[2] + dt*U[2,k]
                        end)
    else
      @addNLConstraints(model, begin
                          Z[1,k] == Z[1,k-1] + dt*U[1,k]
                          Z[2,k] == Z[2,k-1] + dt*U[2,k] 
                        end)
    end
    # state and input constraints
    @addNLConstraint(model, U[1,k]^2+U[2,k]^2<=maxV^2)
    @addConstraint(model, zLB[1] <=Z[1,k] <= zUB[1])
    @addConstraint(model, zLB[2] <=Z[2,k] <= zUB[2])   

    # epigraph variable
    @addNLConstraint(model,dummy_t[k]>=(Z[1,k]-ter_pos[1])^2+(Z[2,k]-ter_pos[2])^2-ter_r^2)
    @addNLConstraint(model,dummy_t[k]>=-((Z[1,k]-ter_pos[1])^2+(Z[2,k]-ter_pos[2])^2-ter_r^2))
  end
  # @addNLConstraint(model,t>=(Z[1,N]-ter_pos[1])^2+(Z[2,N]-ter_pos[2])^2-ter_r^2)
  # @addNLConstraint(model,t>=-((Z[1,N]-ter_pos[1])^2+(Z[2,N]-ter_pos[2])^2-ter_r^2))
  
  @addNLConstraint(model,(Z[1,N]-ter_pos[1])^2+(Z[2,N]-ter_pos[2])^2 <= ter_r^2)

### Solve dummy problem
solve(model)