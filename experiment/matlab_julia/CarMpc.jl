module CarMpc

using JuMP
using Ipopt
using CarMpcUtils
using Polynomials, Interpolations

import JuMP.JuMPArray, JuMP.Variable

export MpcProblem, initializeMPCProblem, updateSolveMpcProblem
# export localizeVehicle, generateReference

type MpcProblem
  # Vehicle
  rob::Robot
  map::SlamMap
  # Constraints
  uLB::Array{Float64,1}
  uUB::Array{Float64,1}
  duLB::Array{Float64,1}
  duUB::Array{Float64,1}
  # Parameters
  nz::Int64
  nu::Int64
  tuning::Tuning
  # JuMP model variables
  model::Model
  Z::Array{Variable,2}
  U::Array{Variable,2}
  # JuMP model parameters
  z0::Array{Float64,1}
  u0::Array{Float64,1}
  ZRef::Array{Float64,2}
  URef::Array{Float64,2}
end

### Set up MPC optimization problem
function initializeMPCProblem(rob::Robot, map:: SlamMap, tuning::Tuning)
  ### Parameters
  N, dt, dtMPC = tuning.N, tuning.dt, tuning.dtMPC
  # P, Q, R, S = tuning.P, tuning.Q, tuning.R, tuning.S
  uLB = [rob.angMin, rob.vMin]
  uUB = [rob.angMax, rob.vMax]
  # uLB = [rob.deltaMin, rob.axMin]
  # uUB = [rob.deltaMax, rob.axMax]
  # duLB = [-rob.deltaRateLim, -rob.axRateLim]
  # duUB = [rob.deltaRateLim, rob.axRateLim]
  # l_r = rob.l_r

  prob_map = map.prob_map
  unocc_cor = map.unocc_cor
  ### Optimization problem
#   model = Model(solver=IpoptSolver(linear_solver="ma27",
#                                      max_iter=100,
#                                      max_cpu_time=dtMPC))
  model = Model(solver=IpoptSolver(max_iter=100, max_cpu_time=100.0, print_timing_statistics="yes")) #max_cpu_time=dtMPC, hessian_approximation ="limited-memory"

  ### Dimensions
  nz = 4
  nu = 2

  ### Dummies
  # (Parameters for optimization problem that can be modified online)
  # ZRef = zeros(nz,N)
  # URef = zeros(nu,N)
  z0 = zeros(nz)
  u0 = zeros(nu)

  ### Variables
  @defVar(model, Z[1:nz,1:N])
  @defVar(model, uLB[i] <= U[i=1:nu,1:N] <= uUB[i])
  # @defVar(model, obj_expr) # expression for objective function
  # @defVar(model, tmp_prob_ind) # POD expression for an individual position
  # @defNLExpr(pond_prod[1:N]) # expression for POND of each position at each planning time
  # @defNLExpr(pond_sum[1:size(map.unocc_cor,2)]) # expression for POND of each position
  # @defNLExpr(pond_sum[1:size(map.unocc_cor,2),1:N]) # expression for POND of each position

  ### Objective function

  left_bd = rob.left_bd # this will be calculated from robot's real time positions
  bottom_bd = rob.bottom_bd # this will be calculated from robot's real time positions
  right_bd = rob.right_bd # this will be calculated from robot's real time positions

  n_pt = 100 #size(unocc_cor,2)
  # for j = 1:n_pt
  #   ### missing part: pos_tar and left_bd should be in same coordinate
  #   tmp_cor = map.unocc_cor[:,j]
  #   pos_tar = tmp_cor

  #   # for k = 1:N      
  #     # if k == 1
  #     #   # tmp_prob_ind = sensorModel([-1;-1], rob, map, tmp_cor, Z[:,k])
  #     #   tmp_prob_ind = (rob.eps_1+rob.eps_2-1)/((1+exp(-(pos_tar[1]-(left_bd[1]+Z[1,k]))))*(1+exp(-(pos_tar[2]-(bottom_bd[2]+Z[2,k]))))*(1+exp(-(pos_tar[1]-(right_bd[1]+Z[1,k])))))+(1-rob.eps_2)
  #     # else
  #     #   tmp_prob_ind = tmp_prob_ind*((rob.eps_1+rob.eps_2-1)/((1+exp(-(pos_tar[1]-(left_bd[1]+Z[1,k]))))*(1+exp(-(pos_tar[2]-(bottom_bd[2]+Z[2,k]))))*(1+exp(-(pos_tar[1]-(right_bd[1]+Z[1,k])))))+(1-rob.eps_2))
  #     # end
  #   # end

  #   # @addNLConstraint(model, prod_constr[k=1:N], pond_prod[k] == (rob.eps_1+rob.eps_2-1)/((1+exp(-(pos_tar[1]-(left_bd[1]+Z[1,k]))))*(1+exp(-(pos_tar[2]-(bottom_bd[2]+Z[2,k]))))*(1+exp(-(pos_tar[1]-(right_bd[1]+Z[1,k])))))+(1-rob.eps_2))
  #   @defNLExpr(pond_prod[k=1:N], (rob.eps_1+rob.eps_2-1)/((1+exp(-(pos_tar[1]-(left_bd[1]+Z[1,k]))))*(1+exp(-(pos_tar[2]-(bottom_bd[2]+Z[2,k]))))*(1+exp(-(pos_tar[1]-(right_bd[1]+Z[1,k])))))+(1-rob.eps_2))
    
  #   # if j == 1
  #   #   obj_expr = tmp_prob_ind
  #   # else
  #   #   obj_expr += tmp_prob_ind
  #   # end

  #   @defNLExpr(pond_sum[j] == prod{pond_prod[l], l = 1:N})
  # end
  @defNLExpr(pond_ind[j=1:n_pt,k=1:N], (rob.eps_1+rob.eps_2-1)/((1+exp(-(unocc_cor[1,j]-(left_bd[1]+Z[1,k]))))*(1+exp(-(unocc_cor[2,j]-(bottom_bd[2]+Z[2,k]))))*(1+exp(-(unocc_cor[1,j]-(right_bd[1]+Z[1,k])))))+(1-rob.eps_2))

  # @setNLObjective(model, Min, sum{R[i]*(U[i,k] - URef[i,k])^2, k=1:N, i=1:nu} +
  #                 sum{Q[i]*(Z[i,k] - ZRef[i,k])^2, k=1:N, i=1:nz} +
  #                 sum{P[i]*(U[i,k+1] - U[i,k])^2, k=1:N-1, i=1:nu} +
  #                 sum{P[i]*(U[i,1] - u0[i])^2, i=1:nu})

  # @addNLConstraint(model, aux == obj_expr)
  # @setNLObjective(model, Min, sum{pond_sum[i],i=1:size(map.unocc_cor,2)})

  @setNLObjective(model, Min, sum{prob_map[unocc_cor[1,i],unocc_cor[2,i]]*prod{pond_ind[i,l],l=1:N},i=1:n_pt})

  ### Constraints
  for k=1:N
    # Dynamics
    if k==1
      @addNLConstraints(model, begin
                          Z[1,k] == z0[1] + dt*z0[4]*cos(z0[3])
                          Z[2,k] == z0[2] + dt*z0[4]*sin(z0[3])
                          Z[3,k] == z0[3] + dt*(U[1,k])
                          Z[4,k] == U[2,k]
                        end)
    else
      @addNLConstraints(model, begin
                          Z[1,k] == Z[1,k-1] + dt*Z[4,k-1]*cos(Z[3,k-1])
                          Z[2,k] == Z[2,k-1] + dt*Z[4,k-1]*sin(Z[3,k-1])
                          Z[3,k] == Z[3,k-1] + dt*(U[1,k])
                          Z[4,k] == U[2,k]
                        end)
    end
  end
  # note: distinguish addNLConstraints with addNLConstraint. former is used for group constraints while the latter is for single constraint.

  # for k=1:N
  #   # Dynamics
  #   if k==1
  #     @addNLConstraints(model, begin
  #                         Z[1,k] == z0[1] + dt*z0[4]*cos(z0[3]+U[1,k])
  #                         Z[2,k] == z0[2] + dt*z0[4]*sin(z0[3]+U[1,k])
  #                         Z[3,k] == z0[3] + dt*z0[4]/l_r*sin(U[1,k])
  #                         Z[4,k] == z0[4] + dt*U[2,k]
  #                       end)
  #   else
  #     @addNLConstraints(model, begin
  #                         Z[1,k] == Z[1,k-1] + dt*Z[4,k-1]*cos(Z[3,k-1]+U[1,k])
  #                         Z[2,k] == Z[2,k-1] + dt*Z[4,k-1]*sin(Z[3,k-1]+U[1,k])
  #                         Z[3,k] == Z[3,k-1] + dt*Z[4,k-1]/l_r*sin(U[1,k])
  #                         Z[4,k] == Z[4,k-1] + dt*U[2,k]
  #                       end)
  #   end
  # end

  # Input rate
  # @addNLConstraint(model, inputRateConsU0[i=1:nu],
  #                  U[i,1] - u0[i] <= dtMPC*duUB[i])
  # @addNLConstraint(model, inputRateConsL0[i=1:nu],
  #                  U[i,1] - u0[i] >= dtMPC*duLB[i])
  # @addNLConstraint(model, inputRateConsU[i=1:nu,k=1:N-1],
  #                U[i,k+1] - U[i,k] <= dt*duUB[i])
  # @addNLConstraint(model, inputRateConsL[i=1:nu,k=1:N-1],
  #                U[i,k+1] - U[i,k] >= dt*duLB[i])

  # Input bounds
  # @addNLConstraint(model, inputRateConsU[i = 1:nu, k=1:N],
  #                  U[i,k] <= duUB[i])
  # @addNLConstraint(model, inputRateConsL[i = 1:nu, k=1:N],
  #                  U[i,k] >= duLB[i])

  ### Solve dummy problem
  solve(model)

  # return MpcProblem(rob, uLB, uUB, duLB, duUB, nz, nu, tuning, model,
                        # Z, U, z0, u0, ZRef, URef)
  return MpcProblem(rob, map, uLB, uUB, nz, nu, tuning, model,
                        Z, U, z0, u0)
end

### Update and solve MPC optimization problem
# function updateSolveMpcProblem(mpc::MpcProblem, ZRef::Array{Float64,2},
#                                URef::Array{Float64,2}, z0::Array{Float64,1},
#                                u0::Array{Float64,1})
function updateSolveMpcProblem(mpc::MpcProblem, z0::Array{Float64,1},
                               u0::Array{Float64,1})
  # Parameters
  nz, nu, dt, N = mpc.nz, mpc.nu, mpc.tuning.dt, mpc.tuning.N

  # Warm start
  UPred = hcat(getValue(mpc.U[:,2:N]), getValue(mpc.U[:,N]))
  ZPred = simRobotModel(mpc.rob, UPred, dt)
  map(setValue, collect(mpc.Z[:,:]), collect(ZPred[:,:]))
  map(setValue, collect(mpc.U[:,:]), collect(UPred[:,:]))

  # Update model parameters
  mpc.z0[:] = z0
  mpc.u0[:] = u0
  mpc.ZRef[:,:] = ZRef
  mpc.URef[:,:] = URef

  # Solve updated problem
  tStart = time()
  solve(mpc.model)
  solveTime = time() - tStart

  return getValue(mpc.Z), getValue(mpc.U), solveTime
end

### Localize vehicle in current lane on map and compute reference speed
# function localizeVehicle(zCurr::Array{Float64,1},Map::TrackMap)
#   # Parameters
#   nodes, wayPointers, nLanes = Map.nodes, Map.wayPointers, Map.nWays
#   speed = Map.speed
#   nNodesPolyFront, nNodesPolyBack, nNodesThreshold = Map.nNodesPolyFront,
#     Map.nNodesPolyBack, Map.nNodesThreshold
#   xEgo, yEgo, psiEgo, vEgo = zCurr

#   # Polynomial fit for each lane
#   coeffsLanes = Inf*ones(4,nLanes)
#   idxClosestPointLanes = round(Int64,zeros(nLanes))
#   for n=1:nLanes
#     # nodes for each lane
#     nodesLane = nodes[:,wayPointers[1,n]:wayPointers[2,n]]
#     nNodesLane = size(nodesLane,2)

#     # distances from current position
#     dists = sqrt(sum((nodesLane-repmat([xEgo,yEgo],1,nNodesLane)).^2,1))

#     # distance and index of closest nodes
#     minDist, minIdx = findmin(dists)
#     idxClosestPointLanes[n] = wayPointers[1,n] + minIdx - 1
#     minIdx = min(minIdx, nNodesLane-1)
#     minIdx = max(minIdx, 2)

#     # direction of motion
#     dirEgo = [cos(psiEgo), sin(psiEgo)]
#     dirPos = nodesLane[:,minIdx+1] - nodesLane[:,minIdx]
#     direction = (dirEgo'*dirPos)[1] < 0 ? -1 : 1

#     # nodes for polynomial fit
#     idxStart = max(1, minIdx - (direction == 1 ? nNodesPolyBack : nNodesPolyFront))
#     idxEnd = min(nNodesLane, minIdx + (direction == 1 ? nNodesPolyFront : nNodesPolyBack))
#     nodesNear = nodesLane[:,idxStart:idxEnd]
#     nNodesNear = size(nodesNear,2)

#     if nNodesNear >= nNodesThreshold
#       # transform points to ego frame
#       R = [cos(psiEgo) -sin(psiEgo); sin(psiEgo) cos(psiEgo)]
#       nodesLocal = R'*(nodesNear - repmat([xEgo,yEgo],1,nNodesNear))

#       # compute least squares fit
#       px = nodesLocal[1,:]'
#       py = nodesLocal[2,:]'
#       H = [ones(nNodesNear) px px.^2 px.^3]
#       coeffsLanes[:,n] = -(H'*H)\H'*py
#     end
#   end

#   # closest lane index
#   minEy, laneIdx = findmin(abs(coeffsLanes[1,:]))
#   coeffsCurrLane = coeffsLanes[:,laneIdx]
#   mapSpeedRef = speed[idxClosestPointLanes[laneIdx]]

#   return coeffsCurrLane, laneIdx, mapSpeedRef
# end

# ### Generate state and input reference trajectory for MPC problem
# function generateReference(zCurr::Array{Float64,1},Map::TrackMap,
#                            tuning::Tuning,mpc::MpcProblem)
#   # parameters
#   dt, N = tuning.dt, tuning.N
#   xEgo, yEgo, psiEgo, vEgo = zCurr
#   nz, nu = mpc.nz, mpc.nu
#   xLocal = Map.xLocal

#   # localize vehicle in lane and get reference speed
#   coeffsCurrLane, laneIdx, vRef = localizeVehicle(zCurr,Map)

#   # longitudinal distance along lane centerline as a function of X
#   yLocal = polyval(Poly(-coeffsCurrLane),xLocal)
#   sLocal = [0.0; cumsum(sqrt(diff(xLocal).^2 + diff(yLocal).^2))]

#   # desired longitudinal positions
#   sDesired = [0.0; cumsum(dt*vRef*ones(N))]

#   # reference points on lane centerline
#   itpX = interpolate((sLocal,), xLocal, Gridded(Linear()))
#   itpY = interpolate((sLocal,), yLocal, Gridded(Linear()))
#   xLaneLocal, yLaneLocal = itpX[sDesired], itpY[sDesired]
#   posRefLocal = hcat(xLaneLocal,yLaneLocal)'

#   # reference points in global frame
#   R = [cos(psiEgo) -sin(psiEgo); sin(psiEgo) cos(psiEgo)]
#   posRefGlobal = repmat([xEgo,yEgo],1,N+1) + R*posRefLocal

#   # heading reference
#   psiRef = psiEgo*ones(N+1,1) + [0.0; atan2(diff(xLaneLocal),diff(yLaneLocal))]

#   # state reference
#   ZRef = [posRefGlobal; psiRef'; vRef*ones(1,N+1)]
#   URef = zeros(nu,N)

#   # outputs
#   return URef, ZRef
# end

end # module
