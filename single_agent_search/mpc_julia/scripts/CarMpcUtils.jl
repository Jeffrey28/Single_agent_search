module CarMpcUtils

using MAT

export Vehicle, Tuning, Sigmoid, TrackMap
export kinematicBicycleDiscrete, doubleIntegratorDiscrete
export linearBicycleDiscrete, simVehicleModel, currentLaneNumber
export constantHeadingDiscrete, updateEgoState!
export rotateFrame, objectToVertices

####################
### Custom Types ###
####################

# Vehicle parameters
type Vehicle
  z::Array{Float64,1}   # state vector: [s, eY, ePsi, V]
  l_f::Float64
  l_r::Float64
  l_w::Float64   # half-width
  mass::Float64
  g::Float64
  Iz::Float64
  deltaMin::Float64
  deltaMax::Float64
  axMin::Float64
  axMax::Float64
  deltaRateLim::Float64
  axRateLim::Float64
  Cf::Float64
  Cr::Float64
  minSpeedTire::Float64
  ayMin::Float64
  ayMax::Float64
  model::Function   # state update model (discrete-time)

  function Vehicle(z::Array{Float64,1}; l_f::Float64 = 1.105,
                   l_r::Float64 = 1.738, l_w::Float64 = 0.8125,
                   mass::Float64 = 1830.0, g::Float64 = 9.81,
                   Iz::Float64 = 1544.0, deltaMin::Float64 = -10*pi/180,
                   deltaMax::Float64 = 10*pi/180, axMin::Float64 = -4.0,
                   axMax::Float64 = 1.0, deltaRateLim::Float64 = 10*pi/180,
                   axRateLim::Float64 = 2.0, Cf::Float64 = 4.0703e4,
                   Cr::Float64 = 6.4495e4, minSpeedTire::Float64 = 5.0,
                   ayMin::Float64 = -0.5, ayMax::Float64 = 0.5,
                   model::Function = x->x)
    new(z, l_f, l_r, l_w, mass, g, Iz, deltaMin, deltaMax, axMin, axMax,
        deltaRateLim, axRateLim, Cf, Cr, minSpeedTire, ayMin, ayMax, model)
  end
end

# MPC tuning parameters
type Tuning
  dt::Float64   # MPC model discretization time
  dtMPC::Float64   # MPC sampling time (computation)
  N::Int64   # Prediction horizon
  nSamples::Int64   # No. of samples for scenario MPC
  Q::Array{Float64,1}   # state
  R::Array{Float64,1}   # input
  P::Array{Float64,1}   # input rate
  vRef::Float64   # reference speed
  S::Array{Float64,1}   # slack
  dSafe::Float64   # longitudinal safety distance
  eYRef::Float64   # reference lateral position
  TTC::Float64   # time-to-collision
  eYSafe::Float64   # lateral safety distance

  function Tuning(; dt::Float64 = 0.2, dtMPC::Float64 = 0.2,
                  N::Int64 = 20, nSamples::Int64 = 20,
                  Q::Array{Float64,1} = [0.0],
                  R::Array{Float64,1} = [0.0],
                  P::Array{Float64,1} = [0.0],
                  vRef::Float64 = 0.0,
                  S::Array{Float64,1} = [0.0], dSafe::Float64 = 5.0,
                  eYRef::Float64 = 0.0, TTC::Float64 = 3.0,
                  eYSafe::Float64 = 0.5)
    new(dt, dtMPC, N, nSamples, Q, R, P, vRef, S, dSafe, eYRef, TTC, eYSafe)
  end
end

# Map
type TrackMap
  nodes::Array{Float64,2}
  dist::Array{Float64,1}
  speed::Array{Float64,1}
  wayPointers::Array{Int64,2}
  nWays::Int64
  GPSPos0::Array{Float64,1}
  laneWidth::Float64
  # tuning parameters for localization and polynomial fitting
  nNodesPolyFront::Int64
  nNodesPolyBack::Int64
  nNodesThreshold::Int64
  # local coordinatses for polynomial fitting
  xLocal::Array{Float64,1}

  function TrackMap(map_file::AbstractString; laneWidth::Float64=3.4, nNodesPolyFront::Int64=25,
               nNodesPolyBack::Int64=5, nNodesThreshold::Int64=5,
               xLocal::Array{Float64,1}=collect(0.0:0.1:100.0))
    file = matopen(map_file)
    inMap = read(file, "Map")
    nWays = round(Int64,inMap["N_Ways"])
    ways = inMap["Ways"]
    wayPointers = round(Int64,zeros(2,nWays))
    for i = 1:nWays
      wayPointers[:,i] = [ways[i][1,1], ways[i][1,end]]
    end
    new(inMap["Nodes"],vec(inMap["Dist"]),vec(inMap["Speed"]),wayPointers,
        nWays,vec(inMap["GPS_Pos0"]),laneWidth,
        nNodesPolyFront,nNodesPolyBack,nNodesThreshold,xLocal)
  end
end

# Sigmoid function
type Sigmoid
  a::Float64
  epsilon::Float64
  tol::Float64
  c::Float64

  function Sigmoid(; a::Float64 = 10.0, epsilon::Float64 = 1e-2,
                   tol::Float64 = 1e-2)
    #     c = 1/a*log((1-tol)/tol)   # under approximation
    c = 1/a*log(tol/(1-tol))   # over approximation
    new(a, epsilon, tol, c)
  end
end

######################
### Vehicle Models ###
######################

# Discrete-time kinematic bicycle model
function kinematicBicycleDiscrete(z::Array{Float64,1}, u::Array{Float64,1},
                                  dtPlant::Float64, veh::Vehicle)
  # z = [x, y, psi, v], u = [beta, ax]
  zNext = zeros(4)
  zNext[1] = z[1] + dtPlant*z[4]*cos(z[3]+u[1])
  zNext[2] = z[2] + dtPlant*z[4]*sin(z[3]+u[1])
  zNext[3] = z[3] + dtPlant*z[4]/veh.l_r*sin(u[1])
  zNext[4] = z[4] + dtPlant*u[2]

  return zNext
end

# Discrete-time double integrator model
function doubleIntegratorDiscrete(z::Array{Float64,1}, u::Array{Float64,1},
                                  dtPlant::Float64)
  # z = [x, y, xDot, yDot], u = [xDDot, yDDot]
  zNext = zeros(4)
  zNext[1] = z[1] + dtPlant*z[3]
  zNext[2] = z[2] + dtPlant*z[4]
  zNext[3] = z[3] + dtPlant*u[1]
  zNext[4] = z[4] + dtPlant*u[2]

  return zNext
end

function doubleIntegratorDiscrete(z::Array{Float64,1}, u::Array{Float64,1},
                                  dtPlant::Float64, veh::Vehicle)
  doubleIntegratorDiscrete(z, u, dtPlant)
end

# Discrete-time constant heading model
function constantHeadingDiscrete(z::Array{Float64,1}, u::Array{Float64,1},
                               dtPlant::Float64)
  # z = [x, y, psi, v], u = [ax]
  zNext = zeros(4)
  zNext[1] = z[1] + (dtPlant*z[4] + 0.5*dtPlant^2*u[1])*cos(z[3])
  zNext[2] = z[2] + (dtPlant*z[4] + 0.5*dtPlant^2*u[1])*sin(z[3])
  zNext[3] = z[3]
  zNext[4] = z[4] + dtPlant*u[1]

  return zNext
end

function constantHeadingDiscrete(z::Array{Float64,1}, u::Array{Float64,1},
                               dtPlant::Float64, veh::Vehicle)
  constantHeadingDiscrete(z, u, dtPlant)
end

# Discrete-time bicycle model with linear tire model
function linearBicycleDiscrete(z::Array{Float64,1}, u::Array{Float64,1},
                               dtPlant::Float64, veh::Vehicle)
  # z = [x, y, psi, xDot, yDot, psiDot], u = [delta, ax]
  psi = z[3]
  xDot = z[4]
  yDot = z[5]
  psiDot = z[6]
  delta = u[1]
  ax = u[2]

  # Slip angles
  alpha_f = -delta + atan((yDot + veh.l_f*psiDot)/max(xDot,veh.minSpeedTire))
  alpha_r = atan((yDot - veh.l_r*psiDot)/max(xDot,veh.minSpeedTire))

  # Tire forces
  Fyf = -veh.Cf*alpha_f
  Fyr = -veh.Cr*alpha_r

  # State update equations
  zNext = zeros(6)
  zNext[1] = z[1] + dtPlant*(xDot*cos(psi) - yDot*sin(psi))
  zNext[2] = z[2] + dtPlant*(xDot*sin(psi) + yDot*cos(psi))
  zNext[3] = z[3] + dtPlant*psiDot
  zNext[4] = z[4] + dtPlant*(yDot*psiDot + ax)
  zNext[5] = z[5] + dtPlant*(-xDot*psiDot + (2/veh.mass)*(Fyf + Fyr))
  zNext[6] = z[6] + dtPlant*((2/veh.Iz)*(a*Fyf - b*Fyr))

  return zNext
end

# Simulate vehicle model using current state and sequence of inputs
function simVehicleModel(f::Function, veh::Vehicle, U::Array{Float64},
                         dt::Float64)
  z0 = veh.z[:]
  nz = size(z0,1)
  N = size(U,2)
  Z = hcat(z0, zeros(nz,N))
  for k=1:N
    Z[:,k+1] = f(Z[:,k], U[:,k], dt, veh)
  end
  return Z[:,2:N+1]
end

function simVehicleModel(ego::Vehicle, U::Array{Float64}, dt::Float64)
  z0 = ego.z[:]
  nz = size(z0,1)
  N = size(U,2)
  Z = hcat(z0, zeros(nz,N))
  for k=1:N
    Z[:,k+1] = ego.model(Z[:,k], U[:,k], dt, ego)
  end
  return Z[:,2:N+1]
end

# Identify current lane
function currentLaneNumber(eY::Float64, laneWidth::Float64)
  floor((eY+0.5*laneWidth)/laneWidth) + 1
end

# Update ego vehicle state
function updateEgoState!(ego::Vehicle, u::Array{Float64,1}, dt::Float64)
  ego.z[:] = ego.model(ego.z, u, dt, ego)
end

###########################
### Object Manipulation ###
###########################

# Rotate points by given angle
function rotateFrame(points::Array{Float64}, angle::Float64)
  R = [cos(angle) -sin(angle); sin(angle) cos(angle)]
  return R*points
end

# Object vertices given center, orientation and dimensions
function objectToVertices(veh::Vehicle)
  l_f, l_r, l_w = veh.l_f, veh.l_r, veh.l_w
  cornersLocal = vcat([l_f, l_f, -l_r, -l_r]', l_w*[1.0, -1.0, -1.0, 1.0]')
  cornersGlobal = repmat(veh.z[1:2], 1, 4) +
    rotateFrame(cornersLocal, veh.z[3])
  return cornersGlobal
end

end # module
