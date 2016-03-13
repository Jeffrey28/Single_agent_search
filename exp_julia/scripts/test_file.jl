using JuMP
using Ipopt


# m = Model(solver=IpoptSolver(max_iter=100, max_cpu_time=1))
m = Model(solver=IpoptSolver(max_iter=100, max_cpu_time=1.0, print_timing_statistics="yes"))
n = 2
@defVar(m, x[1:n])
# @defNLExpr(m, y[i=1:3], sin(x[i]))
# @addNLConstraint(m, myconstr[i=1:3], y[i] <= 0.5)
@defNLExpr(myexpr[i=1:n], sin(x[i]))
# @defNLExpr(myexpr2[1:n])
# @defNLExpr(myexpr[1:n])
# @addNLConstraint(m, myconstr0[j=1:n], myexpr[j] == sin(x[j]))
@addNLConstraint(m, myconstr[j=1:n], myexpr[j] <= 0.5)
# @addNLConstraint(m, myconstr2[k=1:n], myexpr2[k] == myexpr[k]^2)
# @setNLObjective(m, Min, prod{y[1:3]})
# @setNLObjective(m, Min, prod{x[i],i=1:n})
@setNLObjective(m, Min, prod{myexpr[i],i=1:n})
# @setNLObjective(m, Min, prod{myexpr2[i],i=1:n})

solve(m)

# m = Model()
# @defVar(m, 0 <= x <= 1)
# # @defExpr(f1, -x)
# @defNLExpr(f2, -sin(x))

# # fs = [f1, f2]

# # @setObjective(m, Min, fs[1]) #Works
# # solve(m)

# @setNLObjective(m, Min, f2) #Works
# solve(m)

# # @setNLObjective(m, Min, fs[2]) #Works
# # solve(m)