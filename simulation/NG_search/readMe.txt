10/31/17:

This folder is complicated. Robot.m - Robot7.m are all Robot class files. Robot.m is the basic one. Robot2-5.m use Ipopt as the solver, however, this does not work as Matlab keeps complaining model creation error or not solvable etc. From Robot2 to 5, I have tried different reformulation to make it easier for Yalmip, however, Robot5.m still cannot solve the problem. Robot6 and 7 are based on Robot.m and use my own hand-coded SQP. The problem can be solved now. When writing the Julia file to run Ipopt, I should better follow the formulation in Robot6.m first, and then in Robot7.m. 

code files:
gameSim.m: the main file to run the program.
Robot.m: robot class. the planning algorithm uses hand-coded SQP.

— Ipopt
Robot2.m: robot class. same as Robot.m except the planning algorithm uses ipopt. not working actually.
Robot3.m: robot class. same as Robot2.m except that I use epigraph form for obj in ipopt. Now it could work in some cases, but maximum iteration exceeded also happens a lot. 
Robot4.m: robot class. same as Robot3.m except that I add variable gam_var. gamma actually causes much trouble resulting in iterates diverging and other issues. I'm trying to see if adding gam_var will help.
Robot5.m: robot class. same as Robot4.m except that I revert to original objective function without auxiliary variables. The P inverse in the coefficient of Gaussian is assumed constant.

— hand-coded SQP
Robot6.m: robot class. same as Robot.m except that I now modify it to work with nonlinear measurement and target motion model. If this one works, then replace Robot.m.
Robot7.m : robot class. same as Robot6.m except that I now modify it to work with pedestrian motion model. Some states are not directly measurable.

labelResult.m : label the result with decision variable names and constraints. works with Robot.m.
labelResult2.m : same as labelResult.m except that this one contains the auxiliary variables for objective function. works with Robot2.m - Robot4.m.
labelResult3.m : same as labelResult2.m except that some constraints are changed. works with Robot5.m

simSetup.m: work with Robot.m - Robot6.m.
simSetup2.m (not modified yet): will be modified to be used with Robot7.m.

folders:
sim_res stores the simulation results from the program. (not here yet)
figures stores the processed figures that may be used in the paper. (not here yet)
