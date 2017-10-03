code files:
gameSim.m: the main file to run the program
Robot.m: robot class. the planning algorithm uses hand-coded SQP.
Robot2.m: robot class. same as Robot.m except the planning algorithm uses ipopt. not working actually.
Robot3.m: robot class. same as Robot2.m except that I use epigraph form for obj in ipopt. Now it could work in some cases, but maximum iteration exceeded also happens a lot. 
Robot4.m: robot class. same as Robot3.m except that I add variable gam_var. gamma actually causes much trouble resulting in iterates diverging and other issues. I'm trying to see if adding gam_var will help.
Robot5.m: robot class. same as Robot4.m except that I revert to original objective function without auxiliary variables. The P inverse in the coeefficient of Gaussian is assumed constant.
Robot6.m: robot class. same as Robot.m except that I now modify it to work with nonlinear measurement and target motion model. If this one works, then replace Robot.m.


labelResult.m : label the result with decision variable names and constraints. works with Robot.m.
labelResult2.m : same as labelResult.m except that this one contains the auxiliary variables for objective fucntion. works with Robot2.m - Robot4.m.
labelResult3.m : same as labelResult2.m except that soem constraints are changed. works with Robot5.m

folders:
matGeom and signed_distance are from Ashwin, using SQP for path planning.
sim_res stores the simulation results from the program. (not here yet)
figures stores the processed figures that may be used in the paper. (not here yet)
