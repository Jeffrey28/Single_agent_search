code files:
gameSim.m: the main file to run the program
Robot.m: robot class. the planning algorithm uses hand-coded SQP.
Robot2.m: robot class. same as Robot.m except the planning algorithm uses ipopt. not working actually.
Robot3.m: robot class. same as Robot2.m except that I use epigraph form for obj in ipopt. Now it could work in some cases, but maximum iteration exceeded also happens a lot. 
Robot4.m: robot class. same as Robot3.m except that I add variable gam_var. gamma actually causes much trouble resulting in iterates diverging and other issues. I'm trying to see if adding gam_var will help.

folders:
matGeom and signed_distance are from Ashwin, using SQP for path planning.
sim_res stores the simulation results from the program. (not here yet)
figures stores the processed figures that may be used in the paper. (not here yet)
