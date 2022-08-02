# CVPM-MPC 

Welcome to the repository of my Bachelor's thesis at the Department of Electrical and Computer Engineering at TUM. My thesis was about "Constraint Violation Probability Minimization for a Robot Manipulator". To read a more detailed discussion see my thesis and final presentation.

This directory contains the Matlab code, as well as the final [report](https://github.com/sirine90/CVPM-MPC/blob/main/report/Bachelor-Thesis.pdf) and [presentation](https://github.com/sirine90/CVPM-MPC/blob/main/presentation/StudentSlides_Template.pdf).

This is a brief description of some of the functions under the folder "functions":


-"build" is used for building the configuration space and defining the collision-free probabilistic set.

-"rect" defines a rectangular obstacle in the workspace that needs to be transformed to a C-obstacle in c-space.

-"createRobot" defines the robot parameters using the Robotics Toolbox.

-"sys_trajectory" returns the reference trajectory in the configuartion space.

-"mySystem" Class defines the robot's forward dynamics and inverse dynamics computations.

-"myCVPM" Class defines the CVPM class for the robot manipulator. 





