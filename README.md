# Data-Driven-System-Identifiacation-for-LTI-systems
In this repository, algorithms for system identification are developed and implemented for LTI systems.

In SysID.m, a full algorithm for LTI system identification can be found in the first section of the code. It supports system with any order and input & output dimensions, unknown initial condition and random inputs (to make sure it's persistency excitation). It also works for corrupted measurements or even nonlinear systems by capturing its linear dynamics to approximate. 

The second and the third sectoins are for validation, you may examine how accurate the system's dynamics can be reconstructed by the algorithm by comparing to the original input&output data.

You may use LTI_Data_generator.m to randomly generate a discrete LTI system with any order and input&output dimensions to produce its input&output data or even corrupted ones to see the robustness of the algorithm.
Hankel_Matrix.m is a must function to build Hankel Matrices when solving system identification problems.

The algorithm for the whole sytem ID could be broken down into a few steps and there are a few algorithms for solving less general system ID problems with different response scenarios.

FP_cal_N4SID.m : it finds A,C matrices of a system from a response with random inputs and unknown initial conditions by calculating the free response of the system using N4ISD algorithm. No constraints on the LTI system.

HankelfromImpulse.m : It reconstruted A,B,C,D matrices for an impulse response, assuming the system has zero D matrix. 

