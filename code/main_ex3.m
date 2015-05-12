%% OLCAR - Exercise 3: PI2 Learning

close all; clc; 

addpath(genpath(pwd)); % adds folders and subfolders to path

% To generate plots of LQR/ILQG rollout
plotvec = {'quad_pos_noLoad','quad_angles_noLoad','control_input_noLoad','rotor_thrust_noLoad'};
create_pdf = [0 0 0 0 ];    % for which plots should a pdf be created
plot_ind   = [1 2 3 4 ];    % which data to plot on screen


%% Task definition
Task = Task_Design();

%Load the nominal model of the quadcopter
load('Quadrotor_Model.mat','Model');

% Define cost function
Task.cost = Cost_Design( Model.param.mQ, Task );


%% Initial controller design - fill in your ILQC_Design here
[Initial_Controller, Cost_LQR] = LQR_Design(Model, Task);
[ILQC_Controller, Cost] = ...
clear Model;


%% Transform controller from ILQC representation to basis function representation
ReducedController = BaseFcnTrans(ILQC_Controller,Task.n_gaussian);


%% Load the perturbed "real" quadrotor model
load('Quadrotor_Model_perturbed.mat','Model_perturbed'); 


%% Visualization of initial guess with policy in base function representation on perturbed system
Task.noise_insertion_method = '';
Test_sim_out = Sample_Rollout(Model_perturbed,Task,ReducedController);
fprintf('Final Quadcopter Position: xQ = %.3f, yQ = %.3f, zQ = %.3f \n',Test_sim_out.x(1:3,end));
fprintf('Final Quadcopter Velocity: xQ = % .3f, yQ = %.3f, zQ = %.3f \n',Test_sim_out.x(7:9,end));
Visualize(Test_sim_out,Model_perturbed.param,'plot_mode',1);
%Plot_Result(Test_sim_out,Model_perturbed.param,'plots',plotvec(plot_ind),'file',create_pdf(plot_ind),'path',pwd)

 
%% Start PI Learning
t_cpu = cputime; 
[LearnedController,AllCost,AllController] = PIs_Learning(Model_perturbed,Task,ReducedController);
t_cpu = cputime - t_cpu;  
fprintf('CPU time: %f \n',t_cpu);
fprintf('The PI² algorithm took %fs to converge \n\n',t_cpu);


%% Visualization of PI² final trajectory
Task.random=[];
Test_sim_out = Sample_Rollout(Model_perturbed,Task,LearnedController);
fprintf('Final Quadcopter Position: xQ = %.3f, yQ = %.3f, zQ = %.3f \n',Test_sim_out.x(1:3,end));
fprintf('Final Quadcopter Velocity: xQ = %.3f, yQ = %.3f, zQ = %.3f \n',Test_sim_out.x(7:9,end));
Visualize(Test_sim_out,Model_perturbed.param,'plot_mode',1);
%Plot_Result(Test_sim_out,Model_perturbed.param,'plots',plotvec(plot_ind),'file',create_pdf(plot_ind),'path',pwd)


%% Cost plot
figure()
plot(AllCost)
xlabel('PI2 iteration number');
ylabel('Cost')
title('Quadrotor PI Learning Curve')