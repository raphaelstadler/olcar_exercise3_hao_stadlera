function [LQR_Controller, Cost] = LQR_Design(Model,Task)
% LQR_DESIGN Creates an LQR controller to execute the Task
% Controller: Controller structure
%             u(x) = theta([t]) ' * BaseFnc(x)
%             .theta:   Parameter matrix (each column is associated to one
%                       control input)
%             .BaseFnc
%             .time
global LQR_TYPE;

%% Problem 1.1: Find optimal feedback gains according to an LQR controller
% Make use of the elements in Task.cost for linearization points and cost
% functions.

% Define the state input around which the system dynamics should be linearized    
x_lin = Task.cost.x_eq;    % 12x1: equilibrium position x,y,z; roll, pitch, yaw; velocity x,y,z; angular rates roll, pitch, yaw
u_lin = Task.cost.u_eq;    % 4x1: equilibrium forces/torques Fz, Mx, My, Mz % assume that most of the time the quadrotor flies similar to hovering (Fz != 0)

% The linearized system dynamics matrices A_lin B_lin describe the dynamics
% of the system around the equilibrium state. See exercise sheet Eq.(2) for 
% usage example
A_lin = Model.Alin{1}(x_lin, u_lin, Model.param.syspar_vec);
B_lin = Model.Blin{1}(x_lin, u_lin, Model.param.syspar_vec);


% Compute optimal LQR gain (see command 'lqr')
%
% J = Integral {x'Qx + u'Ru + 2*x'Nu} dt
% subject to the system dynamics  dx/dt = Ax + Bu
%
% Solution fullfills: u* = -K*x
K = lqr(A_lin,...           % A
        B_lin,...           % B
        Task.cost.Q_lqr,... % Q
        Task.cost.R_lqr);   % R; Note: N is supposed to be zero if not mentioned implicitely

%% Design the actual controller using the optimal feedback gains K
% The goal of this task is to correctly fill the matrix theta.
LQR_Controller.BaseFnc = @Augment_Base;
LQR_Controller.time    = Task.start_time:Task.dt:(Task.goal_time-Task.dt);
% The quadrotor controller produces inputs u from a combination of 
% feedforward uff and feedback elements as follows:
%
%   u = [Fz, Mx, My, Mz]' = uff + K'(x_ref - x)
%                         = uff + K'x_ref - K'x
%                         = [uff + K'x_ref, -K' ]' * [1,x']'
%                         =        theta'          * BaseFnc

%lqr_type = 'goal_state';  % Choose 'goal_state or 'via_point'
lqr_type = 'via_point';  % Choose 'goal_state or 'via_point'

if ~isempty(LQR_TYPE)
    lqr_type = LQR_TYPE;
end

fprintf('LQR controller design type: %s \n', lqr_type);
Nt = ceil((Task.goal_time - Task.start_time)/Task.dt+1);   
switch lqr_type
    case 'goal_state'
        %% Problem 1.2: Drive system to Task.goal_x with LQR Gain + feedforward
        % Define equilibrium point x_eq correctly and use it to generate
        % feedforward torques (see exercise sheet Eq.(3)).
        
        x_ref = Task.goal_x;   % reference point the controller tries to reach
        
        uff = Task.cost.u_eq;     % duff = uff - u_lin; Since duff = 0, uff = u_lin; same dimensions as u
        % Feedback control integrated in theta (see below): 
        % u_fb = - K * (x - x_ref) = K * (x_ref - x)
        % Calculate theta
        theta = [uff + K*x_ref, -K]'; % 13 x 4
        
        % stack constant theta matrices for every time step Nt
        LQR_Controller.theta = repmat(theta,[1,1,Nt]);
        
    case 'via_point'  
        %% Problem 1.3: Drive system to Task.goal_x through via point p1.
        % (see exercise sheet Eq.(3))
        
        t1 = Task.vp_time;
        p1 = Task.vp1; % steer towards this input until t1
        uff = Task.cost.u_eq;     % duff = uff - u_lin; Since duff = 0, uff = u_lin; same dimensions as u
        
        for t=1:Nt-1
           if t <= t1/Task.dt   % steer towards via point
               temp_goal = p1;
           else                 % steer towareds final goal point
               temp_goal = Task.goal_x;   % reference point the controller tries to reach
           end
           
           % Calculate theta similar as before, by using temp_goal instead of x_ref:
           theta = [uff + K*temp_goal, -K]'; % 13 x 4

           % varying theta matrices for every time step Nt 
           LQR_Controller.theta(:,:,t) = theta;
        end
        
    otherwise
        error('Unknown lqr_type');
        
end

% Calculate cost of rollout
sim_out = Quad_Simulator(Model, Task, LQR_Controller);
Cost = Calculate_Cost(sim_out, Task);

end


function x_aug = Augment_Base(t,x)
% AUGMENT_BASE(t,x) Allows to incorporate feedforward and feedback
%
%   Reformulating the affine control law allows incorporating
%   feedforward term in a feedback-like fashion:
%   u = [Fz, Mx, My, Mz]' = uff + K(x - x_ref)
%                         = uff - K*x_ref + K*x
%                         = [uff - K*x_ref, K ] * [1,x']'
%                         =        theta'       * BaseFnc 

number_of_states = size(x,2);      % multiple states x can be augmented at once
x_aug = [ones(1,number_of_states); % augment row of ones
                    x           ];
end


function Cost = Calculate_Cost(sim_out, Task)
% Calculate_Cost(.) Asses the cost of a rollout sim_out 
%
%   Be sure to correctly define the equilibrium state (X_eq, U_eq) the 
%   controller is trying to reach.
X = sim_out.x;
U = sim_out.u;
Q = Task.cost.Q_lqr;
R = Task.cost.R_lqr;
X_eq = repmat(Task.cost.x_eq,1,size(X,2)-1); % equilibrium state LQR controller tries to reach
U_eq = repmat(Task.cost.u_eq,1,size(U,2));   % equilibrium input LQR controller tries to reach
Ex = X(:,1:end-1) - X_eq;                    % error in state
Eu = U - U_eq;                               % error in input

Cost = Task.dt * sum(sum(Ex.*(Q*Ex),1) + sum(Eu.*(R*Eu),1));
end





