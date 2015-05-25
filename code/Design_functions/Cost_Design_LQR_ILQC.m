function Cost = Cost_Design( m_quad, Task)
%COST_DESIGN Creates a cost function for LQR and for ILQC
%   .Q_lqr 
%   .R_lqr
%
%   A ILQC cost J = h(x) + sum(l(x,u) is defined by:
%   .h - continuous time terminal cost
%   .l - continous time intermediate cost

global ILQC_TYPE;
global ONLY_LQR;
global ILQC_VIAPOINT_2;

% Quadcopter system state x
syms qxQ qyQ qzQ qph qth qps dqxQ dqyQ dqzQ dqph dqth dqps real
x_sym = [ qxQ qyQ qzQ ...    % position x,y,z
          qph qth qps ...    % roll, pitch, yaw
          dqxQ dqyQ dqzQ ... % velocity x,y,z
          dqph dqth dqps ]'; % angular velocity roll, pitch, yaw     
% Quadcopter input u (Forces / Torques)
syms Fz Mx My Mz real; u_sym = [ Fz Mx My Mz ]';
% Time variable for time-varying cost
syms t_sym real;
Cost = struct;


%% LQR cost function
% LQR controller
% Q needs to be symmetric and positive semi definite
Cost.Q_lqr = 5*diag([  1   1    1   ...   % penalize positions
                       3   3    3   ...   % penalize orientations
                       0.1  0.1    2   ...   % penalize linear velocities
                       1   1    1 ]);     % penalize angular velocities         
% R needs to be positive definite
Cost.R_lqr = 10*diag([1 1 1 1]);        % penalize control inputs

%% ILQC cost function
Cost.Qm  = Cost.Q_lqr;                  % it makes sense to redefine these
Cost.Rm  = Cost.R_lqr;
Cost.Qmf = Cost.Q_lqr;

% Reference states and input the controller is supposed to track
gravity = 9.81;
f_hover = m_quad*gravity; % keep input close to the one necessary for hovering
Cost.u_eq = [ f_hover ; zeros(3,1) ];
Cost.x_eq = Task.goal_x;

% alias for easier referencing
x = x_sym; 
u = u_sym;
x_goal = Task.goal_x;

%ilqc_type = 'goal_state'; % Choose 'goal_state or 'via_point'
ilqc_type = 'via_point';

if ~isempty(ILQC_TYPE)
    ilqc_type = ILQC_TYPE;
end
if (isempty(ONLY_LQR) || ...
   (ONLY_LQR == false) )
    fprintf('ILQC cost function type: %s \n', ilqc_type);
end

switch ilqc_type
    case 'goal_state'     
        Cost.h = simplify((x-x_goal)'*Cost.Qmf*(x-x_goal));
        Cost.l = simplify( ...
                           (x-Cost.x_eq)'*Cost.Qm*(x-Cost.x_eq) + ...
                           (u-Cost.u_eq)'*Cost.Rm*(u-Cost.u_eq) );
    case 'via_point'
        %% Problem 2.2: Include a waypoint p1 in the ILQC cost function formulation
        if ~isempty(ILQC_VIAPOINT_2) && ILQC_VIAPOINT_2 == true        
        p1 = Task.vp2;      % p1 = Task.vp2 also try this one
        else
        p1 = Task.vp1;      % p1 = Task.vp2 also try this one
        end
        t1 = Task.vp_time;
        
        % Define an approriate weighting for way points (see script or Eq.(5))
        % Hint: Which weightings must be zero for the algorithm to
        % determine optimal values?
        
        % Q_vp is defined as a symmetric, diagonal matrix
        Q_vp = diag([ 3   3    3   ...   % penalize positions
                      1   1    1   ...   % penalize orientations
                      0   0    0   ...   % DO NOT penalize linear velocities (such that iLQC can determine optimal velocities)
                      0   0    0 ]);     % DO NOT penalize angular velocities  

        % don't penalize position deviations, drive system with final cost
        Cost.Qm(1:3,1:3) = zeros(3); 
              
        % Define symbolic cost function. Use fnct "viapoint(.)" below
        Cost.h = simplify((x-x_goal)'*Cost.Qmf*(x-x_goal));
        Cost.l = simplify( ...
                           (x-Cost.x_eq)'*Cost.Qm*(x-Cost.x_eq) + ...
                           (u-Cost.u_eq)'*Cost.Rm*(u-Cost.u_eq) + ...
                            viapoint(t1, p1, x, t_sym*Task.dt, Q_vp) ); % function call: viapoint(vp_t,vp,x,t,Qm_vp)); Convert time t from timesteps to seconds
    
      otherwise
        error('Unknown ilqc cost function mode');
        
end 

Cost.x  = x;
Cost.u  = u;
Cost.t  = t_sym; 
end


function viapoint_cost = viapoint(vp_t,vp,x,t,Qm_vp)
% WAYPOINT Creates cost depending on deviation of the system state from a
% desired state vp at time t. Doesn't need to be modified. 
% For more details see:
% http://www.adrl.ethz.ch/archive/p_14_mlpc_quadrotor_cameraready.pdf
%
%    vp_t  :   time at which quadrotor should be at viapoint wp
%    vp    :   position where quad should be at given time instance
%    x,t   :   symbolic references to state and time
%    Qm_vp :   The weighting matrix for deviation from the via-point

prec = 3; % how 'punctual' does the quad have to be? The bigger the
          % number, the harder the time constraint
          
viapoint_cost = (x-vp)'*Qm_vp*(x-vp) ...
                *exp(-0.5*prec*(t-vp_t)^2) /sqrt(2*pi/prec);

end

