function delta_theta = PI2_Update(Task,batch_sim_out,batch_cost)
%% PI2_Update
% delta_theta = PI2_Update(Task,batch_sim_out,batch_cost)
% 
%   INPUTS:
%       Task: struct with all relevant Task-specific parameters
%             .goal_time:  end time of simulation
%             .start_time: starting time of simulation
%             .dt:         time step for control inputs
%             .max_iteration_learning: maximum number of PI2 iterations
%             .num_rollouts: Number of rollouts for each batch (incl. repeated ones)
%             .num_reuse:    Number of solutions kept after each batch rollout
%             .n_basefct:    Number of base functions for smoothed-time controller
%             .std_noise:    standard deviation of exploration noise for PI2
%
%       batch_sim_out : a "num_rollouts x 1" struct array which contains
%       all rollout simulation data (of the perturbed system) in batch
%       format. It has the following fields
%              .t           vector of discrete simulation times for every
%              rollout
%              .x           state trajectory for every rollout
%              .u           input trajectory for every rollout
%              .Controller  Controller structure. Contains
%                           .BaseFnc(t, x), which
%                           returns a matrix of dimensions
%                           [num_states+1)*n_basefct, length(t)]
%                           indicating the basis function activation at
%                           every t in vector-fashion (one vector per
%                           time-step).
%                           !!! Output needs to be converted to matrix form using
%                           the function vec2mat().
%                           .theta : parameter matrix of dimensions
%                           [(n_gaussians*num_states+1), num_inputs )]
%              .eps         noise for parameter perturbation, of dimensions
%                           [(n_gaussians*num_states+1),num_inputs, num_timesteps]. 
%                           The noise gets reduced automatically as the 
%                           number of rollouts increases. Requires
%                           conversion using vec2mat()
%   
%        batch_cost :      matrix of dimensions (num_rollouts, num_timesteps)
%                           which holds the cost at every timestep for
%                           every rollout.
%               
%   OUTPUTS
%
%   delta_theta


%% Start your implementation here.
% Hint: start by extracting the following quantities

% Number of control inputs
n_u  = size(batch_sim_out(1).u,1);
    
% The total number of parameters per input
n_theta = size(batch_sim_out(1).Controller.theta,1); % (n+1)*numstates
    
% The length of the Upsilone ("Y") basis function (time-varying basis funciton)
n_gaussian = Task.n_gaussian;
    
% Number of time steps
n_time = (Task.goal_time-Task.start_time)/Task.dt + 1;
    
% Number of rollouts per iteration
n_rollouts = Task.num_rollouts;

n_states = size(batch_sim_out(1).x,1);

% For each control input i and time s:
% - Calculate the Return from starting time s for the kth rollout
% - Calculate alpha from starting time s for the kth rollout

% The return of rollouts
R = batch_cost; % size: n_rollouts x n_time

%% The exponentiated cost calculation is given 
% compute the exponentiated cost with the special trick to automatically
% adjust the lambda scaling parameter
minR = repmat( min(R,[],1), [n_rollouts 1]);
maxR = repmat( max(R,[],1), [n_rollouts 1]);
medR = repmat( median(R,1), [n_rollouts 1]);

% \exp(-1/lambda R)
expR = exp(-10*(R-minR)./(maxR-minR));
alpha = expR./repmat( sum(expR,1), [n_rollouts 1]);

% Used for intermediate calculations to finally come to delta_theta_i
delta_theta_i_s = zeros(n_gaussian, n_states+1, n_time);
% Same but not dependent on time.
delta_theta_i = zeros(n_gaussian, n_states+1);
% Also initialize complete delta_theta
delta_theta = zeros(n_gaussian*(n_states+1),n_u);
%%
for i = 1:n_u % for each input i
    for s = 1:n_time % for each time
        
        % Calculate the time varying parameter increment delta_theta_i_s
        sumDeltaTheta = zeros(n_gaussian,n_states+1);
        for k = 1:n_rollouts   
            x_k = batch_sim_out(k).x;
            t_k = batch_sim_out(k).t;
            
            epsilonVec = batch_sim_out(k).eps(:,i,s);
            epsilonMat = vec2mat(epsilonVec);
            
            upsilonBarVec = batch_sim_out(k).Controller.BaseFnc(t_k(s), x_k); % BaseFnc(t, x) % (n+1)*g x 1
            upsilonBarMat = vec2mat(upsilonBarVec); % g x (n+1)
            upsilonVec = upsilonBarMat(:,1); % g x 1
            
            
            upsEpsMatrix = (upsilonVec*upsilonVec')/(upsilonVec'*upsilonVec) * epsilonMat; % (g x g) * (g x (n + 1))
 
            sumDeltaTheta = sumDeltaTheta + alpha(k,s)*upsEpsMatrix; % scalar * matrix
        end
        
        delta_theta_i_s(:,:,s) = sumDeltaTheta;
    end % for..s
    
    % Time-averaging the parameter vector
    x_1 = batch_sim_out(1).x;
    t_1 = batch_sim_out(1).t;
    upsilonBarVec = batch_sim_out(1).Controller.BaseFnc(t_1, x_1); % BaseFnc(t, x) % g*(n+1) x t
    upsilonBarMat = vec2mat(upsilonBarVec); % g x (n+1) x t

    delta_theta_i = sum((upsilonBarMat .* delta_theta_i_s),3) ./ sum(upsilonBarMat,3);
    
%     for j = 1:size(delta_theta_i_s,2) % For the jth column of delta_theta_i_s matrix    
%         x_1 = batch_sim_out(1).x;
%         t_1 = batch_sim_out(1).t;
%         
%         % Time-averaging the parameter vector
%         
%         % Calculate integrand to assign delta_theta_i_s
%         integrandDeltaTheta = zeros(n_gaussian,1);
%         normalizationVector = zeros(n_gaussian,1);
%         
%         for s = 1:n_time
%             upsilonBarVec = batch_sim_out(1).Controller.BaseFnc(t_1(s), x_1); % BaseFnc(t, x)
%             upsilonBarMat = vec2mat(upsilonBarVec); % g x (n+1)
% 
%             upsilonVec_S = upsilonBarMat(:,1); % g x 1
%             
%             normalizationVector = normalizationVector + upsilonVec_S;
%         
%             integrandDeltaTheta = integrandDeltaTheta + delta_theta_i_s(:,j,s) .* upsilonVec_S;  
%         end
%         
%         delta_theta_i(:,j) = integrandDeltaTheta ./ normalizationVector;        
%     end % for..j

    %% conversion back to vector-style
    delta_theta(:,i) = mat2vec(delta_theta_i);
end


%% additional functions vec2mat and mat2vec that are provided
function mat = vec2mat(vec)
    % Moving the time dimension (2nd dimension) to 3rd dimension
    % if we have just one time instance, this will not change the 
    % dimension of vec
    vec = permute(vec, [1 3 2]);
    
    mat = reshape(vec,[],Task.n_gaussian,size(vec,3));
    mat = permute(mat, [2 1 3]);
end


function vec = mat2vec(mat)
    
    mat = permute(mat, [2 1 3]);
    vec = reshape(mat,[],1,size(mat,3));
end


end
