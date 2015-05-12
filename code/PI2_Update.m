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
n_u  = ...
    
% The total number of parameters per input
n_theta = ...
    
% The length of the Upsilone basis function (time-varying basis funciton)
n_gaussian = ...
    
% Number of time steps
n_time = ...
    
% Number of rollouts per iteration
n_rollouts = ...



% The return of rollouts
R = ...


%% The exponentiated cost calculation is given 
% compute the exponentiated cost with the special trick to automatically
% adjust the lambda scaling parameter
minR = repmat( min(R,[],1), [n_rollouts 1]);
maxR = repmat( max(R,[],1), [n_rollouts 1]);
medR = repmat( median(R,1), [n_rollouts 1]);

% \exp(-1/lambda R)
expR = exp(-10*(R-minR)./(maxR-minR));


%% continue here with computing alpha
alpha = ...


%%
for i = 1:n_u % for each input i
   
    ...
    ...
    ...
  
    %% conversion back to vector-style
    delta_theta(:,i) = mat2vec(delta_theta_i);
    
end



%% additional functions vec2mat and mat2vec that are provided
function mat = vec2mat(vec)
    % Moving the time dimension (2nd dimension) to 3rd dimension
    % if the we have just one time instance, this will not change the 
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
