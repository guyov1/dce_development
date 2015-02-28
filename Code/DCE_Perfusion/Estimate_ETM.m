function [ Sim_Struct ] = Estimate_ETM( Sim_Struct, Verbosity )

if ~strcmp(Verbosity,'None')
    display('-I- Starting Murase/Tofts Comparison...');
end

% Take from struct variables used in local function
min_interval                  = Sim_Struct.min_interval;
time_vec_minutes              = Sim_Struct.time_vec_minutes;
Sim_AIF_delayed_no_noise      = Sim_Struct.Sim_AIF_delayed_no_noise;
Sim_AIF_with_noise            = Sim_Struct.Sim_AIF_with_noise;
SNR_ratio                     = Sim_Struct.SNR_ratio;
One_Iteration_Murase_Tofts    = Sim_Struct.One_Iteration_Murase_Tofts;
Iterate_Murase_Tofts_num_iter = Sim_Struct.Iterate_Murase_Tofts_num_iter;
algorithm_options             = Sim_Struct.algorithm_options;
adjusted_larsson              = Sim_Struct.Adjusted_Larsson_Model;
num_iterations                = Sim_Struct.num_iterations;
Sim_Ct_larss_Murase_noise     = Sim_Struct.Sim_Ct_larss_kernel_noise;

% Check if only one iteration is wanted
if (One_Iteration_Murase_Tofts)
    Iterate_Murase_Tofts_num_iter = 1;
end

% AIF index to choose for iteration
AIF_idx = 1;


% Initiate parameter matrices
Murase_params    = zeros(3,num_iterations);

% Initiate time calculation for first index
tic;
display(sprintf('-I- Starting simulation for %d voxels...',Iterate_Murase_Tofts_num_iter));

%parfor idx = 1 : Iterate_Murase_Tofts_num_iter
for idx = 1 : num_iterations
    
    % Create vectors of A matrix
    A_1 =  cumtrapz(time_vec_minutes,Sim_AIF_with_noise(:,AIF_idx)');
    A_2 = -cumtrapz(time_vec_minutes,Sim_Ct_larss_Murase_noise(:,idx)');
    A_3 =  Sim_AIF_with_noise(:,AIF_idx)';
    A   =  [A_1' A_2' A_3'];
   
    % Create c vector
    C_vec = Sim_Ct_larss_Murase_noise(:,idx);
    
    B     = A \ C_vec;
    %B = pinv(A) * C_vec;
    
    K_trans_Tofts_Murase = B(1) - B(2)*B(3);
    K_ep_Tofts_Murase    = B(2);
    Vp_Tofts_Murase      = B(3);
    
    
    Murase_params(:,idx)    = [K_trans_Tofts_Murase K_ep_Tofts_Murase Vp_Tofts_Murase];
    
end


display(sprintf('Finished simulation for %d voxels...',Iterate_Murase_Tofts_num_iter));
time_finish = toc;
display(sprintf('Took %.2f seconds to finish...',time_finish));
tic;

Sim_Struct.Est_Ktrans_vec = Murase_params(1,:);
Sim_Struct.Est_Kep_vec    = Murase_params(2,:);
Sim_Struct.Est_Vp_vec     = Murase_params(3,:);
Sim_Struct.Est_Ve_vec     = Sim_Struct.Est_Ktrans_vec ./ Sim_Struct.Est_Kep_vec;

if strcmp(Verbosity,'Full')
    display('-I- Finished Murase/Tofts Comparison...');
end

end