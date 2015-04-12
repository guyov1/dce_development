function [ ] = Run_DCE_Perfusion( Subject_name, Subject_Path, Sim_Struct, PefusionOutput, Verbosity )
%Run_DCE_Perfusion Run Full DCE Perfusion Analysis

display('---------------------------------------------------');
display('-I- Started Run_DCE_Perfusion execution at:');
c = clock;
disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
display('---------------------------------------------------');

DCECoregP             = Subject_Path;
WM_mask_absolute_path = [Subject_Path  '\RefT1_WM_830.nii'];
Art_Mask              = [Subject_Path  '\InspectedRepVox.nii'];
Vein_Mask             = [Subject_Path  '\Veins_Mask.nii'];
After_CTC_mat         = [Subject_Path  '\AfterCTC.mat'];
Brain_Extract_path    = [Subject_Path  '\Manual_BrainMask.nii'];

% Brain Mask
Brain_Mask_3D   = loadniidata(Brain_Extract_path);

% Override real data flag and parameters range
Sim_Struct.RealData_Flag              = true;
Sim_Struct.Vb_low                     = 0.1; % [mL/100g]     , they used 3,6,12,18
Sim_Struct.Vb_max                     = 100;
Sim_Struct.Ve_low                     = 0.1; % Must be smaller than Vtis
Sim_Struct.Ve_max                     = 100;
Sim_Struct.LowerBound_Larsson         = [Sim_Struct.Vb_low Sim_Struct.E_low  Sim_Struct.Ve_low];
Sim_Struct.UpperBound_Larsson         = [Sim_Struct.Vb_max Sim_Struct.E_max  Sim_Struct.Ve_max];
Sim_Struct.init_Ve_guess              = 0.1;

% Set parallel processing if needed
Set_Parallel_Processing(Sim_Struct, Verbosity);

% Set output directory for figures/report
Output_directory      =  [PefusionOutput 'Run_Output/'];
% Create directory if does not exist
if ~exist(Output_directory,'dir')
    mkdir(Output_directory);
end

% Add the needed data from the DCE run
display('-I- Loading .mat files from DCE run...');

if(exist(After_CTC_mat,'file'))
    load(After_CTC_mat);
else
    error('-E- AfterCTC.mat does not exist...');
end

% Define needed parameters from DCE data
Sim_Struct.num_time_stamps  = size(CTC2D,2);
Sim_Struct.num_voxels       = size(CTC2D,1);
Sim_Struct.sec_interval     = TimeBetweenDCEVolsFinal;
Sim_Struct.min_interval     = Sim_Struct.sec_interval / 60;
time_vec_minutes            = Sim_Struct.min_interval * ( 0 : Sim_Struct.num_time_stamps - 1 );
Sim_Struct.time_vec_minutes = time_vec_minutes;

%% Create AIF
if exist('AIFFindData_mat','var')
    [AIF_Struct] = chooseAifForRealData(Sim_Struct, CTC2D, Art_Mask, Vein_Mask, Msk2, Output_directory, AIFFindData_mat);
else
    [AIF_Struct] = chooseAifForRealData(Sim_Struct, CTC2D, Art_Mask, Vein_Mask, Msk2, Output_directory);
end

%% Take Ct and AIF and calculate Ht
Ct               = CTC2D(:,:);
num_total_voxels = size(Ct,1);

% Choose the AIF (either parametric or from ICA average)
Chosen_AIF = AIF_Struct.AIF_estimated_ICA; % AIF_paramtertic, transpose(smooth(AIF_estimated_ICA))

% Scale AIF as necessary
Chosen_AIF = double( Sim_Struct.AIF_Scaling_Factor * Chosen_AIF );

[ resultStruct ] = Get_Ht_Deconvolving(Sim_Struct, Chosen_AIF, Ct , Output_directory, Subject_name, Sim_Struct.Force_RealData_Calc, Verbosity);

%% Save and Write results

display('-I- Saving parameters result of Get_Ht_Deconvolving...');

Mat_File_To_Save = [Output_directory 'All_Parameters_Result.mat'];

save(Mat_File_To_Save,'resultStruct','Msk2','WorkingP','PefusionOutput','num_total_voxels','time_vec_minutes','TimeBetweenDCEVolsFinal','time_vec_minutes'...
    ,'Chosen_AIF', 'DCECoregP','Sim_Struct','WM_mask_absolute_path','Subject_Path');

% Write only in case we used the entire brain (say, above 1000 voxels)
if (num_total_voxels > 1000)
    resultToNiiAndReport(resultStruct, time_vec_minutes, CTC2D, Chosen_AIF, Msk2, Brain_Mask_3D, Output_directory, DCECoregP, Sim_Struct, WM_mask_absolute_path);
end

%% Create PDF Report
LogFN = [WorkingP 'Log.mat'];
MakeReport_func(Output_directory, LogFN);
close all;

% Display Ct(t) and fit of a single voxel
% plotSingleCTC(187, 149, 2, CTC_4D, conv_result_IRF_4D)

display('---------------------------------------------------');
display('-I- Finished Run_DCE_Perfusion execution at:');
c = clock;
disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
display('---------------------------------------------------');

end

