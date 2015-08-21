

if Ignore_Delay_Model_Selection
    %%
    NParams          = [0 1 2 3 4];
    num_maps         = length(NParams);
    size_3D          = size(RMS_Larss_no_Delay_3D);
    RMSmaps          = zeros([size_3D num_maps]);
    
    RMSmaps(:,:,:,5) = RMS_Larss_with_Delay_3D;             % 5 Params
    %RMSmaps(:,:,:,8) = RMS_Larss_no_Delay_3D;               % 4 Params
    if useUptakeNoETM
        RMSmaps(:,:,:,4) = RMS_Larss_with_Delay_no_Ve_3D;      % 4 Params -----
    else
        RMSmaps(:,:,:,4) = RMS_Larss_with_Delay_High_F_3D;      % 4 Params -----
    end

    %RMSmaps(:,:,:,6) = RMS_Larss_no_Delay_High_F_3D;        % 3 Params
    RMSmaps(:,:,:,3) = RMS_Larss_with_Delay_no_E_3D;        % 3 Params -----
    %RMSmaps(:,:,:,4) = RMS_Larss_no_Delay_no_E_3D;          % 2 Params
    RMSmaps(:,:,:,2) = RMS_Larss_with_Delay_no_E_High_F_3D; % 2 Params -----
    %RMSmaps(:,:,:,2) = RMS_Larss_no_Delay_no_E_High_F_3D;   % 1 Params
    RMSmaps(:,:,:,1) = RMS_Larss_no_Delay_zero_params_3D;   % 0 Params
    
    [ ChosenByAIC_3D ]                      = Model_Selection( nDataPoints, RMSmaps, NParams, cur_Data_Weight, AIC_Correction);
    
    % F
    F_Model_Selected_3D                     = zeros(size(Flow_with_Delay_3D));
    F_Model_Selected_3D (ChosenByAIC_3D==5) = Flow_with_Delay_3D     (ChosenByAIC_3D==5);
    if useUptakeNoETM
        F_Model_Selected_3D (ChosenByAIC_3D==4) = Flow_with_Delay_3D (ChosenByAIC_3D==4); % Put zero although its Inf
    else
        F_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map          (ChosenByAIC_3D==4); % Put zero although its Inf
    end
    F_Model_Selected_3D (ChosenByAIC_3D==3) = Flow_with_Delay_3D   (ChosenByAIC_3D==3);
    F_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2); % Put zero although its Inf
    F_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % Vb
    Vb_Model_Selected_3D                    = zeros(size(Vb_no_Delay_3D));
    Vb_Model_Selected_3D(ChosenByAIC_3D==5) = Vb_with_Delay_3D             (ChosenByAIC_3D==5);
    if useUptakeNoETM
        Vb_Model_Selected_3D(ChosenByAIC_3D==4) = Vb_with_Delay_no_Ve_3D      (ChosenByAIC_3D==4);
    else
        Vb_Model_Selected_3D(ChosenByAIC_3D==4) = Vb_with_Delay_High_F_3D      (ChosenByAIC_3D==4);
    end
    Vb_Model_Selected_3D(ChosenByAIC_3D==3) = Vb_with_Delay_no_E_3D        (ChosenByAIC_3D==3);
    Vb_Model_Selected_3D(ChosenByAIC_3D==2) = Vb_with_Delay_no_E_High_F_3D (ChosenByAIC_3D==2);
    Vb_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % Ve
    Ve_Model_Selected_3D                    = zeros(size(Ve_no_Delay_3D));
    Ve_Model_Selected_3D(ChosenByAIC_3D==5) = Ve_with_Delay_3D             (ChosenByAIC_3D==5);
    if useUptakeNoETM
        Ve_Model_Selected_3D(ChosenByAIC_3D==4) = zeros_map                (ChosenByAIC_3D==4);
    else
        Ve_Model_Selected_3D(ChosenByAIC_3D==4) = Ve_with_Delay_High_F_3D  (ChosenByAIC_3D==4);
    end
    
    Ve_Model_Selected_3D(ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    Ve_Model_Selected_3D(ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    Ve_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % E
    E_Model_Selected_3D                     = zeros(size(E_no_Delay_3D));
    E_Model_Selected_3D (ChosenByAIC_3D==5) = E_with_Delay_3D              (ChosenByAIC_3D==5);
    if useUptakeNoETM
        E_Model_Selected_3D (ChosenByAIC_3D==4) = E_with_Delay_no_Ve_3D    (ChosenByAIC_3D==4); % Ktrans is not E
    else
        E_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4); % Ktrans is not E
    end
    E_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    E_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    E_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % Ktrans
    Ktrans_Model_Selected_3D                     = zeros(size(Ktrans_no_Delay_3D));
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==5) = Ktrans_with_Delay_3D         (ChosenByAIC_3D==5);
    if useUptakeNoETM
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==4) = E_with_Delay_no_Ve_3D  (ChosenByAIC_3D==4) .* Flow_with_Delay_3D  (ChosenByAIC_3D==4);
    else
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==4) = Ktrans_with_Delay_High_F_3D  (ChosenByAIC_3D==4);
    end
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % Model Selection Map
    MeanFN=[cur_Output_directory filesep map_name];
    Raw2Nii(ChosenByAIC_3D,MeanFN,'float32',DCEFNs{1});
    %%
else % Include Delay
    %%
    NParams          = [0 1 2 2 3 3 4 4 5];
    num_maps         = length(NParams);
    size_3D          = size(RMS_Larss_no_Delay_3D);
    RMSmaps          = zeros([size_3D num_maps]);
    
    RMSmaps(:,:,:,9) = RMS_Larss_with_Delay_3D;             % 5 Params
    RMSmaps(:,:,:,8) = RMS_Larss_no_Delay_3D;               % 4 Params
    if useUptakeNoETM
        RMSmaps(:,:,:,7) = RMS_Larss_with_Delay_no_Ve_3D;      % 4 Params -----
        RMSmaps(:,:,:,6) = RMS_Larss_no_Delay_no_Ve_3D;        % 3 Params
    else
        RMSmaps(:,:,:,7) = RMS_Larss_with_Delay_High_F_3D;      % 4 Params -----
        RMSmaps(:,:,:,6) = RMS_Larss_no_Delay_High_F_3D;        % 3 Params
    end
    RMSmaps(:,:,:,5) = RMS_Larss_with_Delay_no_E_3D;        % 3 Params -----
    RMSmaps(:,:,:,4) = RMS_Larss_no_Delay_no_E_3D;          % 2 Params
    RMSmaps(:,:,:,3) = RMS_Larss_with_Delay_no_E_High_F_3D; % 2 Params -----
    RMSmaps(:,:,:,2) = RMS_Larss_no_Delay_no_E_High_F_3D;   % 1 Params
    RMSmaps(:,:,:,1) = RMS_Larss_no_Delay_zero_params_3D;   % 0 Params
    
    [ ChosenByAIC_3D ]                      = Model_Selection( nDataPoints, RMSmaps, NParams, cur_Data_Weight, AIC_Correction);
    
    % Delay
    AIF_Delay_Model_Selected_3D                    = zeros(size(est_delay_by_AIF_correct_3D));
    AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==9) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==9);
    AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==7) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==7);
    AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==5) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==5);
    AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==3) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==3);
    
    % F
    F_Model_Selected_3D                     = zeros(size(Flow_with_Delay_3D));
    F_Model_Selected_3D (ChosenByAIC_3D==9) = Flow_with_Delay_3D   (ChosenByAIC_3D==9);
    F_Model_Selected_3D (ChosenByAIC_3D==8) = Flow_no_Delay_3D     (ChosenByAIC_3D==8);
    if useUptakeNoETM
        F_Model_Selected_3D (ChosenByAIC_3D==7) = Flow_with_Delay_3D           (ChosenByAIC_3D==7); % Put zero although its Inf
        F_Model_Selected_3D (ChosenByAIC_3D==6) = Flow_no_Delay_3D             (ChosenByAIC_3D==6); % Put zero although its Inf
    else
        F_Model_Selected_3D (ChosenByAIC_3D==7) = zeros_map                    (ChosenByAIC_3D==7); % Put zero although its Inf
        F_Model_Selected_3D (ChosenByAIC_3D==6) = zeros_map                    (ChosenByAIC_3D==6); % Put zero although its Inf
    end
    F_Model_Selected_3D (ChosenByAIC_3D==5) = Flow_with_Delay_3D   (ChosenByAIC_3D==5);
    F_Model_Selected_3D (ChosenByAIC_3D==4) = Flow_no_Delay_3D     (ChosenByAIC_3D==4);
    F_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3); % Put zero although its Inf
    F_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2); % Put zero although its Inf
    F_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    
    % Vb
    Vb_Model_Selected_3D                    = zeros(size(Vb_no_Delay_3D));
    Vb_Model_Selected_3D(ChosenByAIC_3D==9) = Vb_with_Delay_3D             (ChosenByAIC_3D==9);
    Vb_Model_Selected_3D(ChosenByAIC_3D==8) = Vb_no_Delay_3D               (ChosenByAIC_3D==8);
    if useUptakeNoETM
        Vb_Model_Selected_3D(ChosenByAIC_3D==7) = Vb_with_Delay_no_Ve_3D   (ChosenByAIC_3D==7);
        Vb_Model_Selected_3D(ChosenByAIC_3D==6) = Vb_no_Delay_no_Ve_3D     (ChosenByAIC_3D==6);
    else
        Vb_Model_Selected_3D(ChosenByAIC_3D==7) = Vb_with_Delay_High_F_3D  (ChosenByAIC_3D==7);
        Vb_Model_Selected_3D(ChosenByAIC_3D==6) = Vb_no_Delay_High_F_3D    (ChosenByAIC_3D==6);
    end
    Vb_Model_Selected_3D(ChosenByAIC_3D==5) = Vb_with_Delay_no_E_3D        (ChosenByAIC_3D==5);
    Vb_Model_Selected_3D(ChosenByAIC_3D==4) = Vb_no_Delay_no_E_3D          (ChosenByAIC_3D==4);
    Vb_Model_Selected_3D(ChosenByAIC_3D==3) = Vb_with_Delay_no_E_High_F_3D (ChosenByAIC_3D==3);
    Vb_Model_Selected_3D(ChosenByAIC_3D==2) = Vb_no_Delay_no_E_High_F_3D   (ChosenByAIC_3D==2);
    Vb_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % Ve
    Ve_Model_Selected_3D                    = zeros(size(Ve_no_Delay_3D));
    Ve_Model_Selected_3D(ChosenByAIC_3D==9) = Ve_with_Delay_3D             (ChosenByAIC_3D==9);
    Ve_Model_Selected_3D(ChosenByAIC_3D==8) = Ve_no_Delay_3D               (ChosenByAIC_3D==8);
    if useUptakeNoETM
        Ve_Model_Selected_3D(ChosenByAIC_3D==7) = zeros_map      (ChosenByAIC_3D==7);
        Ve_Model_Selected_3D(ChosenByAIC_3D==6) = zeros_map        (ChosenByAIC_3D==6);
    else
        Ve_Model_Selected_3D(ChosenByAIC_3D==7) = Ve_with_Delay_High_F_3D      (ChosenByAIC_3D==7);
        Ve_Model_Selected_3D(ChosenByAIC_3D==6) = Ve_no_Delay_High_F_3D        (ChosenByAIC_3D==6);
    end
    Ve_Model_Selected_3D(ChosenByAIC_3D==5) = zeros_map                    (ChosenByAIC_3D==5);
    Ve_Model_Selected_3D(ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4);
    Ve_Model_Selected_3D(ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    Ve_Model_Selected_3D(ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    Ve_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % E
    E_Model_Selected_3D                     = zeros(size(E_no_Delay_3D));
    E_Model_Selected_3D (ChosenByAIC_3D==9) = E_with_Delay_3D              (ChosenByAIC_3D==9);
    E_Model_Selected_3D (ChosenByAIC_3D==8) = E_no_Delay_3D                (ChosenByAIC_3D==8);
    if useUptakeNoETM
        E_Model_Selected_3D (ChosenByAIC_3D==7) = E_with_Delay_no_Ve_3D    (ChosenByAIC_3D==7); % Ktrans is not E
        E_Model_Selected_3D (ChosenByAIC_3D==6) = E_no_Delay_no_Ve_3D      (ChosenByAIC_3D==6); % Ktrans is not E
    else
        E_Model_Selected_3D (ChosenByAIC_3D==7) = zeros_map                    (ChosenByAIC_3D==7); % Ktrans is not E
        E_Model_Selected_3D (ChosenByAIC_3D==6) = zeros_map                    (ChosenByAIC_3D==6); % Ktrans is not E
    end
    E_Model_Selected_3D (ChosenByAIC_3D==5) = zeros_map                    (ChosenByAIC_3D==5);
    E_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4);
    E_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    E_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    E_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    
    % Ktrans
    Ktrans_Model_Selected_3D                     = zeros(size(Ktrans_no_Delay_3D));
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==9) = Ktrans_with_Delay_3D         (ChosenByAIC_3D==9);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==8) = Ktrans_no_Delay_3D           (ChosenByAIC_3D==8);
    if useUptakeNoETM
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==7) = E_with_Delay_no_Ve_3D  (ChosenByAIC_3D==7) .*  Flow_with_Delay_3D  (ChosenByAIC_3D==7);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==6) = E_no_Delay_no_Ve_3D    (ChosenByAIC_3D==6) .*  Flow_no_Delay_3D    (ChosenByAIC_3D==6);
    else
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==7) = Ktrans_with_Delay_High_F_3D  (ChosenByAIC_3D==7);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==6) = Ktrans_no_Delay_High_F_3D    (ChosenByAIC_3D==6);
    end
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==5) = zeros_map                    (ChosenByAIC_3D==5);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);

    % Model Selection Map
    MeanFN=[cur_Output_directory filesep map_name];
    Raw2Nii(ChosenByAIC_3D,MeanFN,'float32',DCEFNs{1});
    
    
    %%
end