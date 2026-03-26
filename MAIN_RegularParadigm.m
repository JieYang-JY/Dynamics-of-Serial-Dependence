% Main code for the study of serial dependence dynamics - regular paradigm
% ----------------------------------
% This code contains several steps:
% 1) Data preparation
%    - reading data from each session file
%    - calculating RMS factor
%
% 2) Psychometric fits & Probabilistic choice model (logistic regression model)
%    - fitting psychometric function per monkey, modality, and window
%    - modeling data per monkey and modality
%    - smoothing values with Gaussian kernal
%
% 3) Plot figures

%% STEP 1: Data preparation

clear; clc; close all;
codename = 'MAIN_RegularParadigm';
codepath_name = mfilename('fullpath');
currDir = codepath_name(1 : length(codepath_name)-length(codename));
addpath(genpath(currDir))

mod = {'ves'; 'vis'};
monkey = {'B'; 'D'; 'G'; 'J'};

% ordering data files
clear FileNames_B FileNames_D FileNames_G FileNames_J
for k = 1:length(monkey) % per monkey
    data_path = strcat(currDir,'\Data\Regular paradigm\Monkey', monkey{k}, '\');
    data_file = dir(fullfile(data_path, '*_BehData.mat'));
    FileNames = {data_file.name}';

    num_run = [];
    for s = 1 : length(FileNames) % per session
        clear s_Trace;s_Trace = cell2mat(strcat(data_path, FileNames(s)));
        load(s_Trace);
        BehData.FileName

        clear index_r; index_r = find(BehData.FileName=='r');
        clear num_run_s; num_run_s = BehData.FileName(index_r+1 : end-length('.htb')); % str
        num_run(s, 1) = str2num(num_run_s);
    end

    order_sorted = []; [~, order_sorted] = sort(num_run);
    clear FileNames_ordered; FileNames_ordered = FileNames(order_sorted);

    switch k
        case 1 %monkey B
            FileNames_B = FileNames_ordered;
        case 2 %monkey D
            FileNames_D = FileNames_ordered;
        case 3 %monkey G
            FileNames_G = FileNames_ordered;
        case 4 %monkey J
            FileNames_J = FileNames_ordered;
        otherwise
            warning('Monkey Name Not Defined !')
            pause
    end
end

% concatenating sessions
clear all_data_B all_data_D all_data_G all_data_J
for k = 1:length(monkey) % per monkey
    data_path = strcat(currDir,'\Data\Regular paradigm\Monkey', monkey{k}, '\');
    clear FileNames
    switch k
        case 1 %monkey B
            FileNames = FileNames_B;
        case 2 %monkey D
            FileNames = FileNames_D;
        case 3 %monkey G
            FileNames = FileNames_G;
        case 4 %monkey J
            FileNames = FileNames_J;
    end

    clear all_data;
    all_data.allHD = [];
    all_data.currChoice = cell(length(mod), 1);
    all_data.currStim =  cell(length(mod), 1);
    all_data.prevChoice =  cell(length(mod), 1);
    all_data.prevStim =  cell(length(mod), 1);
    all_data.prior =  cell(length(mod), 1);
    for s = 1 : length(FileNames) % per session
        clear s_Trace; s_Trace = cell2mat(strcat(data_path, FileNames(s)));
        load(s_Trace);
        BehData.FileName

        all_data.allHD = [all_data.allHD; BehData.allHD];
        if BehData.StimType==1 % ves
            all_data.currChoice{1} = [all_data.currChoice{1}; BehData.choice(2:end)];
            all_data.currStim{1} = [all_data.currStim{1}; BehData.allHD(2:end)];
            all_data.prevStim{1} = [all_data.prevStim{1}; BehData.allHD(1:end-1)];
            all_data.prevChoice{1} = [all_data.prevChoice{1}; BehData.choice(1:end-1)];
            all_data.prior{1} = [all_data.prior{1}; sign(BehData.allHD(1:end-1))];
        elseif BehData.StimType==2 % vis
            all_data.currChoice{2} = [all_data.currChoice{2}; BehData.choice(2:end)];
            all_data.currStim{2} = [all_data.currStim{2}; BehData.allHD(2:end)];
            all_data.prevStim{2} = [all_data.prevStim{2}; BehData.allHD(1:end-1)];
            all_data.prevChoice{2} = [all_data.prevChoice{2}; BehData.choice(1:end-1)];
            all_data.prior{2} = [all_data.prior{2}; sign(BehData.allHD(1:end-1))];
        else
            warning('StimType Not Defined !')
            pause
        end
    end

    switch k
        case 1 %monkey B
            all_data_B = all_data;
        case 2 %monkey D
            all_data_D = all_data;
        case 3 %monkey G
            all_data_G = all_data;
        case 4 %monkey J
            all_data_J = all_data;
    end
end

%RMS factor
allHD = []; allHD = [all_data_B.allHD; all_data_D.allHD; all_data_G.allHD; all_data_J.allHD];
rms_allHD = []; rms_allHD = rms(allHD)

save_data_path = strcat(currDir,'\Data\');
save(strcat(save_data_path,'dataset_regular'), '*_B', '*_D', '*_G', '*_J', 'rms_allHD');


%% STEP 2: Psychometric fits & Probabilistic choice model (logistic regression model)

clear; clc; close all;
codename = 'MAIN_RegularParadigm';
codepath_name = mfilename('fullpath');
currDir = codepath_name(1 : length(codepath_name)-length(codename));
addpath(genpath(currDir))
addpath(genpath([currDir 'Functions\psignifit-master\private\']))

dataset = strcat(currDir,'\Data\dataset_regular.mat');
load(dataset);

LEFT = -1; RIGHT = 1; prior = [LEFT RIGHT];
mod = {'ves' 'vis'};
monkey = {'B'; 'D'; 'G'; 'J'};

PsychoFits_monks.monkey = monkey;
PsychoFits_monks.modality = mod;
PsychoFits_monks.prior = prior;
PsychoFits_monks.PsychoPerfo.pfit_output = cell(length(monkey), length(mod)); % Psychometric fits: performance
PsychoFits_monks.PsychoPerfo.pseudoR2 = cell(length(monkey), length(mod));
PsychoFits_monks.PsychoPerfo.PSE = cell(length(monkey), length(mod));
PsychoFits_monks.PsychoPerfo.thre = cell(length(monkey), length(mod));
PsychoFits_monks.SerialDepend.pfit_output_pHD = cell(length(monkey), length(mod)); % Psychometric fits: aggregate serial dependence, sorted by prev heading
PsychoFits_monks.SerialDepend.pseudoR2_pHD = cell(length(monkey), length(mod));
PsychoFits_monks.SerialDepend.PSE_pHD = cell(length(monkey), length(mod));
PsychoFits_monks.SerialDepend.deltaPSE_pHD = cell(length(monkey), length(mod));
PsychoFits_monks.SerialDepend.pfit_output_pCho = cell(length(monkey), length(mod)); % Psychometric fits: aggregate serial dependence, sorted by prev choice
PsychoFits_monks.SerialDepend.pseudoR2_pCho = cell(length(monkey), length(mod));
PsychoFits_monks.SerialDepend.PSE_pCho = cell(length(monkey), length(mod));
PsychoFits_monks.SerialDepend.deltaPSE_pCho = cell(length(monkey), length(mod));
LogRegModel_monks.model = cell(length(monkey), length(mod)); % Modeling
LogRegModel_monks.Beta_currStim = cell(length(monkey), length(mod));
LogRegModel_monks.Beta_prevStim = cell(length(monkey), length(mod));
LogRegModel_monks.Beta_prevCho = cell(length(monkey), length(mod));
LogRegModel_monks.deltaBeta = cell(length(monkey), length(mod));
LogRegModel_monks.sm_Beta_currStim = cell(length(monkey), length(mod));
LogRegModel_monks.sm_Beta_prevStim = cell(length(monkey), length(mod));
LogRegModel_monks.sm_Beta_prevCho = cell(length(monkey), length(mod));
LogRegModel_monks.sm_deltaBeta = cell(length(monkey), length(mod));
PsychoFits_monks.avg_pseudoR2 = []; 
PsychoFits_monks.sd_pseudoR2 = []; 

window_size = 1000;
step_size = 500;
for k = 1:length(monkey) % per monkey
    disp(['Analyzing the data of monkey ' monkey{k} ' .'])
    clear all_data
    switch k
        case 1 %monkey B
            all_data = all_data_B;
        case 2 %monkey D
            all_data = all_data_D;
        case 3 %monkey G
            all_data = all_data_G;
        case 4 %monkey J
            all_data = all_data_J;
    end

    for m = 1:length(mod) % per modality
        disp(['Sliding in ' cell2mat(mod(m)) ' ...'])
        pos_windowEnd = window_size; prev_pos_windowEnd = 0; winInd = 0; % RESET
        pfit_output_Pf = cell(0); threshold_Pf = []; PSE_Pf = []; pseudoR2_Pf = [];
        pfit_output_SD_pHD = cell(0); PSE_SD_pHD = []; pseudoR2_SD_pHD = []; deltaPSE_pHD = [];
        pfit_output_SD_pCho = cell(0); PSE_SD_pCho = []; pseudoR2_SD_pCho = []; deltaPSE_pCho = [];
        mdl_glmfit = cell(0); Beta_currStim = []; Beta_prevStim = []; Beta_prevCho = []; deltaBeta = [];
        while prev_pos_windowEnd < length(all_data.currStim{m}) % until the last trial of given stim type
            winInd = winInd+1

            % Psychometric fits
            clear temp_priorHD_type temp_priorCho_type temp_choice temp_anaHD unique_heading;
            temp_priorHD_type = all_data.prior{m}(pos_windowEnd-window_size+1 : pos_windowEnd);
            temp_priorCho_type = all_data.prevChoice{m}(pos_windowEnd-window_size+1 : pos_windowEnd);
            temp_choice = all_data.currChoice{m}(pos_windowEnd-window_size+1 : pos_windowEnd);
            temp_anaHD = all_data.currStim{m}(pos_windowEnd-window_size+1 : pos_windowEnd);
            unique_heading = unique(temp_anaHD);
            %performance ----------------
            psycho_right_Pf = []; fit_data_psycho_cum_Pf = [];
            for h = 1:length(unique_heading) % per unique heading
                clear temp_heading; temp_heading = logical(temp_anaHD==unique_heading(h));
                clear temp_right_choice_trials;temp_right_choice_trials = (temp_heading & temp_choice==RIGHT);
                psycho_right_Pf(h) = 1*sum(temp_right_choice_trials) / sum(temp_heading);
                fit_data_psycho_cum_Pf(h, 1) = unique_heading(h);
                fit_data_psycho_cum_Pf(h, 2) = psycho_right_Pf(h);
                fit_data_psycho_cum_Pf(h, 3) = sum(temp_heading);
            end
            pfit_output_Pf{winInd, 1} = getPsychoFit(fit_data_psycho_cum_Pf(:,:), fit_data_psycho_cum_Pf(:,1));
            threshold_Pf(winInd, 1) = pfit_output_Pf{winInd, 1}.thresh;
            PSE_Pf(winInd, 1) = pfit_output_Pf{winInd, 1}.bias;
            pseudoR2_Pf(winInd, 1) = pfit_output_Pf{winInd, 1}.pseudoR2;
            %aggregate serial dependence ----------------
            psycho_right_SD_pHD = []; fit_data_psycho_cum_SD_pHD = cell(1,length(prior));
            psycho_right_SD_pCho = []; fit_data_psycho_cum_SD_pCho = cell(1,length(prior));
            for p = 1:length(prior) % per prior
                for h = 1:length(unique_heading) % per unique heading
                    clear temp_heading_pHD; temp_heading_pHD = logical(temp_priorHD_type==prior(p) & temp_anaHD==unique_heading(h));
                    clear temp_right_choice_trials_pHD; temp_right_choice_trials_pHD = (temp_heading_pHD & temp_choice==RIGHT);
                    psycho_right_SD_pHD(p,h) = 1*sum(temp_right_choice_trials_pHD) / sum(temp_heading_pHD);
                    fit_data_psycho_cum_SD_pHD{p}(h, 1) = unique_heading(h);
                    fit_data_psycho_cum_SD_pHD{p}(h, 2) = psycho_right_SD_pHD(p,h);
                    fit_data_psycho_cum_SD_pHD{p}(h, 3) = sum(temp_heading_pHD);
                    % ------------------
                    clear temp_heading_pCho; temp_heading_pCho = logical(temp_priorCho_type==prior(p) & temp_anaHD==unique_heading(h));
                    clear temp_right_choice_trials_pCho; temp_right_choice_trials_pCho = (temp_heading_pCho & temp_choice==RIGHT);
                    psycho_right_SD_pCho(p,h) = 1*sum(temp_right_choice_trials_pCho) / sum(temp_heading_pCho);
                    fit_data_psycho_cum_SD_pCho{p}(h, 1) = unique_heading(h);
                    fit_data_psycho_cum_SD_pCho{p}(h, 2) = psycho_right_SD_pCho(p,h);
                    fit_data_psycho_cum_SD_pCho{p}(h, 3) = sum(temp_heading_pCho);
                end
                clear temp_index;temp_index = find(fit_data_psycho_cum_SD_pHD{p}(:,3)==0);
                if temp_index, fit_data_psycho_cum_SD_pHD{p}(temp_index,:) = []; end
                pfit_output_SD_pHD{winInd, p} = getPsychoFit(fit_data_psycho_cum_SD_pHD{p}(:,:), fit_data_psycho_cum_SD_pHD{p}(:,1));
                PSE_SD_pHD(winInd, p) = pfit_output_SD_pHD{winInd, p}.bias;
                pseudoR2_SD_pHD(winInd, p) = pfit_output_SD_pHD{winInd, p}.pseudoR2;
                % ------------------
                clear temp_index;temp_index = find(fit_data_psycho_cum_SD_pCho{p}(:,3)==0);
                if temp_index, fit_data_psycho_cum_SD_pCho{p}(temp_index,:) = []; end
                pfit_output_SD_pCho{winInd, p} = getPsychoFit(fit_data_psycho_cum_SD_pCho{p}(:,:), fit_data_psycho_cum_SD_pCho{p}(:,1));
                PSE_SD_pCho(winInd, p) = pfit_output_SD_pCho{winInd, p}.bias;
                pseudoR2_SD_pCho(winInd, p) = pfit_output_SD_pCho{winInd, p}.pseudoR2;
            end
            deltaPSE_pHD(winInd, 1) = PSE_SD_pHD(winInd, 1) - PSE_SD_pHD(winInd, 2);
            % ------------------
            deltaPSE_pCho(winInd, 1) = PSE_SD_pCho(winInd, 1) - PSE_SD_pCho(winInd, 2);

            % Modeling
            clear temp_currStim temp_prevStim temp_currChoice temp_prevChoice;
            temp_currChoice = temp_choice; temp_currChoice(temp_currChoice==-1) = 0; % probability range: [0, 1]
            temp_currStim = temp_anaHD ./ rms_allHD; % normalizing
            temp_prevStim = all_data.prevStim{m}(pos_windowEnd-window_size+1 : pos_windowEnd) ./ rms_allHD; % normalizing
            temp_prevChoice = all_data.prevChoice{m}(pos_windowEnd-window_size+1 : pos_windowEnd);
            mdl_glmfit{winInd, 1} = fitglm([temp_currStim temp_prevStim temp_prevChoice],temp_currChoice,'Distribution','binomial','Link','logit');
            Beta_currStim(winInd, 1) = table2array(mdl_glmfit{winInd, 1}.Coefficients(2, 1));
            Beta_prevStim(winInd, 1) = table2array(mdl_glmfit{winInd, 1}.Coefficients(3, 1));
            Beta_prevCho(winInd, 1) = table2array(mdl_glmfit{winInd, 1}.Coefficients(4, 1));
            deltaBeta(winInd, 1) = Beta_prevStim(end) - Beta_prevCho(end);

            % Sliding forward
            prev_pos_windowEnd = pos_windowEnd;
            pos_windowEnd = pos_windowEnd + step_size;
            if pos_windowEnd > length(all_data.currStim{m})
                pos_windowEnd = length(all_data.currStim{m});
            end
        end
        
        % Psychometric fits: performance
        PsychoFits_monks.PsychoPerfo.pfit_output{k, m} = pfit_output_Pf;
        PsychoFits_monks.PsychoPerfo.pseudoR2{k, m} = pseudoR2_Pf;
        PsychoFits_monks.PsychoPerfo.PSE{k, m} = PSE_Pf;
        PsychoFits_monks.PsychoPerfo.thre{k, m} = threshold_Pf;
        % Psychometric fits: aggregate serial dependence, sorted by previous heading
        PsychoFits_monks.SerialDepend.pfit_output_pHD{k, m} = pfit_output_SD_pHD; 
        PsychoFits_monks.SerialDepend.pseudoR2_pHD{k, m} = pseudoR2_SD_pHD;
        PsychoFits_monks.SerialDepend.PSE_pHD{k, m} = PSE_SD_pHD;
        PsychoFits_monks.SerialDepend.deltaPSE_pHD{k, m} = deltaPSE_pHD;
        % Psychometric fits: aggregate serial dependence, sorted by previous choice
        PsychoFits_monks.SerialDepend.pfit_output_pCho{k, m} = pfit_output_SD_pCho; 
        PsychoFits_monks.SerialDepend.pseudoR2_pCho{k, m} = pseudoR2_SD_pCho;
        PsychoFits_monks.SerialDepend.PSE_pCho{k, m} = PSE_SD_pCho;
        PsychoFits_monks.SerialDepend.deltaPSE_pCho{k, m} = deltaPSE_pCho;
        % Modeling
        LogRegModel_monks.model{k, m} = mdl_glmfit; 
        LogRegModel_monks.Beta_currStim{k, m} = Beta_currStim;
        LogRegModel_monks.Beta_prevStim{k, m} = Beta_prevStim;
        LogRegModel_monks.Beta_prevCho{k, m} = Beta_prevCho;
        LogRegModel_monks.deltaBeta{k, m} = deltaBeta;

        % Smooth modeling results
        sigma = 6;
        LogRegModel_monks.sm_Beta_currStim{k, m} = getGaussResults(Beta_currStim, sigma);
        LogRegModel_monks.sm_Beta_prevStim{k, m} = getGaussResults(Beta_prevStim, sigma);
        LogRegModel_monks.sm_Beta_prevCho{k, m} = getGaussResults(Beta_prevCho, sigma);
        LogRegModel_monks.sm_deltaBeta{k, m} = getGaussResults(deltaBeta, sigma);
        
        % Spearman's rank correlation
        [LogRegModel_monks.corrSpm_r_currStim{k, m}, LogRegModel_monks.corrSpm_pValue_currStim{k, m}] = ...
            corr(LogRegModel_monks.sm_Beta_currStim{k, m},[1:length(LogRegModel_monks.sm_Beta_currStim{k, m})]','Type','Spearman');
        [LogRegModel_monks.corrSpm_r_prevStim{k, m}, LogRegModel_monks.corrSpm_pValue_prevStim{k, m}] = ...
            corr(LogRegModel_monks.sm_Beta_prevStim{k, m},[1:length(LogRegModel_monks.sm_Beta_prevStim{k, m})]','Type','Spearman');
        [LogRegModel_monks.corrSpm_r_prevCho{k, m}, LogRegModel_monks.corrSpm_pValue_prevCho{k, m}] = ...
            corr(LogRegModel_monks.sm_Beta_prevCho{k, m},[1:length(LogRegModel_monks.sm_Beta_prevCho{k, m})]','Type','Spearman');
        [LogRegModel_monks.corrSpm_r_deltaBeta{k, m}, LogRegModel_monks.corrSpm_pValue_deltaBeta{k, m}] = ...
            corr(LogRegModel_monks.sm_deltaBeta{k, m},[1:length(LogRegModel_monks.sm_deltaBeta{k, m})]','Type','Spearman');
    end
end
%%%%%%%%%%%%%%%%
temp_pseudoR2 = [];
temp_pse_start5 = []; temp_pse_end5 = [];
temp_thre_start5 = []; temp_thre_end5 = [];
temp_dpse_pHD_start3rd = []; temp_dpse_pHD_end3rd = [];
temp_dpse_pCho_start3rd = []; temp_dpse_pCho_end3rd = [];
for m = 1:length(PsychoFits_monks.modality)
    for k = 1:length(PsychoFits_monks.monkey)
        temp_pseudoR2 = [temp_pseudoR2; PsychoFits_monks.PsychoPerfo.pseudoR2{k,m}]; %%%
        temp_pse_start5 = [temp_pse_start5; PsychoFits_monks.PsychoPerfo.PSE{k,m}(1:5)]; %%%
        temp_pse_end5 = [temp_pse_end5; PsychoFits_monks.PsychoPerfo.PSE{k,m}(end-4:end)];
        temp_thre_start5 = [temp_thre_start5; PsychoFits_monks.PsychoPerfo.thre{k,m}(1:5)]; %%%
        temp_thre_end5 = [temp_thre_end5; PsychoFits_monks.PsychoPerfo.thre{k,m}(end-4:end)];
        temp_len = floor(length(PsychoFits_monks.SerialDepend.deltaPSE_pHD{k, m})/3);
        temp_dpse_pHD_start3rd = [temp_dpse_pHD_start3rd; PsychoFits_monks.SerialDepend.deltaPSE_pHD{k, m}(1:temp_len)]; %%%
        temp_dpse_pHD_end3rd = [temp_dpse_pHD_end3rd; PsychoFits_monks.SerialDepend.deltaPSE_pHD{k, m}(end-temp_len+1:end)];
        temp_dpse_pCho_start3rd = [temp_dpse_pCho_start3rd; PsychoFits_monks.SerialDepend.deltaPSE_pCho{k, m}(1:temp_len)]; %%%
        temp_dpse_pCho_end3rd = [temp_dpse_pCho_end3rd; PsychoFits_monks.SerialDepend.deltaPSE_pCho{k, m}(end-temp_len+1:end)];
        for p = 1:length(PsychoFits_monks.prior)
            temp_pseudoR2 = [temp_pseudoR2; PsychoFits_monks.SerialDepend.pseudoR2_pHD{k,m}(:,p); ...
                PsychoFits_monks.SerialDepend.pseudoR2_pCho{k,m}(:,p)];
        end 
    end
end
PsychoFits_monks.avg_pseudoR2 = mean(temp_pseudoR2); 
PsychoFits_monks.sd_pseudoR2 = std(temp_pseudoR2);
fprintf('pseudoR2 (all): mean = %.2f, sd = %.3f, min = %.2f\n\n', mean(temp_pseudoR2), std(temp_pseudoR2), min(temp_pseudoR2))
%--------------------------------------------------------------------------
fprintf('PSE (start): mean = %.2f, sd = %.2f\n', mean(temp_pse_start5), std(temp_pse_start5))
fprintf('PSE (end): mean = %.3f, sd = %.2f\n\n', mean(temp_pse_end5), std(temp_pse_end5))
%--------------------------------------------------------------------------
fprintf('Threshold (start): mean = %.2f, sd = %.2f\n', mean(temp_thre_start5), std(temp_thre_start5))
fprintf('Threshold (end): mean = %.2f, sd = %.2f\n\n', mean(temp_thre_end5), std(temp_thre_end5))
%--------------------------------------------------------------------------
fprintf('\x0394PSE (sorted by prev heading, start): mean = %.2f, sd = %.2f\n', mean(temp_dpse_pHD_start3rd), std(temp_dpse_pHD_start3rd))
fprintf('\x0394PSE (sorted by prev heading, end): mean = %.2f, sd = %.2f\n\n', mean(temp_dpse_pHD_end3rd), std(temp_dpse_pHD_end3rd))
%--------------------------------------------------------------------------
fprintf('\x0394PSE (sorted by prev choice, start): mean = %.2f, sd = %.2f\n', mean(temp_dpse_pCho_start3rd), std(temp_dpse_pCho_start3rd))
fprintf('\x0394PSE (sorted by prev choice, end): mean = %.2f, sd = %.2f\n\n', mean(temp_dpse_pCho_end3rd), std(temp_dpse_pCho_end3rd))


save_results_path = strcat(currDir,'\Results\Regular paradigm\');
save(strcat(save_results_path,'PsychoFits_LogRegModel_regular'), 'PsychoFits_monks', 'LogRegModel_monks');


%% STEP 3: Plot figures

clear; clc; close all;
codename = 'MAIN_RegularParadigm';
codepath_name = mfilename('fullpath');
currDir = codepath_name(1 : length(codepath_name)-length(codename));
addpath(genpath(currDir))

path_PsychoFits_LogRegModel = strcat(currDir,'\Results\Regular paradigm\PsychoFits_LogRegModel_regular.mat');
load(path_PsychoFits_LogRegModel);
prior = PsychoFits_monks.prior % -1: LEFT; 1: RIGHT
mod = PsychoFits_monks.modality
monkey = PsychoFits_monks.monkey
length_filename = length('PsychoFits_LogRegModel_regular.mat');
save_fig_path = path_PsychoFits_LogRegModel(1:(end-length_filename));
color{1} = [0 0 1]; % blue
color{2} = [1 0 0]; % red

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ploting: Figure 2B, 2C & Supplementary figure 4')
FigureIndex1 = 1; figure(FigureIndex1);
set(FigureIndex1,'Position', [10,80 1400,600], 'Name', 'Figure 2B');
sub1 = cell(length(monkey), length(mod));
FigureIndex2 = 2; figure(FigureIndex2);
set(FigureIndex2,'Position', [10,80 1400,600], 'Name', 'Figure 2C');
sub2 = cell(length(monkey), length(mod));
FigureIndex3 = 3; figure(FigureIndex3);
set(FigureIndex3,'Position', [10,80 1400,600], 'Name', 'Supplementary figure 4');
sub3 = cell(length(monkey), length(mod));
for k = 1 : length(monkey)
    for m = 1 : length(mod)
        % beta_currStim ---------------------------
        figure(FigureIndex1)
        sub1{k,m} = subplot(length(mod), length(monkey), k+4*(m-1));
        plot(LogRegModel_monks.sm_Beta_currStim{k, m}, 's', 'MarkerFaceColor', color{m}, 'MarkerEdgeColor', color{m})
        hold on; ylim([0, 14.5]); 
        temp_winlen = length(LogRegModel_monks.sm_Beta_currStim{k, m});
        xlim([0-round(temp_winlen*0.01), temp_winlen+floor(temp_winlen*0.05)]);
        if m==1, title(['Monkey ',monkey{k}]); end
        set(gca,'YTick',[0 7 14]); set(gca,'YTickLabel',{'0', '7', '14'});
        clear legend_txt, legend_txt = '\beta_{currStim}';
        if k==1
            clear leg, leg = legend(legend_txt, 'Location', 'northwest'); legend('boxoff');
            temp_pos = leg.Position;
            leg.Position = [temp_pos(1)-0.05*sub1{k,m}.Position(3) temp_pos(2)-0.1*sub1{k,m}.Position(4) temp_pos(3:4)];
            temp_xlim = []; temp_xlim = xlim;
            if m==1
                txt = 'ves'; temp_text_x1 = temp_pos(1)*0.7; text_y1 = temp_pos(2)*1.02*max(ylim);
            else
                txt = 'vis';
            end
            text(temp_text_x1*temp_xlim(2), text_y1, txt, 'Color', color{m}, 'FontWeight', 'bold', 'FontSize', 11);
        end
        % beta_prev ---------------------------
        figure(FigureIndex2)
        sub2{k,m} = subplot(length(mod), length(monkey), k+4*(m-1));
        if m==1, title(['Monkey ',monkey{k}]); end
        plot(LogRegModel_monks.sm_Beta_prevStim{k, m}, 's', 'MarkerFaceColor', color{m}, 'MarkerEdgeColor', color{m}); hold on; 
        plot(LogRegModel_monks.sm_Beta_prevCho{k, m}, 's', 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', [0.6 0.6 0.6]); hold on;
        legend_txt = cell(1, 2); legend_txt{1} = '\beta_{prevStim}'; legend_txt{2} = '\beta_{prevCho}'; 
        ylim([-1.3, 0.5]); set(gca,'YTick',[-1.2 -0.6 0]); set(gca,'YTickLabel',{'-1.2', '-0.6', '0'});
        xlim([0-round(temp_winlen*0.01), temp_winlen+floor(temp_winlen*0.05)]);
        if m==1, title(['Monkey ',monkey{k}]); end
        if k==1
            clear leg, leg = legend(legend_txt, 'Location', 'southwest'); legend('boxoff');
            temp_pos = leg.Position;
            leg.Position = [temp_pos(1)-0.05*sub2{k,m}.Position(3) temp_pos(2:4)];
            if m==1
                txt = 'ves'; text_y2 = (1-temp_pos(2))*0.85*diff(ylim)+min(ylim); 
            else
                txt = 'vis'; 
            end
            text(temp_text_x1*temp_xlim(2), text_y2, txt, 'Color', color{m}, 'FontWeight', 'bold', 'FontSize', 11);
        end
        plot([min(xlim), max(xlim)], [0 0], '--k'); hold on;
        % deltabeta ---------------------------
        figure(FigureIndex3)
        sub3{k,m} = subplot(length(mod), length(monkey), k+4*(m-1));
        if m==1, title(['Monkey ',monkey{k}]); end
        plot(LogRegModel_monks.sm_deltaBeta{k, m}, 's', 'MarkerFaceColor', color{m}, 'MarkerEdgeColor', color{m}); hold on; 
        ylim([-1.45, 0.7]); set(gca,'YTick',[-1.4 -0.7 0 0.7]); set(gca,'YTickLabel',{'-1.4', '-0.7', '0', '0.7'});
        xlim([0-round(temp_winlen*0.01), temp_winlen+floor(temp_winlen*0.05)]);
        if m==1, title(['Monkey ',monkey{k}]); end
        if k==1
            if m==1
                txt = 'ves'; text_y3 = 0.15*diff(ylim)+min(ylim); 
            else
                txt = 'vis'; 
            end
            text(temp_text_x1*temp_xlim(2), text_y3, txt, 'Color', color{m}, 'FontWeight', 'bold', 'FontSize', 11);
        end
        plot([min(xlim), max(xlim)], [0 0], '--k'); hold on;
    end
end
figure(FigureIndex1)
axes('position',[0.09 0.475, 0.1 0.1]), axis off;
text(0, 0, '\beta [a.u.]', 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 90);
axes('position',[0.45 0.05, 0.1 0.1]), axis off;
text(0, 0, 'Training process (Window #)', 'FontSize', 12, 'FontWeight', 'bold');
saveas(gcf,[save_fig_path,'Figure 2B.png'],'png');
close
figure(FigureIndex2)
axes('position',[0.09 0.475, 0.1 0.1]), axis off;
text(0, 0, '\beta [a.u.]', 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 90);
axes('position',[0.45 0.05, 0.1 0.1]), axis off;
text(0, 0, 'Training process (Window #)', 'FontSize', 12, 'FontWeight', 'bold');
saveas(gcf,[save_fig_path,'Figure 2C.png'],'png');
close
figure(FigureIndex3)
axes('position',[0.09 0.345, 0.1 0.1]), axis off;
text(0, 0, '\Delta\beta (\beta_{prevStim} - \beta_{prevCho}) [a.u.]', 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 90);
axes('position',[0.45 0.05, 0.1 0.1]), axis off;
text(0, 0, 'Training process (Window #)', 'FontSize', 12, 'FontWeight', 'bold');
saveas(gcf,[save_fig_path,'Supplementary figure 4.png'],'png');
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ploting: Figure 1D') % Example psychometric performance
FigureIndex1 = 1; figure(FigureIndex1);
set(FigureIndex1, 'Name', 'Figure 1D');
monkIdx = 3; % Monkey G
winInd = [92, 60]; % the exmple window for each mod
legend_txt = cell(length(mod)*2, 1);
index = []; % illustrating PSE
for m = 1 : length(mod)
    x = []; y = [];
    x = PsychoFits_monks.PsychoPerfo.pfit_output{monkIdx, m}{winInd(m)}.input(1:end, 1);
    y = PsychoFits_monks.PsychoPerfo.pfit_output{monkIdx, m}{winInd(m)}.input(1:end, 2) ./ ...
        PsychoFits_monks.PsychoPerfo.pfit_output{monkIdx, m}{winInd(m)}.input(1:end, 3);
    plot(x, y, 'o', 'Color', color{m}, 'LineWidth', 2, 'MarkerSize', 8); hold on;
    xi = []; pfitcurve = [];
    xi = PsychoFits_monks.PsychoPerfo.pfit_output{monkIdx, m}{winInd(m)}.xi;
    pfitcurve = PsychoFits_monks.PsychoPerfo.pfit_output{monkIdx, m}{winInd(m)}.pfitcurve;
    plot(xi, pfitcurve, '-', 'Color', color{m}, 'LineWidth', 2, 'MarkerSize', 8); hold on;
    if m==1, legend_txt{m*2-1} = 'ves'; else, legend_txt{m*2-1} = 'vis'; end
    legend_txt{m*2} = [''];
    % -------------------------------
    diff = []; diff = abs(pfitcurve - 0.5); 
    index(m) = find(min(diff)==diff);
    if isempty(index(m)); index(m) = find(-min(diff)==diff); end
end
legend(legend_txt, 'Location', 'northwest'); 
for m = 1 : length(mod)
    plot([xi(index(m)), xi(index(m))], [0 0.5], '--', 'Color', color{m}, 'LineWidth', 1); hold on;
end
xlim([min(x),max(x)]); ylim([0, 1]); % configuration
plot([min(xlim), max(xlim)], [0.5 0.5], '--k'); hold on;
text(2, 0.05, 'PSE', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Prop. Rightward Choices'); xlabel('Heading [deg]');
set(gca,'XTick',min(xi):0.5*max(xi):max(xi));
set(gca,'XTickLabel',{num2str(min(xi)),num2str(0.5*min(xi)),'0',num2str(0.5*max(xi)),num2str(max(xi))});
set(gca,'YTick',0:0.5:1);
saveas(gcf,[save_fig_path,'Figure 1D.png'],'png');
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ploting: Supplementary figure 3A') % Serial-dependence example of psycho fits
FigureIndex1 = 1; f = figure(FigureIndex1); f.Position(3) = [1100];
set(FigureIndex1, 'Name', 'Supplementary figure 3A');
monkIdx = 3; % Monkey G
winInd = [92, 60]; % the exmple window for each mod
prior_type = {'left prior'; 'right prior'};
for m = 1 : length(mod)
    ax{m} = axes('Position', [0.07+0.5*(m-1) 0.11 0.4 0.8]);
    legend_txt = cell(length(prior), 1); 
    for p = 1 : length(prior)
        x = []; y = []; % points
        x = PsychoFits_monks.SerialDepend.pfit_output_pHD{monkIdx, m}{winInd(m), p}.input(1:end, 1);
        y = PsychoFits_monks.SerialDepend.pfit_output_pHD{monkIdx, m}{winInd(m), p}.input(1:end, 2) ./ ...
            PsychoFits_monks.SerialDepend.pfit_output_pHD{monkIdx, m}{winInd(m), p}.input(1:end, 3);
        if p==1
            temp_color = color{m};
        else
            temp_color = color{m} + 0.5; temp_color(temp_color>1) = 1;
        end
        plot(x, y, 'o', 'Color', temp_color, 'LineWidth', 2, 'MarkerSize', 8); hold on;
        xi = []; pfitcurve = []; % curve
        xi = PsychoFits_monks.SerialDepend.pfit_output_pHD{monkIdx, m}{winInd(m), p}.xi;
        pfitcurve = PsychoFits_monks.SerialDepend.pfit_output_pHD{monkIdx, m}{winInd(m), p}.pfitcurve;
        plot(xi, pfitcurve, '-', 'Color', temp_color, 'LineWidth', 2, 'MarkerSize', 8); hold on;
        legend_txt{p*2-1} = prior_type{p}; legend_txt{p*2} = [''];
    end
    legend(legend_txt, 'Location', 'southeast'); legend('boxoff');
    xlim([min(xi),max(xi)]); ylim([0,1]);
    ylabel('Prop. Rightward Choices'); xlabel('Heading [deg]');
    set(gca,'XTick',min(xi):0.5*max(xi):max(xi));
    set(gca,'XTickLabel',{num2str(min(xi)),num2str(0.5*min(xi)),'0',num2str(0.5*max(xi)),num2str(max(xi))});
    set(gca,'YTick',0:0.5:1);
    plot([0 0], [0 1], '--k'); hold on; plot([min(xi) max(xi)], [0.5 0.5], '--k'); hold on;
    if m==1, txt = 'ves'; else, txt = 'vis'; end
    text(8.75, 0.29, txt, 'FontWeight', 'bold', 'Color', color{m}, 'FontSize', 11);
end
saveas(gcf,[save_fig_path,'Supplementary figure 3A.png'],'png');
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ploting: Supplementary figure 4A') % Serial-dependence example of psycho fits
FigureIndex1 = 1; f = figure(FigureIndex1); f.Position(3) = [1100];
set(FigureIndex1, 'Name', 'Supplementary figure 3A');
monkIdx = 3; % Monkey G
winInd = [92, 60]; % the exmple window for each mod
prior_type = {'left prior'; 'right prior'};
for m = 1 : length(mod)
    ax{m} = axes('Position', [0.07+0.5*(m-1) 0.11 0.4 0.8]);
    legend_txt = cell(length(prior), 1); 
    for p = 1 : length(prior)
        x = []; y = []; % points
        x = PsychoFits_monks.SerialDepend.pfit_output_pCho{monkIdx, m}{winInd(m), p}.input(1:end, 1);
        y = PsychoFits_monks.SerialDepend.pfit_output_pCho{monkIdx, m}{winInd(m), p}.input(1:end, 2) ./ ...
            PsychoFits_monks.SerialDepend.pfit_output_pCho{monkIdx, m}{winInd(m), p}.input(1:end, 3);
        if p==1
            temp_color = color{m};
        else
            temp_color = color{m} + 0.5; temp_color(temp_color>1) = 1;
        end
        plot(x, y, 'o', 'Color', temp_color, 'LineWidth', 2, 'MarkerSize', 8); hold on;
        xi = []; pfitcurve = []; % curve
        xi = PsychoFits_monks.SerialDepend.pfit_output_pCho{monkIdx, m}{winInd(m), p}.xi;
        pfitcurve = PsychoFits_monks.SerialDepend.pfit_output_pCho{monkIdx, m}{winInd(m), p}.pfitcurve;
        plot(xi, pfitcurve, '-', 'Color', temp_color, 'LineWidth', 2, 'MarkerSize', 8); hold on;
        legend_txt{p*2-1} = prior_type{p}; legend_txt{p*2} = [''];
    end
    legend(legend_txt, 'Location', 'southeast'); legend('boxoff');
    xlim([min(xi),max(xi)]); ylim([0,1]);
    ylabel('Prop. Rightward Choices'); xlabel('Heading [deg]');
    set(gca,'XTick',min(xi):0.5*max(xi):max(xi));
    set(gca,'XTickLabel',{num2str(min(xi)),num2str(0.5*min(xi)),'0',num2str(0.5*max(xi)),num2str(max(xi))});
    set(gca,'YTick',0:0.5:1);
    plot([0 0], [0 1], '--k'); hold on; plot([min(xi) max(xi)], [0.5 0.5], '--k'); hold on;
    if m==1, txt = 'ves'; else, txt = 'vis'; end
    text(8.75, 0.29, txt, 'FontWeight', 'bold', 'Color', color{m}, 'FontSize', 11);
end
saveas(gcf,[save_fig_path,'Supplementary figure 4A.png'],'png');
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ploting: Supplementary figure 2A, 2B, 3B, 4B') % the dynamics of a series of psychometric parameters
FigureIndex1 = 1; figure(FigureIndex1);
set(FigureIndex1,'Position', [10,80 1400,600], 'Name', 'Supplementary figure 2A');
FigureIndex2 = 2; figure(FigureIndex2);
set(FigureIndex2,'Position', [10,80 1400,600], 'Name', 'Supplementary figure 2B');
FigureIndex3 = 3; figure(FigureIndex3);
set(FigureIndex3,'Position', [10,80 1400,600], 'Name', 'Supplementary figure 3B');
FigureIndex4 = 4; figure(FigureIndex4);
set(FigureIndex4,'Position', [10,80 1400,600], 'Name', 'Supplementary figure 4B');
for k = 1 : length(monkey)
    for m = 1 : length(mod)
        % PSE ---------------------------
        figure(FigureIndex1)
        subplot(length(mod), length(monkey), k+4*(m-1));
        plot(PsychoFits_monks.PsychoPerfo.PSE{k, m}, 'o', 'Color', color{m})
        hold on; ylim([-50, 50]);
        temp_winlen = length(PsychoFits_monks.PsychoPerfo.PSE{k, m});
        xlim([0-round(temp_winlen*0.01), temp_winlen+floor(temp_winlen*0.05)]);
        if m==1, title(['Monkey ',monkey{k}]); end
        plot([min(xlim), max(xlim)], [0 0], '--k'); hold on;
        set(gca,'YTick',[-40 0 40]); set(gca,'YTickLabel',{'-40', '0', '40'});
        if k==1
            if m==1, txt = 'ves'; else, txt = 'vis'; end
            text(max(xlim)*0.75, 35, txt, 'FontWeight', 'bold', 'Color', color{m}, 'FontSize', 11);
        end
        % Threshold ---------------------------
        figure(FigureIndex2)
        subplot(length(mod), length(monkey), k+4*(m-1));
        log2_thre = log2(PsychoFits_monks.PsychoPerfo.thre{k, m});
        plot(log2_thre, 'o', 'Color', color{m})
        hold on; ylim([0, 7]);
        xlim([0-round(temp_winlen*0.01), temp_winlen+floor(temp_winlen*0.05)]);
        if m==1, title(['Monkey ',monkey{k}]); end
        plot([min(xlim), max(xlim)], [log2(6) log2(6)], '--k'); hold on;
        set(gca,'YTick',[1 3 5 7]); set(gca,'YTickLabel',{'2', '8', '32', '128'});
        if k==1
            if m==1, txt = 'ves'; else, txt = 'vis'; end
            text(max(xlim)*0.75, 6, txt, 'FontWeight', 'bold', 'Color', color{m}, 'FontSize', 11);
        end
        % deltaPSE, sorted by prev heading ---------------------------
        figure(FigureIndex3)
        subplot(length(mod), length(monkey), k+4*(m-1));
        plot(PsychoFits_monks.SerialDepend.deltaPSE_pHD{k, m}, 'o', 'Color', color{m}); hold on;
        if m==1
            title(['Monkey ',monkey{k}]);
            ylim([-75, 90]); set(gca,'YTick',[-50 0 50]); set(gca,'YTickLabel',{'-50', '0', '50'});
        else
            ylim([-80, 20]); set(gca,'YTick',[-60 -30 0]); set(gca,'YTickLabel',{'-60', '-30', '0'});
        end
        xlim([0-round(temp_winlen*0.01), temp_winlen+floor(temp_winlen*0.05)]);
        plot([min(xlim), max(xlim)], [0 0], '--k'); hold on;
        if k==1
            if m==1, txt = 'ves'; temp_y = -55; else, txt = 'vis'; temp_y = -65; end
            text(max(xlim)*0.75, temp_y, txt, 'FontWeight', 'bold', 'Color', color{m}, 'FontSize', 11);
        end
        % deltaPSE, sorted by prev choice ---------------------------
        figure(FigureIndex4)
        subplot(length(mod), length(monkey), k+4*(m-1));
        plot(PsychoFits_monks.SerialDepend.deltaPSE_pCho{k, m}, 'o', 'Color', color{m}); hold on;
        if m==1
            title(['Monkey ',monkey{k}]);
            ylim([-75, 90]); set(gca,'YTick',[-50 0 50]); set(gca,'YTickLabel',{'-50', '0', '50'});
        else
            ylim([-80, 20]); set(gca,'YTick',[-60 -30 0]); set(gca,'YTickLabel',{'-60', '-30', '0'});
        end
        xlim([0-round(temp_winlen*0.01), temp_winlen+floor(temp_winlen*0.05)]);
        plot([min(xlim), max(xlim)], [0 0], '--k'); hold on;
        if k==1
            if m==1, txt = 'ves'; temp_y = -55; else, txt = 'vis'; temp_y = -65; end
            text(max(xlim)*0.75, temp_y, txt, 'FontWeight', 'bold', 'Color', color{m}, 'FontSize', 11);
        end
    end
end
figure(FigureIndex1)
axes('position',[0.09 0.465, 0.1 0.1]), axis off;
text(0, 0, 'PSE [deg]', 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 90);
axes('position',[0.44 0.05, 0.1 0.1]), axis off;
text(0, 0, 'Training process (Window #)', 'FontSize', 12, 'FontWeight', 'bold');
saveas(gcf,[save_fig_path,'Supplementary figure 2A.png'],'png');
close
figure(FigureIndex2)
axes('position',[0.09 0.4, 0.1 0.1]), axis off;
text(0, 0, 'Threshold [deg]', 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 90);
axes('position',[0.44 0.05, 0.1 0.1]), axis off;
text(0, 0, 'Training process (Window #)', 'FontSize', 12, 'FontWeight', 'bold');
saveas(gcf,[save_fig_path,'Supplementary figure 2B.png'],'png');
close
figure(FigureIndex3)
axes('position',[0.09 0.44, 0.1 0.1]), axis off;
text(0, 0, '\DeltaPSE [deg]', 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 90);
axes('position',[0.44 0.05, 0.1 0.1]), axis off;
text(0, 0, 'Training process (Window #)', 'FontSize', 12, 'FontWeight', 'bold');
saveas(gcf,[save_fig_path,'Supplementary figure 3B.png'],'png');
close
figure(FigureIndex4)
axes('position',[0.09 0.44, 0.1 0.1]), axis off;
text(0, 0, '\DeltaPSE [deg]', 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 90);
axes('position',[0.44 0.05, 0.1 0.1]), axis off;
text(0, 0, 'Training process (Window #)', 'FontSize', 12, 'FontWeight', 'bold');
saveas(gcf,[save_fig_path,'Supplementary figure 4B.png'],'png');
close