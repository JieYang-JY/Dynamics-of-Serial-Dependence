% Main code for the study of serial dependence dynamics - adaptive paradigm
% ----------------------------------
% This code contains several steps:
%
% 1) Data preparation
%    - reading data from each session file
%    - calculating RMS factor
%
% 2) Psychometric fits
%    - fitting psychometric function per monkey, modality, and prior type
%
% 3) Probabilistic choice model - logistic regression model
%    - modeling data per monkey and modality
%    - generalizing individual results with Fisher's method
%
% 4) Probabilistic choice model - mixed-effects logistic regression model
%    - generalizing individual results with mixed-effects model
%
% 5) Plot figures

%% STEP 1: Data preparation

clear; clc; close all;
codename = 'MAIN_AdaptiveParadigm';
codepath_name = mfilename('fullpath');
currDir = codepath_name(1 : length(codepath_name)-length(codename));
addpath(genpath(currDir))

data_path = strcat(currDir,'\Data\Adaptive paradigm\');
data_file = dir(fullfile(data_path,'*_BehData.mat'));
data_postfix = length('_BehData.mat');
FileNames = {data_file.name}';

all_data_B.currChoice = [];%monkey B
all_data_B.prevChoice = [];
all_data_B.currStim = [];
all_data_B.prevStim = [];
all_data_B.modality = [];
all_data_B.allHD = [];
all_data_D.currChoice = [];%monkey D
all_data_D.prevChoice = [];
all_data_D.currStim = [];
all_data_D.prevStim = [];
all_data_D.modality = [];
all_data_D.allHD = [];
all_data_G.currChoice = [];%monkey G
all_data_G.prevChoice = [];
all_data_G.currStim = [];
all_data_G.prevStim = [];
all_data_G.modality = [];
all_data_G.allHD = [];
all_data_J.currChoice = [];%monkey J
all_data_J.prevChoice = [];
all_data_J.currStim = [];
all_data_J.prevStim = [];
all_data_J.modality = [];
all_data_J.allHD = [];
for s = 1 : length(FileNames) % per session
    clear s_Trace;s_Trace = cell2mat(strcat(data_path, FileNames(s)));
    load(s_Trace);
    BehData.FileName

    clear temp_allHD; temp_allHD = BehData.allHD; % for RMS calculation
    clear temp_currStim; temp_currStim = BehData.allHD(BehData.IDinBatch==4);
    temp_currStim = reshape(temp_currStim,[length(temp_currStim),1]); % column shape
    clear temp_currChoice; temp_currChoice = BehData.choice(BehData.IDinBatch==4);
    temp_currChoice = reshape(temp_currChoice,[length(temp_currChoice),1]); % column shape
    temp_prevStim = []; temp_prevChoice = [];
    for b = 1 : length(temp_currChoice) % batches
        % 3 back: avg Stim --------------
        heading_p = BehData.allHD(1+(b-1)*4 : 3+(b-1)*4);
        temp_prevStim = [temp_prevStim; mean(heading_p)];
        % 3 back: mode Cho --------------
        choice_p = BehData.choice(1+(b-1)*4 : 3+(b-1)*4);
        if sum(choice_p) < 0 % previous 3 trials were left
            temp_prevChoice = [temp_prevChoice; -1];
        else % previous 3 trials were right
            temp_prevChoice = [temp_prevChoice; 1];
        end
    end
    clear temp_modality;
    if BehData.unique_StimType_test==1
        temp_modality = ones(length(temp_currChoice),1); % ves
    elseif BehData.unique_StimType_test==2
        temp_modality = ones(length(temp_currChoice),1) * 2; % vis
    else
        warning('Unexpected StimType !')
        pause
    end

    if BehData.MonkeyName=='B'
        all_data_B.currStim = [all_data_B.currStim; temp_currStim];
        all_data_B.currChoice = [all_data_B.currChoice; temp_currChoice];
        all_data_B.prevStim = [all_data_B.prevStim; temp_prevStim];
        all_data_B.prevChoice = [all_data_B.prevChoice; temp_prevChoice];
        all_data_B.modality = [all_data_B.modality; temp_modality];
        all_data_B.allHD = [all_data_B.allHD; temp_allHD];
    elseif BehData.MonkeyName=='D'
        all_data_D.currStim = [all_data_D.currStim; temp_currStim];
        all_data_D.currChoice = [all_data_D.currChoice; temp_currChoice];
        all_data_D.prevStim = [all_data_D.prevStim; temp_prevStim];
        all_data_D.prevChoice = [all_data_D.prevChoice; temp_prevChoice];
        all_data_D.modality = [all_data_D.modality; temp_modality];
        all_data_D.allHD = [all_data_D.allHD; temp_allHD];
    elseif BehData.MonkeyName=='G'
        all_data_G.currStim = [all_data_G.currStim; temp_currStim];
        all_data_G.currChoice = [all_data_G.currChoice; temp_currChoice];
        all_data_G.prevStim = [all_data_G.prevStim; temp_prevStim];
        all_data_G.prevChoice = [all_data_G.prevChoice; temp_prevChoice];
        all_data_G.modality = [all_data_G.modality; temp_modality];
        all_data_G.allHD = [all_data_G.allHD; temp_allHD];
    elseif BehData.MonkeyName=='J'
        all_data_J.currStim = [all_data_J.currStim; temp_currStim];
        all_data_J.currChoice = [all_data_J.currChoice; temp_currChoice];
        all_data_J.prevStim = [all_data_J.prevStim; temp_prevStim];
        all_data_J.prevChoice = [all_data_J.prevChoice; temp_prevChoice];
        all_data_J.modality = [all_data_J.modality; temp_modality];
        all_data_J.allHD = [all_data_J.allHD; temp_allHD];
    else
        warning('Unexpected Monkey Name !')
        pause
    end
end

%RMS factor
allHD = []; allHD = [all_data_B.allHD; all_data_D.allHD; all_data_G.allHD; all_data_J.allHD];
rms_allHD = []; rms_allHD = rms(allHD)

save_data_path = strcat(currDir,'\Data\');
save(strcat(save_data_path,'dataset_adaptive'),'all_data_*', 'rms_allHD');


%% STEP 2: Psychometric fits

clear; clc; close all;
codename = 'MAIN_AdaptiveParadigm';
codepath_name = mfilename('fullpath');
currDir = codepath_name(1 : length(codepath_name)-length(codename));
addpath(genpath(currDir))

dataset = strcat(currDir,'\Data\dataset_adaptive.mat');
load(dataset);

LEFT = -1; RIGHT = 1; prior = [LEFT RIGHT];
mod = {'ves'; 'vis'};
monkey = {'B'; 'D'; 'G'; 'J'};
PsychoFits_monks.monkey = monkey;
PsychoFits_monks.modality = mod;
PsychoFits_monks.prior = prior;
PsychoFits_monks.pfit_output = cell(length(monkey),1);
PsychoFits_monks.pseudoR2 = cell(length(monkey),1);
PsychoFits_monks.PSE = cell(length(monkey),1);
PsychoFits_monks.deltaPSE = cell(length(monkey),1);
for k = 1:length(monkey) % per monkey
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
        otherwise
            warning('Monkey Name Not Defined !')
            pause
    end

    currChoice_mod = cell(length(mod),1);
    currStim_mod = cell(length(mod),1);
    prevChoice_mod = cell(length(mod),1);
    prevStim_mod = cell(length(mod),1);
    for m = 1:length(mod) % per modality
        currChoice_mod{m} = all_data.currChoice(all_data.modality==m);
        currStim_mod{m} = all_data.currStim(all_data.modality==m);
        prevChoice_mod{m} = all_data.prevChoice(all_data.modality==m);
        prevStim_mod{m} = all_data.prevStim(all_data.modality==m);
    end
    prior_mod = cell(1,length(mod));
    for m = 1:length(mod) % per modality
        for j = 1 : length(currStim_mod{m}) % per (test) trial/ batch
            if prevStim_mod{m}(j)<0   % identify whether the prior headings are negative (Left,-1) or not (Right,1)
                prior_mod{m}(length(prior_mod{m})+1, 1) = prior(1);
            else
                prior_mod{m}(length(prior_mod{m})+1, 1) = prior(2);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pfit_output = cell(length(mod),length(prior));
    pseudoR2 = cell(length(mod),length(prior));
    PSE = cell(length(mod),length(prior));
    deltaPSE = cell(length(mod), 1);
    for m = 1:length(mod) % per modality
        if (k==1 && m==1) || k~=1 % monkey B didn't perform this paradigm in vis mod
            psycho_right = [];
            fit_data_psycho_cum = cell(1,length(prior));
            clear temp_prior_mod temp_anaHD unique_heading temp_choice;
            temp_prior_mod = prior_mod{m};
            temp_choice = currChoice_mod{m};
            temp_anaHD = currStim_mod{m};
            unique_heading = unique(currStim_mod{m});
            for p = 1:length(prior) % per prior
                for h = 1:length(unique_heading) % per unique heading
                    clear temp_heading; temp_heading = logical(temp_prior_mod==prior(p) & temp_anaHD==unique_heading(h));
                    clear temp_right_choice_trials;temp_right_choice_trials = (temp_heading & temp_choice==RIGHT);
                    psycho_right(p,h) = 1*sum(temp_right_choice_trials) / sum(temp_heading);
                    fit_data_psycho_cum{p}(h, 1) = unique_heading(h);
                    fit_data_psycho_cum{p}(h, 2) = psycho_right(p,h);
                    fit_data_psycho_cum{p}(h, 3) = sum(temp_heading);
                end
                clear temp_index;temp_index = find(fit_data_psycho_cum{p}(:,3)==0);
                if temp_index % if the heading hasn't been tested in the given prior type,
                    fit_data_psycho_cum{p}(temp_index,:) = []; % then remove it from the input matrix for PsychoFit
                end
                pfit_output{m,p} = getPsychoFit(fit_data_psycho_cum{p}(:,:), fit_data_psycho_cum{p}(:,1));
                pseudoR2{m,p} = pfit_output{m,p}.pseudoR2;
                PSE{m,p} = pfit_output{m,p}.bias;
            end
            deltaPSE{m} = PSE{m,1} - PSE{m,2};
        end
    end
    PsychoFits_monks.pfit_output{k} = pfit_output;
    PsychoFits_monks.pseudoR2{k} = pseudoR2;
    PsychoFits_monks.PSE{k} = PSE;
    PsychoFits_monks.deltaPSE{k} = deltaPSE;
end

save_results_path = strcat(currDir,'\Results\Adaptive paradigm\');
save(strcat(save_results_path,'PsychoFits_adaptive'),'PsychoFits_monks');


%% STEP 3: Probabilistic choice model - logistic regression model

clear; clc; close all;
codename = 'MAIN_AdaptiveParadigm';
codepath_name = mfilename('fullpath');
currDir = codepath_name(1 : length(codepath_name)-length(codename));
addpath(genpath(currDir))

dataset = strcat(currDir,'\Data\dataset_adaptive.mat');
load(dataset);

LEFT = -1; RIGHT = 1; prior = [LEFT RIGHT];
mod = {'ves'; 'vis'};
monkey = {'B'; 'D'; 'G'; 'J'};
LogRegModel_monks.monkey = monkey;
LogRegModel_monks.modality = mod;
LogRegModel_monks.model = cell(length(monkey), length(mod));
LogRegModel_monks.Betas_SE_t_Ps = cell(length(monkey), length(mod));
for k = 1:length(monkey) % per monkey
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
        otherwise
            warning('Monkey Name Not Defined !')
            pause
    end

    currChoice_mod = cell(length(mod),1);
    currStim_mod = cell(length(mod),1);
    prevChoice_mod = cell(length(mod),1);
    prevStim_mod = cell(length(mod),1);
    for m = 1:length(mod) % per modality
        currChoice_mod{m} = all_data.currChoice(all_data.modality==m);
        currChoice_mod{m}(currChoice_mod{m}==-1) = 0; % probability range: [0, 1]
        currStim_mod{m} = all_data.currStim(all_data.modality==m);
        currStim_norm_mod{m} = currStim_mod{m} ./ rms_allHD; % normalizing
        prevChoice_mod{m} = all_data.prevChoice(all_data.modality==m);
        prevStim_mod{m} = all_data.prevStim(all_data.modality==m);
        prevStim_norm_mod{m} = prevStim_mod{m} ./ rms_allHD; % normalizing
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear mdl_glmfit Betas_SE_t_Ps
    for m = 1:length(mod) % per modality
        if (k==1 && m==1) || k~=1 % monkey B didn't perform this paradigm in vis mod
            mdl_glmfit = fitglm([currStim_norm_mod{m} prevStim_norm_mod{m} prevChoice_mod{m}],currChoice_mod{m},'Distribution','binomial','Link','logit');
            Betas_SE_t_Ps = table2array(mdl_glmfit{m}.Coefficients);
            LogRegModel_monks.model{k,m} = mdl_glmfit;
            LogRegModel_monks.Betas_SE_t_Ps{k,m} = Betas_SE_t_Ps;
        end
    end
end

%Fisher's method
FishersMethod.p_FM_prevStim = cell(length(mod),1);
FishersMethod.chisquare_prevStim = cell(length(mod),1);
FishersMethod.p_FM_prevCho = cell(length(mod),1);
FishersMethod.chisquare_prevCho = cell(length(mod),1);
for m = 1:length(mod)
    temp_b_prevStim = [];
    temp_b_prevCho = [];
    temp_SE_prevStim = [];
    temp_SE_prevCho = [];
    temp_t_prevStim = [];
    temp_t_prevCho = [];
    temp_plist_prevStim = [];
    temp_plist_prevCho = [];
    for k = 1:length(monkey)
        if (k==1 && m==1) || k~=1 % monkey B didn't perform this paradigm in vis mod
            temp_plist_prevStim = [temp_plist_prevStim LogRegModel_monks.Betas_SE_t_Ps{k}{m}(3,4)];
            temp_plist_prevCho = [temp_plist_prevCho LogRegModel_monks.Betas_SE_t_Ps{k}{m}(4,4)];
        end
    end
    [FishersMethod.p_FM_prevStim{m}, FishersMethod.chisquare_prevStim{m}] = getFishersP(temp_plist_prevStim);
    [FishersMethod.p_FM_prevCho{m}, FishersMethod.chisquare_prevCho{m}] = getFishersP(temp_plist_prevCho);
end

save_results_path = strcat(currDir,'\Results\Adaptive paradigm\');
save(strcat(save_results_path,'LogRegModel_adaptive'),'LogRegModel_monks', 'FishersMethod');


%% STEP 4: Probabilistic choice model - mixed-effects logistic regression model

clear; clc; close all;
codename = 'MAIN_AdaptiveParadigm';
codepath_name = mfilename('fullpath');
currDir = codepath_name(1 : length(codepath_name)-length(codename));
addpath(genpath(currDir))

dataset = strcat(currDir,'\Data\dataset_adaptive.mat');
load(dataset);

% Pooling data
all_data.currChoice = [all_data_B.currChoice; all_data_D.currChoice; all_data_G.currChoice; all_data_J.currChoice];
all_data.prevChoice = [all_data_B.prevChoice; all_data_D.prevChoice; all_data_G.prevChoice; all_data_J.prevChoice];
all_data.currStim = [all_data_B.currStim; all_data_D.currStim; all_data_G.currStim; all_data_J.currStim];
all_data.prevStim = [all_data_B.prevStim; all_data_D.prevStim; all_data_G.prevStim; all_data_J.prevStim];
all_data.modality = [all_data_B.modality; all_data_D.modality; all_data_G.modality; all_data_J.modality];
all_data.monkey = [ones(length(all_data_B.modality),1); ones(length(all_data_D.modality),1)*2;
    ones(length(all_data_G.modality),1)*3; ones(length(all_data_J.modality),1)*4]; % 1: monkey B // 2: monkey D // 3: monkey G // 4: monkey J

LEFT = -1; RIGHT = 1; prior = [LEFT RIGHT];
mod = {'ves'; 'vis'};
monkey = {'B'; 'D'; 'G'; 'J'};
MixEff_LogRegModel.monkey = monkey;
MixEff_LogRegModel.modality = mod;
MixEff_LogRegModel.model = cell(length(mod),1);
MixEff_LogRegModel.Betas_SE_t_Ps = cell(length(mod),1);

currChoice_mod = cell(length(mod),1);
currStim_mod = cell(length(mod),1);
prevChoice_mod = cell(length(mod),1);
prevStim_mod = cell(length(mod),1);
for m = 1:length(mod) % per modality
    currChoice_mod{m} = all_data.currChoice(all_data.modality==m);
    currChoice_mod{m}(currChoice_mod{m}==-1) = 0; % probability range: [0, 1]
    currStim_mod{m} = all_data.currStim(all_data.modality==m);
    currStim_norm_mod{m} = currStim_mod{m} ./ rms_allHD; % normalizing
    prevChoice_mod{m} = all_data.prevChoice(all_data.modality==m);
    prevStim_mod{m} = all_data.prevStim(all_data.modality==m);
    prevStim_norm_mod{m} = prevStim_mod{m} ./ rms_allHD; % normalizing
    monkey_mod{m} = all_data.monkey(all_data.modality==m);
end
Betas_SE_t_Ps = cell(length(mod), 1);
glme = cell(length(mod), 1);
for m = 1:length(mod)
    mod{m}

    temp_currCho = currChoice_mod{m};
    temp_currStim = currStim_norm_mod{m};
    temp_prevCho = prevChoice_mod{m};
    temp_prevStim = prevStim_norm_mod{m};
    temp_monkey = monkey_mod{m};
    varNames = {'currCho','currStim','prevStim','prevCho','monkey'};
    input = table(temp_currCho,temp_currStim,temp_prevStim,temp_prevCho,temp_monkey, 'VariableNames',varNames);
    glme{m} = fitglme(input,'currCho ~ 1 + currStim + prevStim + prevCho + (currStim|monkey) + (prevStim|monkey) + (prevCho|monkey)', ...
        'Distribution','Binomial','Link','logit','FitMethod','Laplace');
    Betas_SE_t_Ps{m} = [glme{m}.Coefficients(:,2:4) glme{m}.Coefficients(:,6)];

end
MixEff_LogRegModel.model = glme;
MixEff_LogRegModel.Betas_SE_t_Ps = Betas_SE_t_Ps;

save_results_path = strcat(currDir,'\Results\Adaptive paradigm\');
save(strcat(save_results_path,'MixEff_LogRegModel_adaptive'),'MixEff_LogRegModel');


%% STEP 5: Plot figures

clear; clc; close all;
codename = 'MAIN_AdaptiveParadigm';
codepath_name = mfilename('fullpath');
currDir = codepath_name(1 : length(codepath_name)-length(codename));
addpath(genpath(currDir))


path_PsychoFits = strcat(currDir,'\Results\Adaptive paradigm\PsychoFits_adaptive.mat');
load(path_PsychoFits);
disp('ploting: Supplementary figure 5')
prior = PsychoFits_monks.prior % -1: LEFT; 1: RIGHT
mod = PsychoFits_monks.modality
monkey = PsychoFits_monks.monkey
length_filename = length('PsychoFits_adaptive.mat');
save_fig_path = path_PsychoFits(1:(end-length_filename));
color{1} = [0 0 1]; % blue
color{2} = [1 0 0]; % red
prior_type{1}='left prior'; prior_type{2}='right prior';
for k = 1 : length(monkey)
    for m = 1 : length(mod)
        if (k==1 && m==1) || k~=1 % monkey B didn't perform this paradigm in vis mod
            figure(1)
            for p = 1:length(prior) % per prior
                xi = PsychoFits_monks.pfit_output{k}{m, p}.xi;
                pfitcurve = PsychoFits_monks.pfit_output{k}{m, p}.pfitcurve;
                unique_heading = PsychoFits_monks.pfit_output{k}{m, p}.input(:,1);
                prop_Rchoice = PsychoFits_monks.pfit_output{k}{m, p}.input(:,2) ./ PsychoFits_monks.pfit_output{k}{m, p}.input(:,3);
                if p==1
                    temp_color = color{m};
                else
                    temp_color = color{m} + 0.5; temp_color(temp_color>1) = 1;
                end
                plot(unique_heading, prop_Rchoice, 'o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerEdgeColor', temp_color);
                hold on;
                plot(xi, pfitcurve, '-', 'LineWidth', 2, 'Color', temp_color);
                hold on;
                legend_txt{p*2-1} = prior_type{p};
                legend_txt{p*2} = [''];
                hold on;
            end
            xlim([min(xi),max(xi)]); ylim([0,1]);
            ylabel('Prop. Rightward Choices'); xlabel('Heading [deg]');
            title(['Psychometric plot: monkey ',monkey{k},' (',mod{m},')']);
            set(gca,'XTick',min(xi):0.5*max(xi):max(xi));
            set(gca,'XTickLabel',{num2str(min(xi)),num2str(0.5*min(xi)),'0',num2str(0.5*max(xi)),num2str(max(xi))});
            set(gca,'YTick',0:0.5:1);
            plot([0 0], [0 1], '--k'); hold on; plot([min(xi) max(xi)], [0.5 0.5], '--k'); hold on;
            legend(legend_txt, 'Location', 'southeast'); legend('boxoff');
            saveas(gcf,[save_fig_path,'PsychoPlot_monkey',monkey{k},'_',mod{m},'.png'],'png');
            close
        end
    end
end

path_LogRegModel = strcat(currDir,'\Results\Adaptive paradigm\LogRegModel_adaptive.mat');
load(path_LogRegModel);
disp('ploting: Figure 3B')
mod = LogRegModel_monks.modality
monkey = LogRegModel_monks.monkey
length_filename = length('LogRegModel_adaptive.mat');
save_fig_path = path_LogRegModel(1:(end-length_filename));
color{1} = 'b'; color{2} = 'r';
shape{1} = '^'; shape{2} = 'v'; shape{3} = 's'; shape{4} = 'diamond';
figure(1);
counter = 0;
for m = 1 : length(mod) 
    for k = 1 : length(monkey)
        if (k==1 && m==1) || k~=1 % monkey B didn't perform this paradigm in vis mod
            betaSC{m, 1}(k,1:2) = [LogRegModel_monks.Betas_SE_t_Ps{k}{m}(3, 1) LogRegModel_monks.Betas_SE_t_Ps{k}{m}(4, 1)];
            betaSC_err{m, 1}(k,1:2) = [LogRegModel_monks.Betas_SE_t_Ps{k}{m}(3, 2) LogRegModel_monks.Betas_SE_t_Ps{k}{m}(4, 2)];
            plot(betaSC{m, 1}(k,1), betaSC{m, 1}(k,2), shape{k}, 'Color', color{m}, 'LineWidth', 2, 'MarkerSize', 10); hold on;
            counter = counter + 1;
            legend_txt{counter} = [mod{m}, ', monkey ', monkey{k}];
        end
    end
end
xlim([-1.5, 1.5]);ylim([-1.5, 1.5]);
leg = legend(legend_txt, 'Location', 'southeast'); legend('boxoff');
ylabel('\beta_{prevCho} [a.u.]'); xlabel('\beta_{prevStim} [a.u.]');
title('\beta distribution: prevStim VS prevCho');
plot([0 0], [-1.5, 1.5], '--k'); hold on;
plot([-1.5, 1.5], [0 0], '--k'); hold on;
for m = 1 : length(mod) 
    for k = 1 : length(monkey)
        if (k==1 && m==1) || k~=1 % monkey B didn't perform this paradigm in vis mod
            errorbar(betaSC{m, 1}(k,1), betaSC{m, 1}(k,2), betaSC_err{m, 1}(k,1), ...
                'Color', [0.6 0.6 0.6], 'LineWidth', 1); hold on;
            errorbar(betaSC{m, 1}(k,1), betaSC{m, 1}(k,2), betaSC_err{m, 1}(k,2), ...
                'horizontal', 'Color', [0.6 0.6 0.6], 'LineWidth', 1); hold on;
        end
    end
end
saveas(gcf,[save_fig_path,'prevBeta_stimVScho.png'],'png');
close
