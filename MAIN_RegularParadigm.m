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
PsychoFits_monks.PerceptLearn.pfit_output = cell(length(monkey), length(mod)); % Psychometric fits: perceptual learning
PsychoFits_monks.PerceptLearn.pseudoR2 = cell(length(monkey), length(mod));
PsychoFits_monks.PerceptLearn.PSE = cell(length(monkey), length(mod));
PsychoFits_monks.PerceptLearn.thre = cell(length(monkey), length(mod));
PsychoFits_monks.SerialDepend.pfit_output = cell(length(monkey), length(mod)); % Psychometric fits: aggregate serial dependence
PsychoFits_monks.SerialDepend.pseudoR2 = cell(length(monkey), length(mod));
PsychoFits_monks.SerialDepend.PSE = cell(length(monkey), length(mod));
PsychoFits_monks.SerialDepend.deltaPSE = cell(length(monkey), length(mod));
LogRegModel_monks.model = cell(length(monkey), length(mod)); % Modeling
LogRegModel_monks.Beta_currStim = cell(length(monkey), length(mod));
LogRegModel_monks.Beta_prevStim = cell(length(monkey), length(mod));
LogRegModel_monks.Beta_prevCho = cell(length(monkey), length(mod));
LogRegModel_monks.deltaBeta = cell(length(monkey), length(mod));
LogRegModel_monks.sm_Beta_currStim = cell(length(monkey), length(mod));
LogRegModel_monks.sm_Beta_prevStim = cell(length(monkey), length(mod));
LogRegModel_monks.sm_Beta_prevCho = cell(length(monkey), length(mod));
LogRegModel_monks.sm_deltaBeta = cell(length(monkey), length(mod));

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
        disp(['Sliding in ' cell2mat(mod(m)) ' condition...'])
        pos_windowEnd = window_size; prev_pos_windowEnd = 0; winInd = 0; % RESET
        pfit_output_PL = cell(0); threshold_PL = []; PSE_PL = []; pseudoR2_PL = [];
        pfit_output_SD = cell(0); PSE_SD = []; pseudoR2_SD = []; deltaPSE = [];
        mdl_glmfit = cell(0); Beta_currStim = []; Beta_prevStim = []; Beta_prevCho = []; deltaBeta = [];
        while prev_pos_windowEnd < length(all_data.currStim{m}) % until the last trial of given stim type
            winInd = winInd+1;

            % Psychometric fits
            clear temp_prior_type temp_choice temp_anaHD unique_heading;
            temp_prior_type = all_data.prior{m}(pos_windowEnd-window_size+1 : pos_windowEnd);
            temp_choice = all_data.currChoice{m}(pos_windowEnd-window_size+1 : pos_windowEnd);
            temp_anaHD = all_data.currStim{m}(pos_windowEnd-window_size+1 : pos_windowEnd);
            unique_heading = unique(temp_anaHD);
            %perceptual learning ----------------
            psycho_right_PL = []; fit_data_psycho_cum_PL = [];
            for h = 1:length(unique_heading) % per unique heading
                clear temp_heading; temp_heading = logical(temp_anaHD==unique_heading(h));
                clear temp_right_choice_trials;temp_right_choice_trials = (temp_heading & temp_choice==RIGHT);
                psycho_right_PL(h) = 1*sum(temp_right_choice_trials) / sum(temp_heading);
                fit_data_psycho_cum_PL(h, 1) = unique_heading(h);
                fit_data_psycho_cum_PL(h, 2) = psycho_right_PL(h);
                fit_data_psycho_cum_PL(h, 3) = sum(temp_heading);
            end
            pfit_output_PL{winInd, 1} = getPsychoFit(fit_data_psycho_cum_PL(:,:), fit_data_psycho_cum_PL(:,1));
            threshold_PL(winInd, 1) = pfit_output_PL{winInd, 1}.thresh;
            PSE_PL(winInd, 1) = pfit_output_PL{winInd, 1}.bias;
            pseudoR2_PL(winInd, 1) = pfit_output_PL{winInd, 1}.pseudoR2;
            %aggregate serial dependence ----------------
            psycho_right_SD = []; fit_data_psycho_cum_SD = cell(1,length(prior));
            for p = 1:length(prior) % per prior
                for h = 1:length(unique_heading) % per unique heading
                    clear temp_heading; temp_heading = logical(temp_prior_type==prior(p) & temp_anaHD==unique_heading(h));
                    clear temp_right_choice_trials;temp_right_choice_trials = (temp_heading & temp_choice==RIGHT);
                    psycho_right_SD(p,h) = 1*sum(temp_right_choice_trials) / sum(temp_heading);
                    fit_data_psycho_cum_SD{p}(h, 1) = unique_heading(h);
                    fit_data_psycho_cum_SD{p}(h, 2) = psycho_right_SD(p,h);
                    fit_data_psycho_cum_SD{p}(h, 3) = sum(temp_heading);
                end
                clear temp_index;temp_index = find(fit_data_psycho_cum_SD{p}(:,3)==0);
                if temp_index % if the heading hasn't been tested in the given prior type,
                    fit_data_psycho_cum_SD{p}(temp_index,:) = []; % then remove it from the input matrix for PsychoFit
                end
                pfit_output_SD{winInd, p} = getPsychoFit(fit_data_psycho_cum_SD{p}(:,:), fit_data_psycho_cum_SD{p}(:,1));
                PSE_SD(winInd, p) = pfit_output_SD{winInd, p}.bias;
                pseudoR2_SD(winInd, p) = pfit_output_SD{winInd, p}.pseudoR2;
            end
            deltaPSE(winInd, 1) = PSE_SD(winInd, 1) - PSE_SD(winInd, 2);

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

        PsychoFits_monks.PerceptLearn.pfit_output{k, m} = pfit_output_PL; % Psychometric fits: perceptual learning
        PsychoFits_monks.PerceptLearn.pseudoR2{k, m} = pseudoR2_PL;
        PsychoFits_monks.PerceptLearn.PSE{k, m} = PSE_PL;
        PsychoFits_monks.PerceptLearn.thre{k, m} = threshold_PL;
        PsychoFits_monks.SerialDepend.pfit_output{k, m} = pfit_output_SD; % Psychometric fits: aggregate serial dependence
        PsychoFits_monks.SerialDepend.pseudoR2{k, m} = pseudoR2_SD;
        PsychoFits_monks.SerialDepend.PSE{k, m} = PSE_SD;
        PsychoFits_monks.SerialDepend.deltaPSE{k, m} = deltaPSE;
        LogRegModel_monks.model{k, m} = mdl_glmfit; % Modeling
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
    end
end

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
disp('ploting: Figure 2B, 2C & Supplementary figure 4')
for k = 1 : length(monkey)
    for m = 1 : length(mod)
        % beta_currStim
        figure(1)
        plot(LogRegModel_monks.sm_Beta_currStim{k, m}, 's', 'MarkerFaceColor', color{m}, 'MarkerEdgeColor', color{m})
        hold on; ylim([0, 14.5]);
        ylabel('\beta [a.u.] (curr)'); xlabel('Training process (Window #)');
        title(['Dynamics of model fits over training: currStim, monkey ', monkey{k}, ' in ', mod{m}]);
        set(gca,'YTick',[0 7 14]); set(gca,'YTickLabel',{'0', '7', '14'});
        saveas(gcf,[save_fig_path,'ModelFit_currStim_monkey',monkey{k},'_',mod{m},'.png'],'png');
        close
        % beta_prev
        figure(1)
        plot(LogRegModel_monks.sm_Beta_prevStim{k, m}, 's', 'MarkerFaceColor', color{m}, 'MarkerEdgeColor', color{m})
        hold on; legend_txt{1} = 'prevStim'; hold on;
        plot(LogRegModel_monks.sm_Beta_prevCho{k, m}, 's', 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', [0.6 0.6 0.6])
        hold on; legend_txt{2} = 'prevCho'; hold on;
        ylim([-1.3, 0.5]); ylabel('\beta [a.u.] (prev)'); xlabel('Training process (Window #)');
        title(['Dynamics of model fits over training: prevStim & prevCho, monkey ', monkey{k}, ' in ', mod{m}]);
        set(gca,'YTick',[-1.2 -0.6 0]); set(gca,'YTickLabel',{'-1.2', '-0.6', '0'});
        leg = legend(legend_txt, 'Location', 'southwest'); legend('boxoff');
        plot([min(xlim), max(xlim)], [0 0], '--k'); hold on;
        saveas(gcf,[save_fig_path,'ModelFit_prevStim_monkey',monkey{k},'_',mod{m},'.png'],'png');
        close
        
        % deltabeta
        figure(1)
        plot(LogRegModel_monks.sm_deltaBeta{k, m}, 's', 'MarkerFaceColor', color{m}, 'MarkerEdgeColor', color{m})
        hold on; ylim([-1.45, 0.7]);
        ylabel('\Delta\beta (\beta_{prevStim} - \beta_{prevCho}) [a.u.]'); xlabel('Training process (Window #)');
        title(['Dynamics of model fits over training: \Delta\beta, monkey ', monkey{k}, ' in ', mod{m}]);
        set(gca,'YTick',[-1.4 -0.7 0 0.7]); set(gca,'YTickLabel',{'-1.4', '-0.7', '0', '0.7'});
        plot([min(xlim), max(xlim)], [0 0], '--k'); hold on;
        saveas(gcf,[save_fig_path,'ModelFit_prevCho_monkey',monkey{k},'_',mod{m},'.png'],'png');
        close
    end
end

disp('ploting: Figure 1D') % Perceptual-learning example of psycho fits
figure(1)
x = []; y = [];  % ves
x = PsychoFits_monks.PerceptLearn.pfit_output{3, 1}{92}.input(1:end, 1);
y = PsychoFits_monks.PerceptLearn.pfit_output{3, 1}{92}.input(1:end, 2) ./ PsychoFits_monks.PerceptLearn.pfit_output{3, 1}{92}.input(1:end, 3);
plot(x, y, 'o', 'Color', color{1}, 'LineWidth', 2, 'MarkerSize', 8); hold on;
xi = []; pfitcurve_ves = [];
xi = PsychoFits_monks.PerceptLearn.pfit_output{3, 1}{92}.xi;
pfitcurve_ves = PsychoFits_monks.PerceptLearn.pfit_output{3, 1}{92}.pfitcurve;
plot(xi, pfitcurve_ves, '-', 'Color', color{1}, 'LineWidth', 2, 'MarkerSize', 8); hold on;
legend_txt{1} = 'ves';
legend_txt{2} = [''];
x = []; y = [];  % vis
x = PsychoFits_monks.PerceptLearn.pfit_output{3, 2}{60}.input(1:end, 1);
y = PsychoFits_monks.PerceptLearn.pfit_output{3, 2}{60}.input(1:end, 2) ./ PsychoFits_monks.PerceptLearn.pfit_output{3, 2}{60}.input(1:end, 3);
plot(x, y, 'o', 'Color', color{2}, 'LineWidth', 2, 'MarkerSize', 8); hold on;
xi = []; pfitcurve_vis = [];
xi = PsychoFits_monks.PerceptLearn.pfit_output{3, 2}{60}.xi;
pfitcurve_vis = PsychoFits_monks.PerceptLearn.pfit_output{3, 2}{60}.pfitcurve;
plot(xi, pfitcurve_vis, '-', 'Color', color{2}, 'LineWidth', 2, 'MarkerSize', 8); hold on;
legend_txt{3} = 'vis';
legend_txt{4} = ['']; 
legend(legend_txt, 'Location', 'northwest'); 
diff = []; index = []; % illustrating PSE
diff = abs(pfitcurve_ves - 0.5); index = find(min(diff)==diff); 
if isempty(index); index = find(-min(diff)==diff); end
plot([xi(index), xi(index)], [0 0.5], '--', 'Color', color{1}); hold on;
diff = []; index = [];
diff = abs(pfitcurve_vis - 0.5); index = find(min(diff)==diff); 
if isempty(index); index = find(-min(diff)==diff); end
plot([xi(index), xi(index)], [0 0.5], '--', 'Color', color{2}); hold on;
xlim([min(x),max(x)]); ylim([0, 1]); % configuration
plot([min(xlim), max(xlim)], [0.5 0.5], ':k'); hold on;
text(2, 0.05, 'PSE', 'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Prop. Rightward Choices'); xlabel('Heading [deg]');
title('Example psychometric plot: monkey G (Perceptual Learning)');
set(gca,'XTick',min(xi):0.5*max(xi):max(xi));
set(gca,'XTickLabel',{num2str(min(xi)),num2str(0.5*min(xi)),'0',num2str(0.5*max(xi)),num2str(max(xi))});
set(gca,'YTick',0:0.5:1);
saveas(gcf,[save_fig_path,'ExamplePsychoPlot_monkeyG_PL.png'],'png');
close

disp('Supplementary figure 3A') % Serial-dependence example of psycho fits
f = figure(1); f.Position(3) = [1100];
winInd = [92, 60];
prior_type = {'left prior'; 'right prior'};
for m = 1 : length(mod)
    axes('Position', [0.07+0.5*(m-1) 0.11 0.4 0.8])
    for p = 1 : length(prior)
        x = []; y = []; % points
        x = PsychoFits_monks.SerialDepend.pfit_output{3, m}{winInd(m), p}.input(1:end, 1);
        y = PsychoFits_monks.SerialDepend.pfit_output{3, m}{winInd(m), p}.input(1:end, 2) ./ PsychoFits_monks.SerialDepend.pfit_output{3, m}{winInd(m), p}.input(1:end, 3);
        if p==1
            temp_color = color{m};
        else
            temp_color = color{m} + 0.5; temp_color(temp_color>1) = 1;
        end
        plot(x, y, 'o', 'Color', temp_color, 'LineWidth', 2, 'MarkerSize', 8); hold on;
        xi = []; pfitcurve = []; % curve
        xi = PsychoFits_monks.SerialDepend.pfit_output{3, m}{winInd(m), p}.xi;
        pfitcurve = PsychoFits_monks.SerialDepend.pfit_output{3, m}{winInd(m), p}.pfitcurve;
        plot(xi, pfitcurve, '-', 'Color', temp_color, 'LineWidth', 2, 'MarkerSize', 8); hold on;
        legend_txt{p*2-1} = prior_type{p};
        legend_txt{p*2} = [''];
        legend(legend_txt, 'Location', 'southeast'); 
        xlim([min(xi),max(xi)]); ylim([0,1]);
        ylabel('Prop. Rightward Choices'); xlabel('Heading [deg]');
        title([mod{m} ' example']);
        set(gca,'XTick',min(xi):0.5*max(xi):max(xi));
        set(gca,'XTickLabel',{num2str(min(xi)),num2str(0.5*min(xi)),'0',num2str(0.5*max(xi)),num2str(max(xi))});
        set(gca,'YTick',0:0.5:1);
    end
end
saveas(gcf,[save_fig_path,'ExamplePsychoPlot_monkeyG_SD.png'],'png');
close

disp('ploting: Supplementary figure 2 & Supplementary figure 3B')
for k = 1 : length(monkey)
    for m = 1 : length(mod)
        % PSE
        figure(1)
        plot(PsychoFits_monks.PerceptLearn.PSE{k, m}, 'o', 'Color', color{m})
        hold on; ylim([-50, 50]);
        plot([min(xlim), max(xlim)], [0 0], '--k'); hold on;
        ylabel('PSE [deg]'); xlabel('Training process (Window #)');
        title(['Dynamics of PSE over training: monkey ', monkey{k}, ' in ', mod{m}]);
        set(gca,'YTick',[-40 0 40]); set(gca,'YTickLabel',{'-40', '0', '40'});
        saveas(gcf,[save_fig_path,'PsychoPlot_PSE_monkey',monkey{k},'_',mod{m},'.png'],'png');
        close
        % Threshold
        figure(1)
        log2_thre = log2(PsychoFits_monks.PerceptLearn.thre{k, m});
        plot(log2_thre, 'o', 'Color', color{m})
        hold on; ylim([0, 7]);
        plot([min(xlim), max(xlim)], [log2(6) log2(6)], '--k'); hold on;
        ylabel('Threshold [deg]'); xlabel('Training process (Window #)');
        title(['Dynamics of Threshold over training: monkey ', monkey{k}, ' in ', mod{m}]);
        set(gca,'YTick',[1 3 5 7]); set(gca,'YTickLabel',{'2', '8', '32', '128'});
        saveas(gcf,[save_fig_path,'PsychoPlot_Thre_monkey',monkey{k},'_',mod{m},'.png'],'png');
        close
        
        % deltaPSE
        figure(1)
        plot(PsychoFits_monks.SerialDepend.deltaPSE{k, m}, 'o', 'Color', color{m}); hold on;
        if m==1
            ylim([-75, 90]); set(gca,'YTick',[-50 0 50]); set(gca,'YTickLabel',{'-50', '0', '50'});
        else
            ylim([-80, 20]); set(gca,'YTick',[-60 -30 0]); set(gca,'YTickLabel',{'-60', '-30', '0'});
        end
        plot([min(xlim), max(xlim)], [0 0], '--k'); hold on;
        ylabel('\DeltaPSE [deg]'); xlabel('Training process (Window #)');
        title(['Dynamics of \DeltaPSE over training: monkey ', monkey{k}, ' in ', mod{m}]);
        saveas(gcf,[save_fig_path,'PsychoPlot_dPSE_monkey',monkey{k},'_',mod{m},'.png'],'png');
        close
    end
end

