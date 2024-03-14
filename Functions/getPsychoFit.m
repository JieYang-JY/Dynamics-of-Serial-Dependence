function [pfit_output]= getPsychoFit(pfit_input,dir)

%% psychometric fits:
pfit_input(pfit_input(1:end,2)==0,2)=1e-6; %AZ pfit program doesn't like pefect 0s
pfit_input(pfit_input(1:end,2)==1,2)=1-(1e-6); %AZ pfit program doesn't like perfect 1s
pfit_input(:,2) = round(pfit_input(:,2).*pfit_input(:,3));
options.sigmoidName = 'norm';
options.expType = 'equalAsymptote';

result = psignifit(pfit_input,options);

result_std = getStandardParameters(result);
bias = result_std(1);
thr = result_std(2);

%% goodness-of-fit (pseudo-r-square):
% calculate the deviance:
[~, deviance, ~, ~] = getDeviance(result,'');
% calculate the null deviance (based on the getDeviance function from the toolbox psignifit):
pPred = result.psiHandle(result.data(:,1));
% change predicted to the null predicted: The mean of y values.
pPred(1:end) = mean(result.data(:,2) ./ result.data(:,3));
pMeasured = result.data(:,2) ./ result.data(:,3);
loglikelihoodPred = result.data(:,2) .* log(pPred) + (result.data(:,3) - result.data(:,2)) .* log((1 - pPred));
loglikelihoodMeasured = result.data(:,2) .* log(pMeasured) + (result.data(:,3) - result.data(:,2)) .* log((1 - pMeasured));
loglikelihoodMeasured(pMeasured == 1) = 0;
loglikelihoodMeasured(pMeasured == 0) = 0;
devianceResiduals = -2*sign(pMeasured - pPred).*(loglikelihoodMeasured - loglikelihoodPred);
deviance_null = sum(abs(devianceResiduals));
pseudoR2 =(deviance_null - deviance)/deviance_null;


%% results
pfit_output.input = pfit_input;
pfit_output.bias = bias;
pfit_output.thresh = thr;
pfit_output.pseudoR2 = pseudoR2;

% Curve
xi = linspace(min(dir),max(dir),1000);
xi = reshape(xi,[length(xi),1]);
pfit_output.xi = xi;
y = result.psiHandle(xi);
y = reshape(y,[length(y),1]);
pfit_output.pfitcurve = y;


