function residuals = residualize(variables, covariates)

% This function fits a linear model between a covariate and a variable, and
% returns the residuals from that linear model. This is iterated across
% every variable in "variables".

for vars = 1:length(variables(1,:))
    mdl = fitlm(covariates, variables(:,vars));
    variables(:,vars) = mdl.Residuals.raw;
end

residuals = variables;