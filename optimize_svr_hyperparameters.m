function results = optimize_svr_hyperparameters(X,y)

%  UNFINISHED

% This code optimizes the hyperparameters for an SVR model by selecting
% parameters which maximize prediction accuracy and reproducibility. The
% procedure is the same as that used in Zhang et al., 2014


% removing NaN's from the data

X(isnan(y),:) = [];
y(isnan(y)) = [];



C_combs = [1, 10, 20, 30, 40, 50, 60, 70, 80];
gamma_combs = [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30];

%C_combs = 0.1:0.05:1;
%gamma_combs = 0.5:0.1:1.5;


num_combinations = length(C_combs)*length(gamma_combs);

pred_accuracy = zeros(num_combinations,1);
reproducibility = zeros(num_combinations,1);
C_array = zeros(num_combinations,1);
gamma_array = zeros(num_combinations,1);

comb_num = 1;

disp(num_combinations)

for C = C_combs
    for gamma = gamma_combs
        [pred_accuracy(comb_num), reproducibility(comb_num)] ...
             = get_param_criteria(X,y,C,gamma);
        C_array(comb_num) = C; 
        gamma_array(comb_num) = gamma;
        
        comb_num = comb_num + 1
    end
end

combine_parameters = (pred_accuracy + reproducibility)./2;
[M,I] = max(combine_parameters);

optim_C = C_array(I);
optim_gamma = gamma_array(I);
optim_pred_accuracy = pred_accuracy(I);
optim_reproducibility = reproducibility(I);

results.distributions = struct('C_combs', C_combs, 'gamma_combs', gamma_combs, ...
    'pred_accuracy', pred_accuracy, 'reproducibility', reproducibility);
results.optim = struct('optim_C', optim_C, 'optim_gamma', optim_gamma, ...
    'optim_pred_accuracy', optim_pred_accuracy, 'optim_reproducibility', optim_reproducibility);


end





function [prediction_accuracy_index, reproducibility_index] = get_param_criteria(X,y,C,gamma)

% Step 1: generate cross validation indices -- we will use these to seperate X and
% y into training and test sets


N = length(y); % number of observations
k = 5; % Zhang et al., 2014 does five fold cross validation
m = length(X(1,:)); % number of features


Indices = crossvalind('Kfold', N, k);

prediction_accuracy = zeros(k,1);
betas = zeros(m, k); % we will store beta maps/graphs to generate reproducibility index

for i = 1:k
    X_train = X(Indices~=i,:);
    X_test = X(Indices==i,:);
    y_train = y(Indices~=i);
    y_test = y(Indices==i);

    Mdl = fitrsvm(X_train,y_train,'KernelFunction','rbf', 'BoxConstraint', C, 'KernelScale', gamma);
    y_predictions = predict(Mdl,X_test);
    corr_matrix = corrcoef(y_predictions, y_test);  % corrcoef outputs a 2x2 matrix but we just want the correlation coefficient
    prediction_accuracy(i) = corr_matrix(2,1);
    
    
    w = Mdl.Alpha' * Mdl.SupportVectors; % beta values
    beta_scale = 10/max(abs(w));
    betas(:,i) = w'*beta_scale; 
end

% Mdl = fitrsvm(X,y, 'KernelFunction','rbf', 'BoxConstraint', C, 'KernelScale', gamma);
% XVMdl = crossval(Mdl,'Kfold',k); % this is a 5-fold
% predicted = kfoldPredict(XVMdl);

% predicted = 1 - predicted;

%corr_matrix = corrcoef(y, predicted);
%prediction_accuracy_index = corr_matrix(2,1);
    
prediction_accuracy_index = mean(prediction_accuracy);
    
beta_corr_matrix = tril(corrcoef(betas));
beta_corr_matrix = beta_corr_matrix - eye(k); % get rid of diagnol of 1's
reproducibility_index = sum(sum(beta_corr_matrix)) / nchoosek(k,2); % sum of the correlations of the betas, normalized by number of correlations run, gives reproducibility index

end
