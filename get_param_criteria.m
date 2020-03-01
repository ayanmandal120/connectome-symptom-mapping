function [prediction_accuracy_index, reproducibility_index] = get_param_criteria(X,y)

% UNFINISHED

% Step 1: Data preprocessing
% we will remove features from X if values are missing for a certain number of
% subjects

X(isnan(y),:) = [];
y(isnan(y)) = [];


N = length(y); % number of observations


cutoff = floor(N * 0.1);

non_zero_observations = X > 0;
num_non_zero_tracts = sum(non_zero_observations);
X(:,num_non_zero_tracts<=cutoff) = [];


% Step 2: generate cross validation indices -- we will use these to seperate X and
% y into training and test sets



k = 5; % Zhang et al., 2014 does five fold cross validation



Indices = crossvalind('Kfold', N, k);

m = length(X(1,:)); % number of features
prediction_accuracy = zeros(k,1);
betas = zeros(m, k); % we will store beta maps/graphs to generate reproducibility index


for i = 1:k
    X_train = X(Indices~=i,:);
    X_test = X(Indices==i,:);
    y_train = y(Indices~=i);
    y_test = y(Indices==i);
    
    Mdl = fitrsvm(X_train,y_train,'KernelFunction','rbf');
    y_predictions = predict(Mdl,X_test);
    corr_matrix = corrcoef(y_predictions, y_test);  % corrcoef outputs a 2x2 matrix but we just want the correlation coefficient
    prediction_accuracy(i) = corr_matrix(2,1);
    
    w = Mdl.Alpha' * Mdl.SupportVectors; % beta values
    beta_scale = 10/max(abs(w));
    betas(:,i) = w'.*beta_scale; 
end

%Mdl = fitrsvm(X,y, 'KernelFunction','rbf');
%XVMdl = crossval(Mdl,'Kfold',k); % this is a 5-fold
%predicted = kfoldPredict(XVMdl);

disp(Mdl.ModelParameters.BoxConstraint)
disp(Mdl.ModelParameters.KernelScale)

%corr_matrix = corrcoef(y, predicted);
%prediction_accuracy_index = corr_matrix(2,1);

prediction_accuracy_index = mean(prediction_accuracy);

beta_corr_matrix = tril(corrcoef(betas));
beta_corr_matrix = beta_corr_matrix - eye(k); % get rid of diagnol of 1's
reproducibility_index = sum(sum(beta_corr_matrix)) / nchoosek(k,2); % sum of the correlations of the betas, normalized by number of correlations run, gives reproducibility index

end
