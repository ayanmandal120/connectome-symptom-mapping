function [rvalues, pvalues] = get_r_and_p_values(X,y)


% This function returns r values and p values for correlations between
% every variable in X and y

num_features = length(X(1,:));

rvalues = zeros(num_features,1);
pvalues = zeros(num_features,1);

for ii = 1:num_features

    FA_tract = X(:,ii);
    
    [r,p] = corrcoef(y, FA_tract,'Rows','complete');
    
    rvalues(ii) = r(1,2);
    pvalues(ii) = p(1,2);
    
end
    