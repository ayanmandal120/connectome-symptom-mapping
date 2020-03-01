function analysis = connectome_symptom_mapping()

% Written by Ayan Mandal (asm253@georgetown.edu) 8/23/18

% This code can use either mass bivariate correlations or support vector
% regression to determine associations between white matter structures 
% and behavior. The code is specifically designed to use the data structures
% of connectome data given to the Turkeltaub lab from the Medaglia lab. 

% The function will output a struct  with two fields corresponing to the 
% input and output of the analysis. Input gives you information about what 
% went into the analysis (i.e. what behavior was analyzed, how many permutations were
% performed, multivariate or univariate analysis, etc.). Output contains
% the results of the analysis (i.e. which tracts were significant, what does
% the adjacency matrix of significant tracts look like, etc.)


% Step 1: Decisions about the analysis

FA_correction_on_tracts = false; % regress total FA out of each individual FA vector
control_by_healthy_tract = true; % calculate z scores for patient tracts via mean and standard deviation of control tract
SVR = true; % if true, use SVR connectome sympmtom mapping, if false, use mass univariate connectome symptom mapping
permutations = 1000;
scale = 'scale60';
sig_threshold = 0.005; % set your p value to determine signifcance
[~,~, y_column] = xlsread('DTI_Behaviors.xlsx', 'B1:B40'); % read in the dependent variable from the Excel spreadsheet
percent_nonzero_FAs = 0.9; % what percentage of controls need to have the tract for it to be counted as legitimate?



% Step 2: Load appropriate data structures

load('COGNEW_TURKELTAUB_DORISDUKE_STROKE_NETWORKS');
load('COGNEW_TURKELTAUB_DORISDUKE_CONTROL_NETWORKS');

ROI_file = ['ROI_', scale, '.mat'];

load(ROI_file); 

num_of_ROIs = numel(ROI);

matrix_3D_stroke = make_3D_DTI_matrix(Stroke_Data, num_of_ROIs, scale); 
matrix_3D_control = make_3D_DTI_matrix(Control_Data, num_of_ROIs, scale); 


ROI_matrix = make_ROI_matrix(ROI);


% Step 3: Get y (behavior) and X (matrix of predictors)

y_label = y_column{1};
y = cell2mat(y_column(2:end));

% form X from matrix_3D but store coordinates of each column of X

X = []; % X will be our matrix of features (or vectors of FA values for each tract)
indices_of_tracts = [];
Tracts = [];


tract_index = 0;
tracts_to_skip = mask_tracts(matrix_3D_control, percent_nonzero_FAs); % this function finds tracts that should be excluded given the condition that a tract should have a certain amount of nonzero values in controls to be considered valid   

 for ii = 1:num_of_ROIs
     for jj = (ii+1):num_of_ROIs
         
         tract_index = tract_index+1;
         
         if any(tract_index==tracts_to_skip)
             continue
         end
         
         individual_tract = squeeze(matrix_3D_stroke(ii,jj,:));
         
         if control_by_healthy_tract % calculate z scores for patient tracts via mean and standard deviation of control tract
             healthy_tract = squeeze(matrix_3D_control(ii,jj,:));
             individual_tract = (individual_tract - mean(healthy_tract)) ./ std(healthy_tract);
         end
         
         
         X = [X, individual_tract];
         indices_of_tracts = [indices_of_tracts; ii, jj];
         Tracts = [Tracts; ROI_matrix(ii,jj)];
     end
 end
 
% Remove rows (subjects) in X for which there are no observations in y 

X(isnan(y),:) = [];
y(isnan(y)) = [];

y = y ./ max(y); % dependent variable between 0 and 1



Total_FA = sum(X,2);

[~, p_value_matrix] = corrcoef(y, Total_FA);
FA_corr_p_value = p_value_matrix(2,1);

disp(FA_corr_p_value)

if FA_corr_p_value < 0.1 % if there's a trending relationship between total_FA and behavior, control for total FA
    y = residualize(y, Total_FA);
    FA_correction_on_behavior = true;
else
    FA_correction_on_behavior = false;
end

if FA_correction_on_tracts
    X = residualize(X, Total_FA); 
end


% Step 4: Calculate tract-wise p values
% (analagous to voxelwise p values if comparing to lesion symptom mappping)


% We will calucalte beta values and calculate p values via permutation testing if we're doing SVR connectome symptom mapping
% Or, we will calculate r values and p values if we're doing mass univariate connectome symptom mapping

if SVR == true
   
    Mdl = fitrsvm(X,y, 'KernelFunction', 'rbf');
    
    w = Mdl.Alpha' * Mdl.SupportVectors; % compute beta values
    beta_scale = 10/max(abs(w));
    results = w'*beta_scale;
    
    perms = zeros(length(results),permutations);
    
    for tt = 1:permutations % generate random beta graphs
        y_rand = y(randperm(length(y)));
        %    Mdl_rand = fitrsvm(X, y_rand, 'KernelFunction', 'rbf', 'BoxConstraint', optim_C, ...
        %        'KernelScale',optim_gamma);
        
        Mdl_rand = fitrsvm(X, y_rand, 'KernelFunction', 'rbf');
        w_rand = Mdl_rand.Alpha' * Mdl_rand.SupportVectors;
        perms(:,tt) = w_rand'*beta_scale; % distributions of beta values are stored in perms
        disp(tt)
    end
    
    pvalues = mean(perms>results,2);
    
else
    [results, pvalues] = get_r_and_p_values(X,y); % mass univariate approach
end


sig_indices = and(results>0, pvalues<sig_threshold); % one-tailed test, so we want betas/r values to be greater than zero

sig_tracts = [Tracts(sig_indices), num2cell([results(sig_indices), pvalues(sig_indices)])];

sig_matrix = get_sig_matrix(results, pvalues, indices_of_tracts, num_of_ROIs, sig_threshold);

% Step 5: Calculate component-wise p values
% (analagous to cluster-wise family error rate if comparing to lesion
% symptom mapping)

comp_p_dist = zeros(permutations,1);

for tt = 1:permutations
    
    if SVR == false 
        y_rand = y(randperm(length(y)));
        [pseudo_results, pseudo_pvalues] = get_r_and_p_values(X,y_rand);
    else
        pseudo_results = perms(:,tt); % no need to recalculate beta distributions when doing SVR since we already did that for the voxelwise p values
        perm_dist = perms;
        perm_dist(:,tt) = [];
        pseudo_pvalues = mean(perm_dist>pseudo_results,2);
    end
        
    pseudo_sig_matrix = get_sig_matrix(pseudo_results, pseudo_pvalues, indices_of_tracts, num_of_ROIs, sig_threshold);
    pseudo_comp_sizes = get_num_edges_components(pseudo_sig_matrix);
    comp_p_dist(tt) = max(pseudo_comp_sizes);
    disp(tt)
end
 
histogram(comp_p_dist)

comp_sizes = get_num_edges_components(sig_matrix);

comp_sizes(comp_sizes==0) = [];
comp_sizes = comp_sizes';

comp_p_values = [];

for cc = 1:length(comp_sizes)
    comp_p_values = [comp_p_values; mean(comp_p_dist>comp_sizes(cc))];
end


analysis.output = struct('sig_tracts',{sig_tracts},'sig_matrix', ... 
    sig_matrix, 'comp_p_values', comp_p_values, 'comp_p_dist', comp_p_dist);

analysis.input = struct('behavior', y_label, 'SVR', SVR, 'FA_corr_p_value', FA_corr_p_value, 'FA_correction_on_behavior', ...
    FA_correction_on_behavior, 'FA_correction_on_tracts', FA_correction_on_tracts, ...
    'control_by_healthy_tract', control_by_healthy_tract, ...
    'permutations', permutations, 'sig_threshold', sig_threshold, 'scale', scale, 'num_tracts', ...
    length(Tracts), 'percent_nonzero_FAs', percent_nonzero_FAs);