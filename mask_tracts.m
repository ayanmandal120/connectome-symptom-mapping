function tracts_to_skip = mask_tracts(matrix_3D_controls, percent_nonzero_FAs)

% returns a vector of tract indices to skip

% percent_nonzero_FAs = percentage of controls that need to have a nonzero FA value in a tract for the tract to be included

num_of_ROIs = length(matrix_3D_controls(:,1));

N = length(matrix_3D_controls(1,1,:)); % number of controls


cutoff = N * percent_nonzero_FAs;

tract_index = 1;
tracts_to_skip = [];


% FA_controls = [];

 for ii = 1:num_of_ROIs
     for jj = (ii+1):num_of_ROIs
         individual_tract = squeeze(matrix_3D_controls(ii,jj,:));
         individual_tract_bin = logical(individual_tract);
         
         
         if sum(individual_tract_bin)<cutoff
             tracts_to_skip = [tracts_to_skip, tract_index];
         end
         tract_index = tract_index+1;
     end
 end
disp(tract_index)
disp(length(tracts_to_skip)/tract_index)
 