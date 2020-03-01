function ROI_matrix = make_ROI_matrix(ROI)

% This function converts a string vector of ROI's into a matrix where each
% element (of index i and j)represents the connection between ROI i and ROI j

ROI_matrix = cell(length(ROI));

for ii = 1:length(ROI)
    for jj = 1:length(ROI)
        ROI_matrix{ii,jj} = [ROI{ii}, '_', ROI{jj}];
    end
end
