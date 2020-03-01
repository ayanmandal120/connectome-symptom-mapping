function matrix_3D = make_3D_DTI_matrix(data_structure, num_of_ROIs, scale)

% This function returns a 3 dimensional matrix. Adjacency matrices comprise
% the first two dimensions -- the indices in these matrices indicate the FA
% value between tracts connecting ROI i and ROI j for a specific subject.
% Matrices for each participant are "stacked" upon one another, so the
% length of the third dimension corresponds to the number of partipants in
% the sample.

S = data_structure;

ID = fieldnames(S);

matrix_3D = zeros(num_of_ROIs, num_of_ROIs, length(ID));

for nn = 1:length(ID)
    
    if strcmpi('scale33',scale)   
        matrix = S.(ID{nn}).TIME1.DTI.LAUSANNE_scale33.FA.connectivity;
    elseif strcmpi('scale60',scale)
        matrix = S.(ID{nn}).TIME1.DTI.LAUSANNE_scale60.FA.connectivity;
    elseif strcmpi('scale125',scale)
        matrix = S.(ID{nn}).TIME1.DTI.LAUSANNE_scale125.FA.connectivity;
    elseif strcmpi('scale250',scale)
        matrix = S.(ID{nn}).TIME1.DTI.LAUSANNE_scale250.FA.connectivity;
    elseif strcmpi('scale500',scale)
        matrix = S.(ID{nn}).TIME1.DTI.LAUSANNE_scale500.FA.connectivity;
    end
       
    while length(matrix(:,1)) ~= num_of_ROIs % some patients are missing an ROI
        matrix = [matrix, zeros(length(matrix(:,1)),1)]; % zero pad until right size
        matrix = [matrix; zeros(1,length(matrix(1,:)))];
    end
    matrix_3D(:,:,nn) = matrix;
end
