function sig_matrix = get_sig_matrix(results, pvalues, indices_of_tracts, num_of_ROIs, sig_threshold)

sig_matrix = zeros(num_of_ROIs);

for kk = 1:length(results)
    if pvalues(kk) < sig_threshold && results(kk) > 0
        sig_indices = indices_of_tracts(kk,:);
        sig_matrix(sig_indices(1),sig_indices(2)) = 1;
        sig_matrix(sig_indices(2),sig_indices(1)) = 1;
    end
end
