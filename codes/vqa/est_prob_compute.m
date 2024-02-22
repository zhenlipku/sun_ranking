function rrr = est_prob_compute(Mean,Sigma,all_pairs_D)
D_mean = all_pairs_D*Mean;
D_sigma = all_pairs_D*Sigma*all_pairs_D';
rrr = normcdf(D_mean./sqrt(diag(D_sigma)));
rrr = sum(max(transpose([1-rrr rrr])));
end