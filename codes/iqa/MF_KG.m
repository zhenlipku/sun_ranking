function out = MF_KG(x,y,oi,oj,n,all_pairs_D)
x_t = zeros(1,n);
x_t(oi) = 1;
x_t(oj) = -1;
x_temp = [x;x_t];
y_temp = [y,1];
out_temp = est_para_mf(x_temp,y_temp);
Mean_temp = out_temp.Mean;
Sigma_temp = out_temp.Sigma;
out = est_prob_compute(Mean_temp,Sigma_temp,all_pairs_D);

end