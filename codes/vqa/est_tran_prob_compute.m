function rrr = est_tran_prob_compute(Mean,Sigma,oi,oj,n)
D = zeros(1,n);
D(oi) = 1/sqrt(2);
D(oj) = -1/sqrt(2);
rrr = normcdf(D*Mean*(1+D*Sigma*D')^(-0.5));
end