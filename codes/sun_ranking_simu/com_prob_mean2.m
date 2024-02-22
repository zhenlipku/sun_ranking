function rr = com_prob_mean2(data,Delta,G)

len = length(data(:,1));
% Generate random samples from the corresponding SUN distribution with the 
% provided Gamma and Delta inputs.

randomn = sundraw_new2(Delta,G); 
[a,~] = size(randomn);
index1 = zeros(a,len);
index2 = zeros(a,len);
for ll = 1:len
    temp1 = data(ll,1);
    temp2 = data(ll,2);
    index1(:,ll) = (1:a)'+(temp1-1)*a;
    index2(:,ll) = (1:a)'+(temp2-1)*a;
end
out1 = sum((randomn(index1) - randomn(index2))>0)/a;
out2 = mean(randomn);
rr.prob1 = out1;
rr.prob2 = 1-rr.prob1;
rr.mean = out2;
 end