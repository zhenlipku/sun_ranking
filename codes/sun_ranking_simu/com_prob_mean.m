function rr = com_prob_mean(data,Delta,G,trials)
len = length(data(:,1));
randomn = sundraw(Delta,G);
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
parfor ii =1:trials
%ramdomn_temp = sundraw(Delta,G);
out1 = out1+com_prob(Delta,G,index1,index2,a);%sum((ramdomn_temp(index1) - ramdomn_temp(index2))>0)/a;
out2 = out2+com_mean(Delta,G);
%out2 = out2+mean(randomn_temp);
end
rr.prob1 = out1/(trials+1);
rr.prob2 = 1-rr.prob1;
rr.mean = out2/(trials+1);
 end