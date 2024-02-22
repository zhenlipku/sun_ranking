function out = com_prob(Delta,G,index1,index2,a)
randomn_temp = sundraw(Delta,G);
out = sum((randomn_temp(index1) - randomn_temp(index2))>0)/a;
end