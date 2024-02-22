function out_matrix = compair_matrix(vec1,vec2,n,data)
out_matrix = zeros(n,n);
lout = vec1-vec2;
for ss = 1:length(lout)
    if lout(ss) >=0 
        out_matrix(data(ss,1),data(ss,2)) = 1;
    else
        out_matrix(data(ss,2),data(ss,1)) = 1;
    end
end

end