function out = if_forms_a_rank(A)
[n,~] = size(A);
out.score = sum(A,2);
out.logic = (sum( sort(out.score) == (0:(n-1))') == n);
end