function out = rank_ver(com_matrix)
    [n,~] = size(com_matrix);
    A = com_matrix >= 0.5;
    out.score = sum(A,2);
    if (sum( sort(out.score) == (0:(n-1))') == n)
        out.logical_val = 1;
    else
        out_temp = if_loop_exist(A);
        out.score = out_temp.score;
        out.logical_val = out_temp.logical_val;
    end
end