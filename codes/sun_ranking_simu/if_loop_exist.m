function rank = if_loop_exist(A_temp)
    [n,~] = size(A_temp);
    rank.score = zeros(1,n);
    rank.logical_val = 1;
    for ii = 1:(n-1)     % to find whether loops exist
    id_temp = find(sum(A_temp,2) == n-ii);
        if ~isempty(id_temp)
            rank.score(id_temp(1)) = n-ii;
            A_temp(id_temp(1),:) = zeros(1,n);
            A_temp(:,id_temp(1)) = zeros(1,n);            
        else
            rank.logical_val = 0;
        end
    end   
end