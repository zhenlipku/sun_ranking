function rrr = tran_prob_compute(oi,oj,x,y) 
% Determine the transition probability as outlined in the computational details provided in the appendix
    [~,n] = size(x);
    x_t = zeros(1,n);
    x_t(oi) = 1;
    x_t(oj) = -1;
    x_new = [x;x_t];
    y_new = [y,1];
    out = generateG(x_new,y_new);
    [t,~] = size(out.G);
    
    if t<=25
        rrr = mvncdf(zeros(1,t),zeros(1,t),out.G);
    else
        rrr= mvncdf_fast(-Inf*ones(1,t),zeros(1,t),out.G,1000).prob;
    end
end