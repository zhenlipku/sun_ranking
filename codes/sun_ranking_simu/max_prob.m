function out = max_prob(x,y,oi,oj,data)

        % The summation of the maximized summation of pred_T(i>j)
        
        % Compute the parameters of the corresponding SUN distribution
        % using a sample where item oi is preferred to oj
       [t,n] = size(x);
       t = t+1;
       x_t = zeros(1,n);
       x_t(oi) = 1;
       x_t(oj) = -1;
       x = [x;x_t];
       y_t = 1;
       y = [y,y_t];
       D = diag(2*y-1)*x;
       D = D./sqrt(2);
       s_vec = zeros(1,t);
       for ii = 1:t
           s_vec(ii) = sqrt(D(ii,:)*D(ii,:)'+1);
       end
    
       Sinv = diag(1./s_vec);
       Delta = D'*Sinv;
       G = Sinv*(D*D'+diag(ones(1,t)))*Sinv;
       %trials_temp = 12;
       
       
       rr = com_prob_mean2(data,Delta,G); % compute all the Pred_T() via monte carol integration
       
       out_matrix = zeros(2,n*(n-1)/2);
       out_matrix(1,:) = rr.prob1;
       out_matrix(2,:) = rr.prob2;
       out = sum(max(out_matrix));
end