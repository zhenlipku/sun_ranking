function out = generateG(x,y)
    D = diag(2*y-1)*x;
    D = D./sqrt(2);
    [t,~] = size(x);
    s_vec = zeros(1,t);

        for ii = 1:t
            s_vec(ii) = sqrt(D(ii,:)*D(ii,:)'+1);
        end
        
S = diag(s_vec);
Sinv = inv(S);
out.Delta = D'*Sinv;
out.G = Sinv*(D*D'+diag(ones(1,t)))*Sinv;
end