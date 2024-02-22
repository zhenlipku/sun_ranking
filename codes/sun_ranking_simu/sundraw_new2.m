function out = sundraw_new2(Delta,G)
n_size = 5000;  % number of samples
[n,~] = size(Delta);
[t,~] = size(G);
if t < 2
    dd = 1;
else
    dd = t;
end
if t<2
    sigma_new1 = diag(ones(1,n)) - Delta/G*Delta';
%sigma_new1 = D'*inv(D*D'+diag(ones(1,dd)))*D;
     V0 = mvnrnd(zeros(1,n),sigma_new1,n_size);
     pd = makedist('Normal');
     s =  truncate(pd,0,inf);
     V1 = random(s,n_size,1);
     V1 = (Delta/G)*V1';
     V1 = V1';
     out = V0+V1;
else

sigma_new1 = diag(ones(1,n)) - Delta/G*Delta';
%sigma_new1 = D'*inv(D*D'+diag(ones(1,dd)))*D;
V0 = mvnrnd(zeros(1,n),sigma_new1,n_size);

%sigma_new2 = Sinv*(D*D'+diag(ones(1,dd)))*Sinv;
sigma_new2 = G;

l = zeros(1,dd);
u = repelem(Inf,dd);

V1 = mvrandn(l,u,sigma_new2,n_size);
%V1 = D'/(D*D'+diag(ones(1,dd)))*S*V1;
V1 = (Delta/G)*V1;
V1 = V1';
out = V0+V1;
end
end