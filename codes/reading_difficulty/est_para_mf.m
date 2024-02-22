function est_para = est_para_mf(x,y)

% This is a matlab version of of the mean-field variational Bayes implementation 
% available at https://github.com/augustofasano/Probit-PFMVB.

tolerance = 10^-2;
maxIter = 10^4;
[m,n] = size(x);

if n <= m
    Omega = eye(n);
    invOmega = Omega;
    V = (eye(n)/(x'*x+invOmega));
    diagV = diag(V);
    VXt = V*x';
    H = x*VXt;
else
    XXt = x*x';
    invOmZ = eye(m)/(eye(m)+XXt);
    VXt = x'*invOmZ; 
    diagV = 1 - sum(VXt'.*x);
    H = XXt*invOmZ;
    V = (eye(n)/(x'*x+eye(n)));
end

XVVXt = VXt'*VXt ;
signH = H;
signH(y==0,:) = -signH(y==0,:);

meanZ = zeros(m,1);
elbo = -Inf;
diff = 1;
nIter=0;

while diff > tolerance && nIter < maxIter
    elboOld = elbo;
    mu = H*meanZ;
    meanZ = mu +(2*y-1)'.* (normpdf(mu)./((1-erf(-(2*y-1)'.*mu./sqrt(2)))./2) );
    elbo = -0.5*(meanZ'*XVVXt*meanZ) + sum(log((1-erf(-signH*meanZ./sqrt(2)))./2));
   % elbo = -(meanZ'*invOmZ*meanZ - sum((meanZ.^2)'.*diagInvOmZ))/2 - sum(meanZ'.*mu./sigma2) + sum((mu.^2)./sigma2)/2 + sumLogPhi;
    diff = abs(elbo-elboOld);
    nIter = nIter+1;
end

est_para.Mean = VXt*meanZ;
est_para.Sigma = V;
 

end

 