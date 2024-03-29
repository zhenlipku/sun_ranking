clear
N = 50;
n = N;
trials = 25;
T = 200;
accuracy_sun = zeros(trials,T);

%% Get the pairwise structure of the data

all_pair = dlmread('data\all_pair.txt'); % load comparison data
score = dlmread('data\doc_info.txt'); % load ground-truth score

s = 1:1:n;   % select the first 50 items
score = score(:,2);
score = score(s);
for l = 1:length(all_pair(:,1))
    if ismember(all_pair(l,2),s) && ismember(all_pair(l,3),s)
        ;
    else
        all_pair(l,1) = 0;
    end
end
all_pair = all_pair(all_pair(:,1) ~= 0,:);
all_pair(:,1) = [];
data_ref = all_pair;

for s = 1:trials
    rng(12345+s,'philox')  % set random seed
    all_pairs = combnk(1:n, 2); % generate all the possible pairs for ranking
    data = all_pairs;
    if data(1) ~= 1
        data =  flipud(data);
    end
    all_pairs_D = zeros(n*(n-1)/2,n);
    for ll = 1:n*(n-1)/2
        all_pairs_D(ll,data(ll,1)) = 1;
        all_pairs_D(ll,data(ll,2)) = -1;
    end
    
    data = data_ref;
    info  = zeros(n,n);
    data_temp = data;
    for ii = 1:length(data_temp(:,1))
        if data_temp(ii,1) > data_temp(ii,2)
            temp_value = data_temp(ii,1);
            data_temp(ii,1) = data_temp(ii,2);
            data_temp(ii,2) = temp_value;
        end
    end
    
    for ii = 1:length(data_temp(:,1))
        info(data_temp(ii,1),data_temp(ii,2)) = info(data_temp(ii,1),data_temp(ii,2))+1;  
    end
    t = 1;
    x = [];
    y = [];
    m = length(data_ref(:,1));
    
    init = randsample(m,1); % randomly select a pair for comparison initially
    i = data(init,1);
    j = data(init,2);
    info(i,j) = info(i,j) - 1;
    data(init,:) = [0 0];
    
    x_t = zeros(1,n);
    x_t(i) = 1;
    x_t(j) = -1;
    x = [x;x_t];
    y_t = 1;
    y = [y,y_t];
    D = diag(2*y-1)*x./sqrt(2); % generate parameters of the SUN distribution
    s_vec = zeros(t,1);
    for ii = 1:t
        s_vec(ii) = sqrt(D(ii,:)*D(ii,:)'+1);
    end
    Sinv = diag(1./s_vec);
    Delta = D'*Sinv;
    G = Sinv*(D*D'+diag(t))*Sinv;
    M1 = sqrt(2/pi)*Delta*ones(1,1); % the posterior with only one observation is skew normal, the mean of which is as shown.
    dist = 0;   % compute the kendall's tau
    for i = 1:n
        for j = i+1:n
            if score(i) == score(j)
                dist = dist+1;
            elseif  (score(i) - score(j))*(M1(i) - M1(j)) > 0
                dist = dist + 1;
            end
        end
    end
    accuracy_sun(s,t) = 2 * dist / (n * (n-1));
    for t = 2:T
        
        para_temp = est_para_mf(x,y);  % using MFVI to estimate the posterior
        Mean = para_temp.Mean;
        Sigma = para_temp.Sigma;
        outcomes = zeros(n,n);
        
        tic
        for i = 1:n
            for j = i+1:n
                if info(i,j)>0
                    est_trans_prob = est_tran_prob_compute(Mean,Sigma,i,j,n);
                    outcomes(i,j) = est_trans_prob*MF_KG(x,y,i,j,n,all_pairs_D)+(1-est_trans_prob) * MF_KG(x,y,j,i,n,all_pairs_D);
                end
            end
        end
        toc
        
        
        [maxOutcome,ind] = max(outcomes(:));
        [xx,yy] = ind2sub(size(outcomes),ind);
        Lia = ismember(data, [xx yy;yy xx],'rows');
        
        k = find(Lia);
        if sum(Lia) > 1
            r = randsample(k,1);
        else
            r = k;
        end
        
        i = data(r,1);
        j = data(r,2);
        
        y_t = 1;
        x_t = zeros(1,n);
        x_t(i) = 1;
        x_t(j) = -1;
        x = [x;x_t];
        
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
        [nrow,~] = size(D);
        
        if rem(t,10) == 0
            out = est_para_mf(x,y);
            out_score = out.Mean;
            output.logic = 1;
            
            dist = 0;   % compute the kendall's tau
            for i = 1:n
                for j = i+1:n
                    if score(i) == score(j)
                        dist = dist+1;
                    elseif  (score(i) - score(j))*(out_score(i) - out_score(j)) > 0
                        dist = dist + 1;
                    end
                end
            end
            accuracy_sun(s,t) = 2 * dist / (n * (n-1));
            %accuracy_sun(s,t) = corr(out_score,theta,'type','Kendall');
            fprintf('Iter=%d, logic = %d ,ACCURACYT1=%f\n', t,output.logic,accuracy_sun(s,t))
        else
            out = est_para_mf(x,y);
            out_score = out.Mean;
            dist = 0;   % compute the kendall's tau
            for i = 1:n
                for j = i+1:n
                    if score(i) == score(j)
                        dist = dist+1;
                    elseif  (score(i) - score(j))*(out_score(i) - out_score(j)) > 0
                        dist = dist + 1;
                    end
                end
            end
            accuracy_sun(s,t) = 2 * dist / (n * (n-1));
        end
        
    end
end