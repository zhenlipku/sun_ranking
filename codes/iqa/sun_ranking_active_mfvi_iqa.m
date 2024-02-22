clear
%rng(123,'philox') %set random seed
n = 16; % number of items
trials = 25;
T = 120; % the total budget
accuracy_sun = zeros(trials,T);
ref_num=15;
for ref = 1:ref_num
    str=strcat('PC-IQA dataset\','data',num2str(ref),'.mat');
    load(str) % load data
    score=Hodgerank(data_ref)'; % Compute the estimated ground-truth score using HodgeRank
    data_ref_temp=data_ref;
    
    Z = zeros(n,n);
    for i=1:length(data_ref)
        Z(data_ref(i,1),data_ref(i,2)) = Z(data_ref(i,1),data_ref(i,2))+1;
        data_ref_temp(i,1)=min(data_ref(i,:));
        data_ref_temp(i,2)=max(data_ref(i,:));
    end
    edge=unique(data_ref_temp,'rows');
    w = zeros(1,length(edge));
    p0 = w;
    for i=1:length(edge)
        w(i) = Z(edge(i,1),edge(i,2)) + Z(edge(i,2),edge(i,1));
        p0(i) = Z(edge(i,1),edge(i,2))/w(i);
    end
    for s = 1:trials
        
        rng(12345+s,'philox')    % set random seed
        %l1 = zeros(1,n*(n-1)/2);
        %l2 = zeros(1,n*(n-1)/2);
        t = 1;
        x = [];
        y = [];
        
        theta = score;
        all_pairs = combnk(1:n, 2); % all the possible pairs for comparison
        [m,~] = size(all_pairs);
        data = all_pairs;
        if data(1) ~= 1
            data =  flipud(data);
        end
        all_pairs_D = zeros(n*(n-1)/2,n);
        for ll = 1:n*(n-1)/2
            all_pairs_D(ll,data(ll,1)) = 1;
            all_pairs_D(ll,data(ll,2)) = -1;
        end
        init = randsample(m,1); % randomly select a pair for comparison initially
        i = data(init,1);
        j = data(init,2);
        
        x_t = zeros(1,n);
        x_t(i) = 1;
        x_t(j) = -1;
        x = [x;x_t];
        
        r = rand; % start comparison
        
        split = p0(init) ;
        if r <= split
            y_t = 1;
            %       l1(init) = l1(init)+1;
        else
            y_t = 0;
            %        l2(init) = l2(init)+1;
        end
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
                if (theta(i) - theta(j))*(M1(i) - M1(j)) > 0
                    dist = dist + 1;
                end
            end
        end
        accuracy_sun(s,t) = 2 * dist / (n * (n-1));
        
        for t = 2:T
            para_temp = est_para_mf(x,y); % using MFVI to estimate the posterior
            Mean = para_temp.Mean;
            Sigma = para_temp.Sigma;
            
            outcomes = zeros(m,1); % carry out the KG method
            
            parfor k = 1:m
                i = data(k,1);
                j = data(k,2);
                est_trans_prob = est_tran_prob_compute(Mean,Sigma,i,j,n);  % the estimated transitional probability
                outcomes(k) = est_trans_prob*MF_KG(x,y,i,j,n,all_pairs_D)+(1-est_trans_prob) * MF_KG(x,y,j,i,n,all_pairs_D);
                
            end
            
            
            maxval = max(outcomes);
            rows = find(outcomes == maxval);
            if length(rows) == 1
                row = rows;
            else
                row = randsample(rows,1);
            end
            
            i = data(row,1);
            j = data(row,2);
            
            
            r = rand;
            %split = theta(i)/(theta(i)+theta(j));
            split = p0(row) ;
            if r <= split
                y_t = 1;
                %l1(row) = l1(row)+1;
            else
                y_t = 0;
                %l2(row) = l2(row)+1;
            end
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
            %M1 = generatemean(Delta,G);
            
            if rem(t,10) == 0  % Calculate Kendall's Tau for every 10 samples
                
                out = est_para_mf(x,y);
                out_score = out.Mean;
                output.logic = 1;
                accuracy_sun(s,t) = corr(out_score,theta','type','Kendall');
                fprintf('Iter=%d, logic = %d ,ACCURACYT1=%f\n', t,output.logic,accuracy_sun(s,t))
            else
                out = est_para_mf(x,y);
                out_score = out.Mean;
                accuracy_sun(s,t) = corr(out_score,theta','type','Kendall');
            end
            
            
            
        end
    end
    save(strcat('sun_active_iqa_trials25_rng12345_refnum_',num2str(ref),'.mat'))
end
