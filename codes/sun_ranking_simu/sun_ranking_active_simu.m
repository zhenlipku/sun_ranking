clear
n = 10; % number of items
trials = 25; % number of trials
alpha0 = ones(1,n);
T = 80; % the total budget
accuracy_sun = zeros(trials,T); % store the kendall's tau
for s = 1:trials
    s
    rng(12345+s,'philox')   %set random seed
    t = 1;
    x = [];
    y = [];
    %score = drchrnd(alpha0,1); % dirichlet distributed ground-truth score
    score =2*n*rand(1,n);  % normal distributed ground-truth score
    theta = score;
    
    all_pairs = combnk(1:n, 2); % all the possible pairs for comparison
    [m,~] = size(all_pairs);
    data = all_pairs;
    if data(1) ~= 1
        data =  flipud(data);
    end
    
    init = randsample(m,1); % randomly select a pair for comparison initially
    i = data(init,1);
    j = data(init,2);
    
    
    x_t = zeros(1,n);
    x_t(i) = 1;
    x_t(j) = -1;
    x = [x;x_t];
    
    r = rand; % start comparison
    % split = theta(i)/(theta(i)+theta(j));  % BTL model
    split = normcdf((score(i)-score(j))) ;  % TM model
    if r <= split
        y_t = 1;
    else
        y_t = 0;
    end
    y = [y,y_t];  % store the comparison outcomes
    
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
    accuracy_sun(s,t) = 2*accuracy_sun(s,t) - 1;
    
    for t = 2:T
        if t<26
            Denominator0 = mvncdf(zeros(1,t-1),zeros(1,t-1),G); % the denominator of the corresponding SUN distribution
        else
            Denominator0=   mvncdf_fast(-Inf*ones(1,t-1),zeros(1,t-1),G,1000).prob;  % when the dimension of multivariate normal CDF larger than 25
        end
        outcomes = zeros(m,1);  %Pairs that never participate in comparison are expected to yield the same results in the KG policy.
        for k = 1:m
            if  sum(sum(Delta(data(k,:),:) == 0)) == 2*(t-1)
                i = data(k,1);
                j = data(k,2);
                trans_prob = tran_prob_compute(i,j,x,y)/Denominator0;
                if  trans_prob >1
                    trans_prob = 1;
                end
                outcomes_special = trans_prob*max_prob(x,y,i,j,data)+ (1-trans_prob) *max_prob(x,y,j,i,data);
                break
            end
        end
        tic
        parfor k = 1:m
            i = data(k,1);
            j = data(k,2);
            if sum(sum(Delta(data(k,:),:) == 0)) == 2*(t-1)
                outcomes(k) = outcomes_special;  %Pairs that never participate in comparison are expected to yield the same results in the KG policy.
            else
                trans_prob = tran_prob_compute(i,j,x,y)/Denominator0;
                if  trans_prob >1
                    trans_prob = 1;
                end
                outcomes(k) = trans_prob*max_prob(x,y,i,j,data)+ (1-trans_prob) *max_prob(x,y,j,i,data);
            end
        end
        toc
        
        maxval = max(outcomes);
        rows = find(outcomes == maxval);
        if length(rows) == 1
            row = rows;
        else
            row =  randsample(rows,1);
        end
        
        %row = randsample((n-1)*n/2,1);    % uniformly select a pair
        
        i = data(row,1);
        j = data(row,2);
        
        
        r = rand;
        %split = theta(i)/(theta(i)+theta(j));  % BTL model
        split = normcdf((score(i)-score(j))) ;  % TM model
        if r <= split
            y_t = 1;
        else
            y_t = 0;
        end
        
        % compute parameters of SUN distribution
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
        
        if rem(t,10) == 0  % Calculate Kendall's Tau for every 10 samples
            trials_temp = 12;  % Times for Monte Carlo integration
            rr = com_prob_mean(data,Delta,G,trials_temp);  % compute all the Pred_T() via monte carol integration
            
            % To verify if the upper bound of C_T will form a rank.
            com_matrix = zeros(n,n);
            oj_series = data(:,2);
            oi_series = data(:,1);
            com_matrix((oj_series-1)*n+oi_series) = rr.prob1;
            com_matrix((oi_series-1)*n+oj_series) = rr.prob2;
            
            %The condition logical_val == 0 does not guarantee the existence of a loop, as numerical computational errors may occur when $\mathrm{Pred}_T(i>j)$ theoretically equals 0.5.
            
            output = rank_ver(com_matrix);    % ranking verification
            if output.logical_val == 1
                out_score = output.score;
            else
                out = generateG(x,y);  % the function to generate parameters of a SUN distribution with X_T and y_T
                out_score = rr.mean';
            end
            
            output.logic = output.logical_val;
            % output.logic = output.logic;
            accuracy_sun(s,t) = corr(out_score,theta','type','Kendall');
            fprintf('Iter=%d, logic = %d ,ACCURACYT1=%f\n', t,output.logic,accuracy_sun(s,t))
        end
        
        
    end
end
