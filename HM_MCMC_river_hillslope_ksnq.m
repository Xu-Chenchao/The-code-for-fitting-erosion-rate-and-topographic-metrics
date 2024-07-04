function [params,like_data,mMAP,accrate] = HM_MCMC_river_hillslope_ksnq(data,param_start,priors,burn_in,n_iter)

% code will predict mean slope, ksn and ksnq given an erosion rate (E),
    % stream power erodibility (K) and precip modified erodibility (Kp =
    % K/(p^m)), drainage area exponent (m), slope exponent(n), critical
    % slope angle (Sc) diffusivity (D), and ratio of rock to soil density
    % (B).
    %
    % Note:
    % E needs to be in m/yr
    % Sc in m/m
    %
    % Author: Sean F. Gallen
    % contact: sean.gallen@colostate.edu
    % date modified: 10.09.2023

    % unpack data and uncertainty
    E = data(:,1);
    S = data(:,2);
    Ser = data(:,3);
    ksn = data(:,4);
    ksner = data(:,5);
    ksnq = data(:,6);
    ksnqer = data(:,7);


    % unpack parameters
    % free parameters
    K = 10^param_start(1);
    %Kp = 10^param_start(2);
    n = param_start(2);
    Sc = param_start(3);
    D = 10^param_start(4);
    current = param_start(1:4);

    % fixed parameters
    m = n*param_start(5);
    B = param_start(6);

    % steplength (currently based on priors
    p_steps = range(priors,2)'.*5e-3;
    
    % initialize parameters matix
    params = nan(burn_in+n_iter,4);

    % allocate memory to catch likelihoods, acceptance rate
    like_data = nan(burn_in+n_iter,3);

    % make sure new matlab session use different random numbers
    rng shuffle

    % make initial model
    params(1,:) = current'; 

    % run initial model
    [S_mod,ksn_mod,ksnq_mod] = river_nonlinear_hillslope_forward_model(E,K,K,m,n,Sc,D,B);

    % stack and observed and modeled data and errors
    obs = [ksnq;S];
    o_err = [ksnqer;Ser];
    mod = [ksnq_mod;S_mod];

    % calculate the log-likelihood of the model given the data
    lpcurrent = logprior_uniform(params(1,:),priors);
    llcurrent = mod_loglikelihood(obs,mod,o_err);

    current = params(1,:);
    lMAP = -Inf;
    mMAP = current;
    nacc = 0;

   h = waitbar(0,'Please wait...');
   for i = 2:burn_in+n_iter
        % pick candidate parameters
        candidate = current+p_steps.*randn(1,4);

        % update free parameters
        K = 10^candidate(1);
        %Kp = 10^candidate(2);
        n = candidate(2);
        Sc = candidate(3);
        D = 10^candidate(4);
        m = n*param_start(5);

        % run model with candidate parameters
        [S_mod,ksn_mod,ksnq_mod] = river_nonlinear_hillslope_forward_model(E,K,K,m,n,Sc,D,B);

        % stack and observed and modeled data and errors
        obs = [ksnq;S];
        o_err = [ksnqer;Ser];
        mod = [ksnq_mod;S_mod];

        % calculate model, priors and step for log of acceptance ratio
        lpcandidate = logprior_uniform(candidate,priors);
        llcandidate = mod_loglikelihood(obs,mod,o_err);
        
        % transition probabilities
        lr1 = logproposal(candidate,current,p_steps);
        lr2 = logproposal(current,candidate,p_steps);

        % add all the probabilities together for the acceptance ratio
    logalpha = lpcandidate + llcandidate + lr1 - lpcurrent - llcurrent - lr2;
    
    % Take the minimum of the log(alpha) and 0.
    if (logalpha > 0)
        logalpha = 0;
    end
    
    % Generate a U(0,1) random number and take its logarithm.
    logt = log(rand());
    
    % Accept or reject the step.
    if (logt < logalpha)
        
        % Accept the step.
        current = candidate;
        %if i > burn_in
            nacc = nacc + 1;
        %end
        
        % Update the MAP solution if this one is better.
        if ((lpcandidate + llcandidate) > lMAP)
            lMAP = lpcandidate + llcandidate;
            mMAP = candidate;
        end
        
        lpcurrent = lpcandidate;
        llcurrent = llcandidate;
      
    else
        % reject the step
    end

        % store the chain paths and other data
    params(i, :) = current;
    like_data(i,1) = llcurrent;
    like_data(i,2) = lpcurrent; 
    %  
    %if i > burn_in
        accrate = nacc / (burn_in+n_iter);
        like_data(i,3) = nacc/(i)*100;
    %end
    waitbar(i/(burn_in+n_iter),h)
   end
    close(h)
end

%% likelihood functions used in the model
function l = mod_loglikelihood(obs,mod,sig)
    err = obs-mod;
    l = (-1/2).*sum((err./sig).^2);
end

function lp = logprior_uniform(candidate,prior_bounds)
    
    log_test = zeros(size(candidate));
    for i = 1:length(candidate)
        if candidate(i) >= prior_bounds(i,1) &&  candidate(i) <= prior_bounds(i,2)
            log_test(i) = 1;
        end
    end
    
    if sum(log_test) == length(candidate)
        lp=0;
    else
        lp = -inf;
    end
end


function [lr] = logproposal(p_1,p_2,step)
    lr = (-1/2)*sum((p_1-p_2).^2./step.^2);
end
