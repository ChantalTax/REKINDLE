function [res, X] = Robust_Linear_Reweighting(S0,B,par)

% REKINDLE reweighting as published in Tax et al., MRM 2015
% Copyright Chantal Tax (taxc@cardiff.ac.uk) and Alexander Leemans
% (alexander@isi.uu.nl)

% Initialization
crit_all = 1;
cn_all = 0; 
iter_max_in = par.iter1;
iter_max = par.iter2;
con = par.con;

% Step 1: Initial LLS fit
w = ones(size(S0));
LogS0 = log(S0);
W = diag(w);
X = (B'*W*B)\(B'*W)*LogS0;

while crit_all==1
    
    % Initialization
    cn_all = cn_all + 1;
    X_p_all = X;
    cn = 0;
    crit = 1;
    
    % Step 2: Compute a robust estimate for homoscedastic regression using IRLS.
    while crit==1

        % Initialization
        cn = cn+1;
        X_p = X;
        
        fit = B*X; 
        % a. Calculate the residuals e in the linear domain
        res = (LogS0-fit); 
        if all(res==0)
            break;
        end
        % b. Obtain an estimate of the dispersion of the residuals by
        % calculating the median absolute deviation (MAD).
        C = 1.4826*median(abs(res-median(res))) ;
        % c. Recompute the weights according to Eq. [13].
        w = 1./(1 + (res/C).^2).^2;
        W = diag(w);
        % d. Perform WLLS fit with new weights
        X = (B'*W*B)\(B'*W)*LogS0;
        % e. Check convergence
        if all(abs(X-X_p) <= con*max(abs(X),abs(X_p))) || cn == iter_max_in
            crit=0;
        end

    end

    % Step 3: Transform variables for heteroscedasticity
    fit = B*X;
    LogS0_2 = LogS0./(exp(-fit));
    B2 = B./repmat((exp(-fit)),[1,size(B,2)]);

    % Initialization
    crit = 1;
    cn = 0;
    
    % Step 4: Initial LLS fit in * domain
    w = ones(size(S0));
    W = diag(w);
    X = (B2'*W*B2)\(B2'*W)*LogS0_2;  %

    % Step 5: Compute a robust estimate for homoscedastic regression using IRLS.
    while crit==1

        % Initialization
        cn = cn+1;
        X_p = X;
        
        fit = B2*X; 
        % a. Calculate the residuals e* in the linear domain
        res = (LogS0_2-fit); 
        if all(res==0)
            break;
        end
        % b. Obtain an estimate of the dispersion of the residuals by
        % calculating the median absolute deviation (MAD).
        C = 1.4826*median(abs(res-median(res))) ; 
        % c. Recompute the weights according to Eq. [13].
        w = 1./(1 + (res/C).^2).^2;
        W = diag(w);

        % d. Perform WLLS fit with new weights
        X = (B2'*W*B2)\(B2'*W)*LogS0_2;
        % e. Check convergence
        if all(abs(X-X_p) <= con*max(abs(X),abs(X_p))) || cn == iter_max_in
            crit=0;
        end

    end

    %  Step 6: Check convergence overall loop
    if all(abs(X-X_p_all) <= con*max(abs(X),abs(X_p_all))) || cn_all == iter_max
            crit_all=0;
    end
end

% Step 7: Identify and exclude outliers
fit = B2*X; % In the first iteration, this is the first fit
res = (LogS0_2-fit);




