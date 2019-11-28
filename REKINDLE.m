function [res, DWIB0, DT, varargout] = REKINDLE(DWI, bval, g, mask, DKI, par)

% REKINDLE algorithm as published in Tax et al., MRM 2015
% Copyright Chantal Tax (taxc@cardiff.ac.uk) and Alexander Leemans
% (alexander@isi.uu.nl)

% Input
% DWI: X x Y x Z x N matrix with DWI signals
% bval: N x 1 matrix with b-values
% g: N x 3 matrix with gradient unit directions
% mask: X x Y x Z matrix with mask
% DKI: 0 for DTI fit, 1 for DKI fit
% par: has fields 
%      con: convergence value
%      iter1: maximum iterations inner loop
%      iter2: maximum iterations outer loop

% Output
% res: X x Y x Z x N matrix with residuals
% DWIB0: X x Y x Z matrix estimated b=0 signal 
% DT: cell with the diffusion tensor elements
% varargout{1}: cell with the kurtosis tensor elements

% Note that in the output estimates outliers are downweigthed by the iterative
% procedure, but not removed

Vm = sum(mask(:));
DWI_m = vec(DWI);
DWI_m(DWI_m<=0) = eps;

b = repmat(bval,[1 6]).*[g(:,1).^2 2*g(:,1).*g(:,2) 2*g(:,1).*g(:,3) ...
    g(:,2).^2 2*g(:,2).*g(:,3) g(:,3).^2];
    
if DKI 
    b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
        4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
        4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
        4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
        6*(g(:,1).^2).*(g(:,2).^2) ...
        6*(g(:,1).^2).*(g(:,3).^2) ...
        6*(g(:,2).^2).*(g(:,3).^2) ...
        12*g(:,2).*g(:,3).*(g(:,1).^2) ...
        12*g(:,1).*g(:,3).*(g(:,2).^2) ...
        12*g(:,1).*g(:,2).*(g(:,3).^2)];
    b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;
    b_final = [ones(size(b,1),1) -b b_kurt];
    X = zeros(22,Vm);
else
    b_final = [ones(size(b,1),1) -b];
    X = zeros(7,Vm);
end
res_m = nan(size(DWI_m));

parfor i=1:Vm
    [res_m(:,i), X(:,i)] = Robust_Linear_Reweighting(DWI_m(:,i),b_final,par);
end

res = unvec(res_m,mask);

DWIB0 = nan(size(mask));
DWIB0(mask) = exp(X(1,:));

dummy = nan(size(mask));
for i = 1:6
    DT{i} = dummy;
    DT{i}(mask) = X(i+1,:);
end

if DKI
    for i = 7:21
        KT{i-6} = dummy;
        KT{i-6}(mask) = X(i+1,:)'./((DT{1}(mask)+DT{4}(mask)+DT{6}(mask)).^2);
    end
    varargout{1} = KT;
end

dwib0 = DWIB0(mask);
y = max(dwib0(~isinf(dwib0)));
dwib0(dwib0>y) = y;
DWIB0(mask) = dwib0;

end