function [conv,cons,X] = Newton(phi,dphi,M,MSUBS_F,MSUBS_DF,tol,maxiter,seedmax,size,print,specifycrossing,specifysolution)
global pushInf eta
rseed = 0;
while (rseed < seedmax)
    X = size.*randn(M+2,1);
    if nargin >= 11
        if ~isempty(specifycrossing)
            X(end) = specifycrossing;
        end
    end
    if nargin == 12
        if ~isempty(specifysolution)
            X(1:end-1) = specifysolution;
        end
    end
    %X(end)=0.5002;
    F = function_F_Stumpf(X,M,MSUBS_F,phi,eta);
    iters = 0;
    while (norm(F)>tol) && (iters<maxiter) && norm(F)<100
        X = X - function_DF_Stumpf(X,M,MSUBS_DF,phi,dphi,eta)\F;
        if sum(isnan(X))>0
            error('There is a NaN.');
        end
        F = function_F_Stumpf(X,M,MSUBS_F,phi,eta);
        iters = iters + 1;
        if mod(iters,20)==0 && print==1
            disp(['Iteration ',num2str(iters),', seed ',num2str(rseed), ', Newton residual: ',num2str(norm(F)),'.']);
        end
    end
    rseed = rseed + 1;
    if print==1 && norm(F)>tol && mod(rseed,100)==0
        disp(['Bad seed. Re-seeding. (',num2str(rseed),').']);
        size = size + 0.1*(-1 + 2*rand(1));
    end
    if norm(F)<=tol
        break
    end
end
if norm(F)<=tol
    disp(['Newton converged to within tolerance after ',num2str(iters-1),' iterates.']);
    disp(['Final residual: ',num2str(norm(F)),'.']);
    conv = 1;
    % Constraint check.
    pushInf = 1;
    F = function_F_Stumpf(X,M,MSUBS_F,phi,eta);
    if sum(isinf(F))>0
        cons = 0;
    else
        cons = 1;
    end
    pushInf = 0;
else
    error(['Newton did not converge.']);
end

end