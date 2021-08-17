function [fval,sM_grid] = EvaInterp(X,MSUB)
% Evaluates a linear interpolant X on M+1 gridpoints [0,1/M,...,1] on a
% submesh such that [j/M,(j+1)/M] is given MSUB further sub-points.
% Evaluates at all sub-points. Compatible with INTLAB.
% Inputs: X (double / intval) having 2*M rows. MSUB (int64). 
% Outputs: fval, value of interplant. sM_grid, points at which interplant
% is evaluated.
% Note: Code is vectorized such that F([X1,X2,...,Xk]) = [F(X1),...,F(Xk)] 
% for F(::) = evalinterp(::,M,MSUB,SIDE).
M = size(X,1)/2;
a = X(1:M,:);
b = X(M+1:2*M,:);
fval = zeros(M*MSUB,1).';
if isintval(X)
    sgrid = (1:MSUB).'/(intval(M*(MSUB+1)));
    Mgrid = (0:M-1).'/intval(M);
    fval = intval(fval);
else
    sgrid = (1:MSUB).'/(M*(MSUB+1));
    Mgrid = (0:M-1).'/M;
end
x_grid = kron(Mgrid,ones(MSUB,1));
sM_grid = reshape((sgrid.'+Mgrid).',[MSUB*M,1]);
a_grid = kron(a,ones(MSUB,1));
b_grid = kron(b,ones(MSUB,1));
fval = a_grid.*(sM_grid-x_grid) + b_grid;
end