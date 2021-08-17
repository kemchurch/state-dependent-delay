function [a,b,X,Operator] = VerifyInterp(y,M)
% Computes linear interpolants yi=ai(x-xi) + bi on the grid
% [0,1/M,2/M,...,1] for xi = (i-1)/M. Generates a "fast" interpolation
% operator mapping y->X = [a;b]. 
b = y(1:M,:);
%if isintval(y)
%    a = (y(2:M+1,:)-y(1:M,:))*intval(M);
%else
    a = (y(2:M+1,:)-y(1:M,:))*M;
%end
X = [a;b];
Operator = -eye(M,M+1) + [zeros(M,1),eye(M,M)];
%if isintval(y)
%    Operator = intval(Operator);
%end
Operator = [M*Operator; [eye(M,M),zeros(M,1)]];
end