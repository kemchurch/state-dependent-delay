function F = function_F_Stumpf(X,M,MSUB,phi,eta)
u = X(1:end-1); tau = X(end);
% \psi component
Fpsi = zeros(M+1,1);
Fpsi_integral = zeros(M,1);
[~,~,U] = VerifyInterp(u,M);
[Su,s] = EvaInterp(U,MSUB);
y = f(Su,eval1spline(phi,h(tau*s,Su,eta)));
Fpsi(1) = -u(1) + eval1spline(phi,0);
for j=1:M
    yj = y(1+(j-1)*MSUB:j*MSUB);
    sj = [(j-1)/M,j/M];
    Fpsi_integral(j) = Quad(yj,sj);
    Fpsi(j+1) = - u(j+1) + eval1spline(phi,0) + tau*sum(Fpsi_integral(1:j));
end
% \R component
FR = h(tau,u(end),eta);
% Gather components.
F = [Fpsi;FR];
end

function f = f(x,y)
global alpha
f = alpha*(x-y) - x.*(abs(x));
end

function h = h(t,x,eta)
global mu
h = t - exp(-abs(x).*sin(eta*x).^2 - mu*x.^2);
end