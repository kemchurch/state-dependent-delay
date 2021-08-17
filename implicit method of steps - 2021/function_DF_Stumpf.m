function DF = function_DF_Stumpf(X,M,MSUB,phi,dphi,eta)
global mu
u = X(1:end-1); tau = X(end);
% \psi component
% Du_k
DFpsi_u = zeros(M+1,M+1);
DFpsi_u_integrals = zeros(M,M+1);
[~,~,U,Mat_S] = VerifyInterp(u,M);
[Su,s] = EvaInterp(U,MSUB);
I = eye(M+1,M+1);
SI = EvaInterp(Mat_S*I,MSUB);
Df1 = D1f(Su,eval1spline(phi,h(tau*s,Su,eta)));
Df2 = D2f(Su,eval1spline(phi,h(tau*s,Su,eta)));
Dh1 = D1h(tau*s,Su);
Dh2 = D2h(tau*s,Su,eta);
for j=1:M
    Suj = Su(1+(j-1)*MSUB:j*MSUB);
    SIj = SI(1+(j-1)*MSUB:j*MSUB,:);
    D_1fj = Df1(1+(j-1)*MSUB:j*MSUB);
    D_2fj = Df2(1+(j-1)*MSUB:j*MSUB);
    %D_1hj = Dh1(1+(j-1)*MSUB:j*MSUB);
    D_2hj = Dh2(1+(j-1)*MSUB:j*MSUB);
    hj = h(tau*s(1+(j-1)*MSUB:j*MSUB),Suj,eta);
    sj = [(j-1)/M,j/M];
    DFpsi_u_integrals(j,:) = Quad(...
        D_1fj.*SIj + D_2fj.*eval1spline(dphi,hj).*D_2hj.*SIj, sj);
    DFpsi_u(j+1,:) = tau*sum(DFpsi_u_integrals(1:j,:),1);
end
DFpsi_u = DFpsi_u - eye(M+1,M+1);
% Dtau
y = f(Su,eval1spline(phi,h(tau*s,Su,eta)));
DFpsi_tau = zeros(M+1,1);
DFpsi_tau_integrals = zeros(M,1);
for j=1:M
    yj = y(1+(j-1)*MSUB:j*MSUB);
    Suj = Su(1+(j-1)*MSUB:j*MSUB);
    D_2fj = Df2(1+(j-1)*MSUB:j*MSUB);
    D_1hj = Dh1(1+(j-1)*MSUB:j*MSUB);
    hj = h(tau*s(1+(j-1)*MSUB:j*MSUB),Suj,eta);
    sj = [(j-1)/M,j/M];
    DFpsi_tau_integrals(j) = Quad(...
        yj + tau*D_2fj.*eval1spline(dphi,hj).*D_1hj.*s(1+(j-1)*MSUB:j*MSUB),sj );
    DFpsi_tau(j+1) =  sum(DFpsi_tau_integrals(1:j));
end
% \R component
DFR_u = zeros(1,M+1);
DFR_u(M+1) = D2h(1,u(end),eta);
DFR_R = D1h(1,u(end));
% Gather
DF = [DFpsi_u,DFpsi_tau];
DF = [DF;[DFR_u,DFR_R]];
end

function f = f(x,y)
global alpha
f = alpha*(x-y) - x.*(abs(x));
end

function h = h(t,x,eta)
global mu
h = t - exp(-abs(x).*sin(eta*x).^2 - mu*x.^2);
end

function D1f = D1f(x,y)
global alpha
D1f = alpha - 2*abs(x);
end

function D2f = D2f(x,y)
global alpha
D2f = -alpha+0*x;
end

function D1h = D1h(t,x)
D1h = 1 + 0*x;
end

function D2h = D2h(t,x,eta)
global mu
D2h = ((sign(x).*sin(eta*x).*(sin(eta*x)+2*eta*x.*cos(eta*x))) + 2*mu*x).*exp(-abs(x).*sin(eta*x).^2 - mu*x.^2);
end