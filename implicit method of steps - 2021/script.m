clc
clf
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:deval:NonuniqueSolution')
% Model parameter; adjust asd needed.
global alpha eta mu
alpha = 3;
eta = 6; mu = 0.5;
%eta = 10;mu = 0.1;
% Integration parameters
use_data = 2; 
talkative = 1;
N_STEPS = 14;
M = 40;
MSUBS_F = 100;
MSUBS_DF = 10;
tol = 5E-15;
max_iter = 1000;
max_seed = 1;   
max_crossfail = 1;
% Define the initial condition
dom_PHI = [-1;0];   % Can be adjusted if need be.
phi = @(x)2 + 0*x;  % Note: include 0*x for function overloading.
dphi = @(x)0*x;     % See above.
% Get guesses for crossing times and initial data - if applicable
[data,get_taus] = generate_data(@(t,x,y)f(x,y),@(t,y)h(t,y),phi,N_STEPS);
if use_data == 0
    tau_guess = []; X_guess = [];
elseif use_data == 1
    tau_guess = [get_taus(1);diff(get_taus)];   X_guess = [];
elseif use_data == 2
    tau_guess = [get_taus(1);diff(get_taus)];
    get_taus = [0;get_taus];
    X_guess = zeros(M+1,N_STEPS);
    for j=1:N_STEPS
        X_guess(:,j) = deval(data,linspace(get_taus(j),get_taus(j+1),M+1)).';
    end
end
% Interpolate the initial condition
subs = M+1;
[a_phi,b_phi,~,~] = VerifyInterp(phi(linspace(dom_PHI(1),dom_PHI(2),subs)).',subs-1);
[a_dphi,b_dphi,~,~] = VerifyInterp(dphi(linspace(dom_PHI(1),dom_PHI(2),subs)).',subs-1);
PHI = spline1(a_phi,b_phi,dom_PHI);
DPHI = spline1(a_dphi,b_dphi,dom_PHI);
% Begin
global pushInf
pushInf = 0;        % Used in Newton to check constraints.
X = zeros(M+2,N_STEPS);
tau = 0;
TAU = 0;
guess_scale = zeros(M+2,1);
SOL = PHI;
SOL_DERIVATIVE = DPHI;
for k=1:N_STEPS
    disp('---------------------');
    disp(['Computing step ',num2str(k),'.']);
    converged = 0;
    constraint = 0;
    Newton_attempt = 1;
    while ((constraint==0) | (converged==0)) & Newton_attempt <= max_crossfail
        if isempty(tau_guess)
            TAU_GUESS = [];
        else
            TAU_GUESS = tau_guess(k);
        end
        if isempty(X_guess)
            X_GUESS = [];
        else
            X_GUESS = X_guess(:,k);
        end
        [converged,constraint,X(:,k)] = Newton(PHI,DPHI,M,MSUBS_F,...
            MSUBS_DF,tol,max_iter,max_seed,guess_scale,1,TAU_GUESS,X_GUESS);
        Newton_attempt = Newton_attempt + 1;
        if converged==1 && constraint == 0
            disp('Constraints violated.');
        end
    end
    if (constraint==0) || (converged==0)
        error('Newton failed to converge or constraints are violated.');
    end
    guess_scale = abs(X(:,k));
    disp('Processing data.');
    % Get psi and dpsi on [0,1]
    [a_psi,b_psi,~,~] = VerifyInterp(X(1:M+1,k),M);
    PSI = spline1(a_psi,b_psi,[0;1]);
    psi = eval1spline(PSI,(0:M)'/M);
    dpsi = f(eval1spline(PSI,(0:M)'/M),...
        eval1spline(PHI,h(X(M+2,k)*(0:M)'/M,psi)));
    [a_dpsi,b_dpsi,~,~] = VerifyInterp(dpsi,M);
    % Update computed solution
    TAU = [TAU;TAU(end)+X(M+2,k)];
    [SOL,w] = fuse(SOL,spline1(a_psi,b_psi,[TAU(end-1);TAU(end)]));
    if w>1E-14
        disp(['Fused SOLUTION segments contain a gap of size ',num2str(w),'.']);
    end
    % Update derivative of solution
    if k==1
        errflag = 0;    % dphi is discontinuous on first iterate; expected.
    else
        errflag = 1;    % in general, dphi should be continuous.
    end
    [SOL_DERIVATIVE,w] = fuse(SOL_DERIVATIVE,...
        spline1(a_dpsi,b_dpsi,[TAU(end-1);TAU(end)]),errflag);
    if w>1E-14
        disp(['Fused SOLUTION DERIVATIVE segments contain a gap of size ',num2str(w),'.']);
    end
    if errflag == 0
        disp('(Expected; first segment)');
    end
    % Shift data back
    PHI = shift(SOL,-TAU(end));
    DPHI = shift(SOL_DERIVATIVE,-TAU(end));
end
warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:deval:NonuniqueSolution')

subplot(1,2,1)
plot(SOL)
subplot(1,2,2)
t = linspace(0,TAU(5),10000);
ht = h(t,eval1spline(SOL,t));
plot(t,ht,'k')
hold on
for j=1:4
    plot([TAU(1),TAU(j+1)],[TAU(j),TAU(j)],'r--');
end

function f = f(x,y)
global alpha
f = alpha*(x-y) - x.*(abs(x));
end

function h = h(t,x)
global mu eta
h = t - exp(-abs(x).*sin(eta*x).^2 - mu*x.^2);
end