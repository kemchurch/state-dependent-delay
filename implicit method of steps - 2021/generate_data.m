function [sol,tau] = generate_data(f,h,ic,steps)
VF = @(t,y,Z) f(t,y,Z);
LAG = @(t,y) h(t,y);
tau=zeros(steps,1);
% First step
opts = ddeset('Events',@(t,y,Z)crossing(t,y,Z,h,0));
sol = ddesd(@(t,y,Z)VF(t,y,Z),@(t,y)LAG(t,y),ic,[0,100],opts);
tau(1) = sol.x(end);
% Remaining steps
for k=2:steps
    opts = ddeset('Events',@(t,y,Z)crossing(t,y,Z,h,tau(k-1)));
    sol = ddesd(@(t,y,Z)VF(t,y,Z),@(t,y)LAG(t,y),sol,[sol.x(end),sol.x(end)+100],opts);
    tau(k) = sol.x(end);
end
end

function [value, isterminal, direction] = crossing(t,y,Z,h,cval)
value = h(t,y) - cval;
isterminal = 1;
direction = +1;
end