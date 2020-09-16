function RIOTS_example_ode(r0,b0,tend) 
%% First ODE Ssystem considered in Section 5.4
% Use this code and RIOTS_ABM_extended.m to reproduce Figure 14

% r0 = initial rioter density
% b0 = initial bystander density
% tend = total run time

% Sets up time run
trun = [0 tend];
% 2x1 column vector u0 contains initial conditions
u0 = [r0; b0];

% Solving the system
[t,u] = ode45(@f,trun,u0);
% Plot trajectories for 2-d system
figure(801)
plot(t,u(:,1),'m','linewidth',4)
hold on
plot(t,u(:,2),'c','linewidth',4)
hold off
xlabel('Time, t')
ylabel('r(t), b(t)')
legend({'r(t)' 'b(t)'})


%% ODE System
function dudt = f(t,u) 
% 2x1 column vector u contains variables r (as R below) and b (as B below)
    R=u(1);
    B=u(2);
    dR = -(R.^2 .* B)+(R .* B.^2) ;
    dB = (R.^2 .* B)-(R .* B.^2) ;
    
dudt = [dR; dB];
end

end