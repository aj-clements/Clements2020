function RIOTS_contlim_spatiallyuniform(kr,kd,r0,b0,tend)
%% Spatially-uniform continuum limit considered in Section 4.2.1
% Use this code and RIOTS_ABM_spatiallyuniform.m to reproduce Figure 11

% kr = recruitment rate
% kd = defection rate
% r0 = initial rioter density
% b0 = initial bystander density
% tend = total run time

%% Growth rate equations

function dudt=f(t,u)
    dr=kr*u(2)*(1-(1-u(1))^4)-kd*u(1)*(1-(1-u(2))^4);
    db=-kr*u(2)*(1-(1-u(1))^4)+kd*u(1)*(1-(1-u(2))^4);
dudt=[dr;db];
end

%% Initial Conditions

u0=[r0;b0];

%% Simulation time run

trun=[0 tend];

%% Solving the ODE system

[t,u]=ode45(@f,trun,u0);

%% Plotting trajectories

figure(501)
plot(t,u(:,1),'m--','linewidth',4)
hold on
plot(t,u(:,2),'c--','linewidth',4)
hold off
xlabel('Time, t')
ylabel('r(t), b(t)')
legend({'r(t)' 'b(t)'})

end