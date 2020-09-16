function RIOTS_preliminary_ode(kr,kd,m,n,r0,b0,tend) 
%% Preliminary ODE system considered in Section 2.1
% Use this code to reproduce Figures 1, 2, 3
% Use this code and RIOTS_ABM_extended.m to reproduce Figures 15, 16

% kr = recruitment rate
% kd = defection rate
% m = 'number of rioters needed for a recruitment'
% n = 'number of bystanders needed for a defection'
% r0 = initial rioter density
% b0 = initial bystander density
% tend = total run time

%% System solution and plots

% Sets up time run
trun = [0 tend];
% 2x1 column vector u0 contains initial conditions
u0 = [r0; b0];
% Density total
K=r0+b0;
% Domain for phase portrait
r=linspace(0,K,1000);
% Growth rate function for phase portrait
v=r.*(K-r).*(kr*r.^(m-1)-kd*(K-r).^(n-1));
vmark=r0*(K-r0)*(kr*r0^(m-1)-kd*(K-r0)^(n-1));
% Solving the system
[t,u] = ode45(@f,trun,u0);
% Plot trajectories for 2-d system and phase portrait for 1-d system
figure(101)
plot(t,u(:,1),'r','linewidth',4) % Rioter trajectory
hold on
plot(t,u(:,2),'b','linewidth',4) % Bystander trajectory
hold off
xlabel('Time, t')
ylabel('r(t), b(t)')
legend({'r(t)' 'b(t)'})
figure(102)
plot(r,v,'k','linewidth',4) % Rioter density growth rate
hold on
plot(r,zeros(1,1000),'k--') % Emphasises the r-axis on the plot
plot(r0,vmark,'r*','MarkerSize',30) % Marker positioned at initial rioter density
hold off
xlabel('Rioter density, r')
ylabel('dr/dt')

%% ODE System
function dudt = f(t,u) 
% 2x1 column vector u contains variables r (as R below) and b (as B below)
    R=u(1);
    B=u(2);
    dR = (kr * R.^m .* B)-(kd * R .* B.^n) ;
    dB = -(kr * R.^m .* B)+(kd * R .* B.^n) ;
    
dudt = [dR; dB];
end

end