function RIOTS_preliminary_pde(Dr,Db,kr,kd,M,N,tend)
%% Preliminary PDE system considered in Section 2.2
% Use this code to reproduce Figures 4, 5

% Dr = effective diffusivity of rioters
% Db = effective diffusivity of bystanders
% kr = recruitment rate
% kd = defection rate
% M = 'number of rioters needed for a recruitment' (same as 'm' for ode system)
% N = 'number of bystanders needed for a defection' (same as 'n' for ode system)
% tend = total run time

%% Growth rate Equations
% Writes the growth rates for the continuum limit system
% u(1) is (mean) rioter occupancy
% u(2) is (mean) bystander occupancy

m=0;

    function [c,f,s]=pdefun(x,t,u,dudx)
        c=[1;1];
        f=[Dr*dudx(1);Db*dudx(2)];
        s=[kr*u(2)*u(1)^(M)-kd*u(1)*u(2)^(N);kd*u(1)*u(2)^(N)-kr*u(2)*u(1)^(M)];
    end

%% Initial Conditions (comment out conditions not needed)

% Initial conditions with a rioter strip and a bystander strip

    function u0=icfun(x)
        if x>=80 && x<=90
            u0=[1;0];
        elseif x>=110 && x<=120
            u0=[0;1];
        else
            u0=[0;0];
        end
    end

% Initial conditions with a central rioter strip and two bystander strips

%    function u0=icfun(x)
%        if x>=60 && x<=80||x>=120 && x<=140
%            u0=[0;1];
%        elseif x>=90 && x<=110
%            u0=[1;0];
%        else
%            u0=[0;0];
%        end
%    end

%% Boundary Conditions
% Neumann (zero flux) boundary conditions

    function [pl,ql,pr,qr]=bcfun(xl,ul,xr,ur,t)
        pl=[0;0];
        ql=[1;1];
        pr=[0;0];
        qr=[1;1];
    end

%% Spatial and time domains
% Defines 1d x-space 0<=x<=200 and time to run until 'tend'

xval=linspace(0,200,1000);
tval=linspace(0,tend,1000);

%% Solve the PDE
% Solves the continuum limit equations

sol=pdepe(m,@pdefun,@icfun,@bcfun,xval,tval);

r=sol(:,:,1); % Captures rioter distribution solution
b=sol(:,:,2); % Captures bystander distribution solution

%% Plot solutions in x-space

figure(201)
plot(xval,r(end,:),'r','linewidth',4) % Rioter density
hold on
plot(xval,b(end,:),'b','linewidth',4) % Bystander density
xlim([0,200])
ylim([0,1])
xlabel('x')
ylabel('r(x,t), b(x,t)')
legend({'r(x,t)' 'b(x,t)'})
hold off
end
