function RIOTS_contlim_spatiallydependent(kmr,kmb,kr,kd,tend)
%% Spatially-dependent continuum limit considered in Section 4.2.2
% Use this code and RIOTS_ABM_spatiallydependent.m to reproduce Figures 12, 13

% kmr = rioter motility rate
% kmb = bystander motility rate
% kr = recruitment rate
% kd = defection rate
% tend = total run time

%% Growth rate Equations
% Writes the growth rates for the continuum limit system
% u(1) is (mean) rioter occupancy
% u(2) is (mean) bystander occupancy

Dr=kmr/4; % Effective rioter diffusivity
Db=kmb/4; % Effective bystander diffusivity

m=0;

    function [c,f,s]=pdefun(x,t,u,dudx)
        c=[1;1];
        f=[Dr*(1-u(2))*dudx(1)+Dr*u(1)*dudx(2);Db*(1-u(1))*dudx(2)+Db*u(2)*dudx(1)];
        s=[kr*u(2)*(1-(1-u(1))^4)-kd*u(1)*(1-(1-u(2))^4);-kr*u(2)*(1-(1-u(1))^4)+kd*u(1)*(1-(1-u(2))^4)];
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

figure(601)
plot(xval,r(end,:),'m--','linewidth',4)
hold on
plot(xval,b(end,:),'c--','linewidth',4)
xlim([0,200])
ylim([0,1])
xlabel('x')
ylabel('r(x,t), b(x,t)')
legend({'r(x,t)' 'b(x,t)'})
hold off
end
