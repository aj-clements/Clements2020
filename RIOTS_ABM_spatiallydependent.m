function RIOTS_ABM_spatiallydependent(kmr,kmb,kr,kd,Ttot,P)
%% Runs P simulations of ABM for riots model for spatially-dependent ICs
% Use this code to reproduce Figures 8, 9
% Use this code and RIOTS_contlim_spatiallydependent.m to reproduce Figures 12, 13

% kmr = movement rate for rioters (unbiased in direction)
% kmb = movement rate for bystanders (unbiased in direction)
% kr = recruitment rate
% kd = defection rate
% Ttot = Overall simulation time
% P = Number of identically-prepared simulations to run

%% Setup and Initial Conditions

X=200; % Number of x-ordinates
Y=20; % Number of y-ordinates

Pmr=kmr; % Moving prob/rate for rioters
Pmb=kmb; % Moving prob/rate for bystanders
Pr=kr; % Recruitment prob/rate
Pd=kd; % Defection rate/prob

CCR=zeros(Y,X); % Empty lattice for tracking rioters (0s will denote empty sites)
CCB=zeros(Y,X); % Empty lattice for tracking bystanders

% Initial conditions with a rioter strip and a bystander strip (comment out if not needed)

CCR(:,81:90)=1; % Initial conditions (no spatial-uniformity) Left end of avenue full of rioters (denoted by 1s)
CCB(:,111:120)=1; % Initial conditions (no spatial-uniformity) Right side of avenue full of bystanders (denoted by 1s)

% Initial conditions with a central rioter strip and two bystander strips (comment out if not needed)

% CCR(:,91:110)=1; % Initial conditions 2 (no spatial-uniformity) Central stripe of rioters
% CCB(:,61:80)=1; % Initial conditions 2 (no spatial-uniformity) Bystander crowd to left
% CCB(:,121:140)=1; % Initial conditions 2 (no spatial-uniformity) Bystander crowd to right

QRstart=0; % Counts the total number of rioters initially
QBstart=0; % Counts the total number of bystanders initially
for y=1:Y
    for x=1:X % Counts over every site on the lattice
        if CCR(y,x)==1 % There's a rioter in that site
            QRstart=QRstart+1; % Counts a rioter
        elseif CCB(y,x)==1 % There's a bystander in that site
            QBstart=QBstart+1; % Counts a bystander
        else
        end
    end
end % QRstart and QBstart now have the total numbers of rioters and bystanders

CRstart=CCR; % Initial profile of rioters
CBstart=CCB; % Initial profile of bystanders

%% Nearest neighbour index structure

XX=zeros(Y,X); % To contain row index i for each (i,j)
YY=zeros(Y,X); % To contain column index j for each (i,j)

for i=1:Y
    for j=1:X
            XX(i,j)=i;
            YY(i,j)=j;
    end
end

S=(X*Y+1)*ones(X*Y,4); % Set up matrix for neighbours, (XY+1) where less neighbours

for i=1:Y
    for j=1:X
        Z=(XX-XX(i,j)).^2+(YY-YY(i,j)).^2-1; % Metric for distance from (i,j)
        V=find(Z<=1e-1 & Z>-1); % Chooses sites neighbouring (i,j)
        S(i+(j-1)*Y,1:length(V))=V; % Lists neighbour indexes for (i,j) indexed by J, along row J
    end
end

S(S==0)=X*Y+1; % Buffer

%% Running Simulations

for p=1:P % Runs a simulation for each p from same initial conditions
    
    T=0; % Time starts at zero
    j=1; % Index to be used for steps
    
    tau=0; % Actually, to be a vector of cumulative timesteps
    QR=QRstart; % Actually, to be a vector of total numbers of rioters
    QB=QBstart; % Actually, to be a vector of total numbers of bystanders
    
    CR=CCR; % To manipulate the lattice
    CB=CCB;
    
    JJ=randi([1,X*Y],1,10000); % Vector (length 10000) of [1,XY] randm integers
    U1=rand(1,10000); % First vector (length 10000) of Unif(0,1)s
    U2=rand(1,10000); % Second vector (length 10000) of Unif(0,1)s
    U3=rand(1,10000); % Third vector (length 10000) of Unif(0,1)s
    U4=rand(1,10000); % Fourth vector (length 10000) of Unif(0,1)s
    y=1; % Arbitrary index
    QRend=QRstart; % To be updated after step
    QBend=QBstart; % To be updated after step
    
    %% Giant Gillespie Algorithm
    
    while T<Ttot % While simulation has time left to run
        prop=(Pmr+Pd)*QR(j)+(Pmb+Pr)*QB(j); % Computes propensity of the whole system
        tau(j+1)=tau(j)-log(U1(y))/prop; % Next entry of tau is time after timestep
        
        % First pick agent, then pick agent event later
        if U4(y)<(Pmr+Pd)*QR(j)/((Pmr+Pd)*QR(j)+(Pmb+Pr)*QB(j))
            % Agent chosen will be a rioter
            Agent=1;
            Ind = find(CR); % Finds lattice sites occupied by rioters
            J=Ind(randi(length(Ind))); % Pick a rioter
            v=S(J,:); % Contains indexes of neighbours of J
        else
            % Agent chosen will be a bystander
            Agent=2;
            Ind = find(CB); % Finds lattice sites occupied by bystanders
            J=Ind(randi(length(Ind))); % Pick a bystander
            v=S(J,:); % Contains indexes of neighbours of J
        end
        
        QR(j+1)=QR(j); % Creates a new element for defining later
        QB(j+1)=QB(j); % Creates a new element for defining later
        
        if Agent==1 && U2(y)<Pmr/(Pmr+Pd)  % If rioter chosen, in the event a rioter moves
            Jnext=v(randi([1,4])); % Chooses the direction to move with no bias
            if  Jnext<(X*Y+1) && CR(Jnext)==0 && CB(Jnext)==0 % Checks site is in matrix and not occupied
                CR(Jnext)=1; % Neighbouring site now contains rioter
                CR(J)=0; % Site originally selected now empties
            else % Movement aborted
            end
            
        elseif Agent==2 && U2(y)<Pmb/(Pmb+Pr) % If bystander chosen, in the event a bystander moves
            Jnext=v(randi([1,4])); % Chooses the direction to move with no bias
            if  Jnext<(X*Y+1) && CR(Jnext)==0 && CB(Jnext)==0 % Checks site is in matrix and not occupied
                CB(Jnext)=1; % Neighbouring site now contains bystander
                CB(J)=0; % Site originally selected now empties
            else % Movement aborted
            end
            
        elseif Agent==2 % If bystander chosen, in the event a bystander recruits
            if sum(CR(v(v~=X*Y+1)))>0 % If at least one rioter neighbours
                CB(J)=0; % Site originally selected loses bystander...
                CR(J)=1; % ... who turns into a rioter
                QR(j+1)=QR(j+1)+1; % Gain rioter
                QB(j+1)=QB(j+1)-1; % Lose bystander
            else % Bystander not near rioters and so not recruited
            end
            
        elseif Agent==1 % If rioter chosen, in the event a rioter defects
            if sum(CB(v(v~=X*Y+1)))>0 % If at least one bystander neighbours
                CR(J)=0; % Site originally selected loses rioter...
                CB(J)=1; % ... who turns into a bystander
                QR(j+1)=QR(j+1)-1; % Lose rioter
                QB(j+1)=QB(j+1)+1; % Gain bystander
            else % Rioter not near bystanders and so does not defect
            end
             
        end
        
        T=tau(j+1); % Updating time for the next iteration
        QRend=QR(j+1); % Updating number of rioters after step
        QBend=QB(j+1); % Updating number of bystanders after step
        
        y=y+1;
        if y==10001
            U1=rand(1,10000);
            U2=rand(1,10000);
            U3=rand(1,10000);
            U4=rand(1,10000);
            y=1;
        end
        
        j=j+1; % Increase step counter to proceed to next step
        
        % Agent snapshots
                 if mod(j,1000)==0
                    figure(401)
                    spy(CR,'r')
                    hold on
                    spy(CB,'b')
                    delete(findall(findall(gcf,'Type','axe'),'Type','text'))
                    hold off
                    pause(0.01)
                 end
         CRend=CR;
         CBend=CB;
    end
    
    %% Avergaing in each column
    
    CRmeanoverYs(p,:)=mean(CRend,1); % Averaging over the rows to get averaged rioter occupancy at each x
    CBmeanoverYs(p,:)=mean(CBend,1); % Averaging over the rows to get averaged bystander occupancy at each x
    
end

CRmeanoverYsPs=mean(CRmeanoverYs,1); % Average over simulations
CBmeanoverYsPs=mean(CBmeanoverYs,1); % Average over simulations
figure(402)
plot(0.5:1:199.5,CRmeanoverYsPs,'r','linewidth',4) % Plot averaged rioter occupancy
hold on
plot(0.5:1:199.5,CBmeanoverYsPs,'b','linewidth',4) % Plot averaged bystander occupancy
xlim([0,200])
ylim([0,1])
xlabel('x')
ylabel('C_{r}^{P}(x,t), C_{b}^{P}(x,t)')
legend({'C_{r}^{P}(x,t)' 'C_{b}^{P}(x,t)'})
hold off

end