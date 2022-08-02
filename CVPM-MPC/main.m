%% Initialize
addpath functions
clear all

%% Init
Ts = 0.1;
T_sim = 10;
steps = floor(T_sim/Ts);

time = tic;     % Time evaluation of method
%Horizon
Nmpc  = 10;     % Horizon for MPC
Ncvpm = 7;     % Horizon for CVPM
% System 
A = [0 0 1 0;
     0 0 0 1;
     0 0 0 0;
     0 0 0 0];
B =    [0 0;     
        0 0;
        1 0;
        0 1];
 
A=eye(4)+Ts*A;

     
B =  Ts*B ;     

G = eye(4);

% Dimensions
nx = size(A,1);         % Dimension of state vector
nu = size(B,2);         % Dimension of input vector
nw = size(G,2);         % Dimension of distrubance vector

% System uncertainty
Sigmaw = 0.2*eye(nw);

% initial state
x0 = [200; 150;0;0];


Q  = [1 0 0 0;
      0 1 0 0;
      0 0 0 0;
      0 0 0 0];
  
R  = 0.0001*eye(nu);
%% Trajectories
X_ref  = sys_trajectory(x0,10,Ts,steps,Nmpc);    % reference trajectorie for the controlled system

%% Constraints
Px = 360*Polyhedron.unitBox(nx);         % State Constraints  E*x_k<e   k in [0;N-1]
Pu = 10* Polyhedron.unitBox(nu);         % Input Constraints  F*u_k<f   k in [0;N-1]
Pw = 1*Polyhedron.unitBox(nw);      % Distrubance Constraint
%% CVPM Constraint

[cspace,K]=build([4,4])
P=K;
T=100*Polyhedron.unitBox(2);
Pxp=P*T;
Pxp=Pxp.intersect(Px);
Pxp12=Pxp.projection(1:2);
% occur at timestep
kp=50;

%% Variable declaration for logging
X_log = [x0 zeros(nx,steps-1)];
U_log = zeros(nu,steps);
TAU_log = zeros(nu,steps);
[A_,B_] = liftedModel(A,B,Nmpc);


%% CVPM
time_CVPM = myTime('CVPM');
CVPM = myCVPM('A',A,'B',B,'G',G,'Nmpc',Nmpc,'Ncvpm',Ncvpm,'Q',Q,'R',R,'Pu',Pu,'Pw',Pw,'Sigma',Sigmaw,'Stability','Riccati');
% set state set as propabilistic constraint
% for CVPM.setConstraint(Px, 'stable') Pt= ? not found;
CVPM.setConstraint(Px);
CVPM.data_X_ref=X_ref;

%% System
Rob=createRobot();
sys= mySystem(Rob, G,  Px, Pu,  Pw,  Sigmaw,Ts);
%% Simulation
time_CVPM.offline
x=x0;
for i = 1:steps
    if i == kp
        CVPM.setConstraint(Pxp.intersect(Px));
        CVPM.Pxp12=Pxp12;
        sys.Px=CVPM.Pxp;
    end
        u = CVPM.do(x,X_ref(:,i:i+Nmpc-1));
    % Log
    X_log(:,i) = x;
    U_log(:,i) = u;
    
    % System
    sys.option='uni';
    x=sys.do(x,u,i);
    time_CVPM.step
    %plot
    figure(7)
    hold on
    plot(X_log(1,1:i),X_log(2,1:i),'o',X_ref(1,:),X_ref(2,:),'b');
    Pxp12.plot('color','blue','alpha',0.3);
    set(gcf,'color','w');
    hold off
    %axis([-1,5,-1,5]);
    xlabel('x_1/x_{ref,1}')
    ylabel('x_2/x_{ref,2}')
    title(['Simulation, time: ',num2str(Ts*i),' s'])  ;
end
time_CVPM.online
X_log=X_log*pi/180;
figure(2);
hold on
plot(Rob, X_log(1:2,:)','tile1color',[1 1 1],'notiles','top');
Ap=[1,0;
    -1,0;
    0,1;
    0,-1];
bp=[7;1;12;-4];
P = Polyhedron('A',Ap,'b',bp);
P.plot
hold off
% plot
CVPM.subplot(5);
CVPM.plot(6);



