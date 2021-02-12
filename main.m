clc
clear
close all

%% parameters definition
global d1 d3 d5 d7 use_tau_f %KUKA-LWR parameters
global Kp Kv kd %control parameters
global tol Ts
%Robot parameters for the KUKA-LWR
d1=0.0;     % distance from base frame to second joint frame
d3=0.4;     % distance from second joint frame to fourth joint frame
d5=0.39;    % distance from fourth joint frame to sixth joint frame
d7=0.0675;   % distance from sixth joint frame to EE frame; EE in the tip of KUKA without auxiliary addition
use_tau_f=false; % boolean telling whether tau_f should be used in the torque calculation or not
%% Desired Trajectory
Ts = 0.005; %sampling time
T = 20; %Duration of the path (s)
t = 0:Ts:T; %instance of time
numSteps = size(t,2); %number of steps
tau_s = t/T; %normalized time
% Define the circlular trajectory in xz-plane centered at (x0,z0), with radius 'R',
% initial phase 'phi', and offset from the y-axis 'y0'.
R= 0.35; %raduis (m)
x0 = 0;y0 = -0.4;z0 = 0.4; % (m)
phi = 20; %initial phase
dyt_0 = 0;
z = zeros(numSteps,1);x = zeros(numSteps,1);
dzs = zeros(numSteps,1);dys = zeros(numSteps,1);
dzt = zeros(numSteps,1);dyt = zeros(numSteps,1);
Pd = zeros(3,numSteps,1);dPd = zeros(3,numSteps,1);
for k = 1:numSteps
    s = 2*pi*(tau_s(k));
    ds = 2*pi/T;
    z(k) = z0 + R*cos(s+phi);
    x(k) = x0 + R*sin(s+phi);
    dzs(k) = -R*sin(s+phi);
    dxs(k) = R*cos(s+phi);
    dzt(k) = dzs(k)*ds;
    dxt(k) = dxs(k).*ds;
    Pd(:,k) = [x(k);y0;z(k)];
    dPd(:,k) = [dxt(k);dyt_0;dzt(k)];
end
%% Set simulation parameters
Tstop = T; %Duration of simulation (s)
numSimSteps = size(0:Ts:Tstop,2);
% Initialization
q = zeros(6,numSimSteps);dq = zeros(6,numSimSteps);
ddq = zeros(6,numSimSteps);tau_stack_SJS = zeros(6,numSimSteps);
Pee = zeros(4,numSimSteps);Vee = zeros(3,numSimSteps);
Aee = zeros(3,numSimSteps);A_elbow = zeros(3,numSimSteps);
P_elbow = zeros(4,numSimSteps);V_elbow = zeros(3,numSimSteps);
F = zeros(3,numSimSteps);tau = zeros(6,numSimSteps);

IK = 0; %Flag for finding initial config. via Inverse Kinematics
%(0 = defualt initial config., 1 = find a possible new initial config.)

%% Assigning the constraints
% Joints limits
QMax = deg2rad([165; 115; 165; 115; 165; 115]);
QMin = -QMax;
VJntMax = deg2rad([100; 15; 15; 130; 130; 180]);
VJntMin = -VJntMax;
AJntMax = deg2rad([300; 300; 300; 300; 300; 300]);
AJntMin = -AJntMax;

% Cartesian limits
axs = 2; %Index of the limiting axis (y-axis)
XMax = [100;100;100];
XMin = [-100;-0.3;-100];
VXMax = [100;0.03;100];
VXMin = -[100;100;100];
AXMax = [100;0.2;100];
AXMin = -AXMax;

%% Assigning the control parameters
Kp = 8*diag([100,100,100]);
Kv = 2*sqrt(Kp);
kd = norm(Kp)/80;

% Tasks structure
task = struct('J',{},'Jdot',{},'F_ast',{});

tol = 1e-4; %Tolerance

%% Inverse kinematics based on gradient method
a = QMin;b = QMax;
Q(:,1) = (b-a).*rand(6,1) + a; %generate random inital config.
robot=KUKA_LWR(d1,d3,d5,d7,use_tau_f,Q(:,1),dq(:,1));
PEE = robot.computeEnd_effectorPosition;
i = 1;alpha = 0.3;
if (IK)
    while (norm(PEE(1:3)-Pd(:,1))>0.001 || any(Q(:,i)>QMax) || any(Q(:,i)<QMin) || any(PElbow(1:3)>XMax) || any(PElbow(1:3)<XMin) )
        robot=KUKA_LWR(d1,d3,d5,d7,use_tau_f,Q(:,i),dq(:,1));
        PElbow = robot.computeElbowPosition;
        PEE = robot.computeEnd_effectorPosition
        J = robot.computeJacobian;
        Q(:,i+1) = Q(:,i) + alpha*J'*(Pd(:,1)-PEE(1:3));
        i = i + 1;
        if i>10000 %search again
            Q(:,1) = (b-a).*rand(6,1) + a;
            i=1;
            robot=KUKA_LWR(d1,d3,d5,d7,use_tau_f,Q(:,i),dq(:,1));
            PEE = robot.computeEnd_effectorPosition;
        end
    end
end
%% Set the initial configuration
if (~IK)
    q(:,1) = [-0.625659084174681;-1.00333763134685;-2.19333615638075;0.598679120275943;1.71848417978296;1.69603819481226];
else
    q(:,1) = Q(:,i);
end

robot = KUKA_LWR(d1,d3,d5,d7,use_tau_f,q(:,1),dq(:,1));
Pee(:,1) = robot.computeEnd_effectorPosition; %initial position of the end-effector
Vee(:,1) = robot.computeEnd_effectorVelocity; %initial velocity of the end-effector
P_elbow(:,1) = robot.computeElbowPosition; %initial position of the elbow
V_elbow(:,1) = robot.computeElbowVelocity; %initial velocity of the elbow

%% The Main Loop
for k = 1:numSimSteps
    robot = KUKA_LWR(d1,d3,d5,d7,use_tau_f,q(:,k),dq(:,k));
    [M,S,n,J,Jdot] = robot.computeMSnJJdot; %dynamic model
    J_elbow = robot.computeElbowJacobian; %compute jacobian of the elbow
    dJ_elbow = robot.computeElbowJacobiandot; %compute derivative of elbow jacobian
    tasks = genTasks(J,Jdot,Pd(:,k),Pee(1:3,k),dPd(:,k),Vee(:,k),dq(:,k)); %generate tasks
    [tasksOutputSCS] = scs(M, n, J_elbow, dJ_elbow,dq(:,k), P_elbow(1:3,k), V_elbow(:,k), XMin, XMax, VXMin, VXMax, AXMin, AXMax, tasks);
    [tau_stack_SJS(:,k)] = sjs(M, n, J, q(:,k),dq(:,k), QMin, QMax, VJntMin, VJntMax, AJntMin, AJntMax, tasksOutputSCS);
    ddq(:,k) = M^(-1)*(tau_stack_SJS(:,k) - n); %compute the acceleration of the computed torque
    dq(:,k+1) = dq(:,k) + ddq(:,k)*Ts; %integrate to obtain joint velocity
    q(:,k+1) = q(:,k) + dq(:,k)*Ts + 0.5*ddq(:,k)*Ts.^2; %integrate to obtain joint position
    Aee(:,k) = J*ddq(:,k) + Jdot*dq(:,k); %compute the acceleration of the end-effector
    A_elbow(:,k) = J_elbow*ddq(:,k) + dJ_elbow*dq(:,k); %compute the acceleration of the elbow
    robot = KUKA_LWR(d1,d3,d5,d7,use_tau_f,q(:,k+1),dq(:,k+1)); 
    Pee(:,k+1) = robot.computeEnd_effectorPosition; %compute end-effector position
    Vee(:,k+1) = robot.computeEnd_effectorVelocity; %compute end-effector velocity
    P_elbow(:,k+1) = robot.computeElbowPosition; %compute elbow position
    V_elbow(:,k+1) = robot.computeElbowVelocity; %compute elbow velocity
end
%% Plot results
results_plot;