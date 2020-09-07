clc
clear classes
clear all
%close all
%%  Define parameters
dt = 0.001; %sampling time
T = 20; %Duration of the path
tSteps = 0:dt:T; %instant of time
nSteps = size(tSteps,2); %number of steps
tau_s = tSteps/T; %normalized time

%% Define the cubic trajectory
% A = [-0.1342,-0.1761,0.7726]; % for q=deg2rad([0 30 -30 45 60 15])';
% B = [2,-0.1761,0.7726];
% %B = [-0.1342,-0.1761,0.7726]; % for q=deg2rad([0 30 -30 45 60 15])';
% vin = [0,0,0];
% vf = [0,0,0]; 
% dp = A; cp = vin; bp = (-3*(A-B)-(vin+vf)*T)/power(T,2); ap = (2*(A-B)+(vin+vf)*T)/power(T,3);
% for i = 1:nSteps
%     X_d(i,:) = ap*power(tSteps(i),3) + bp*power(tSteps(i),2) + cp*tSteps(i) + dp;
%     dX_d(i,:) = 3*ap*power(tSteps(i),2) + 2*bp*tSteps(i) + cp;
% end

%% Define the circlular trajectory centered at (x0,y0), with radius 'R',
% initial phase 'phi', and offset from the Z-axis 'z_0'.
R= 0.35; %raduis
x0 = 0.3;y0 = 0.3;z_0 = 0.2;
phi = 0; %initial phase
dzt_0 = 0;
x = zeros(nSteps,1);y = zeros(nSteps,1);
dxs = zeros(nSteps,1);dys = zeros(nSteps,1);
dxt = zeros(nSteps,1);dyt = zeros(nSteps,1);
X_d = zeros(nSteps,3);dX_d = zeros(nSteps,3);ddX_d = zeros(nSteps,3);
for i = 1:nSteps
    s = 2*pi*(tau_s(i));
    ds = 2*pi/T;
    x(i) = x0 + R*cos(s+phi);
    y(i) = y0 + R*sin(s+phi);
    dxs(i) = -R*sin(s+phi);
    dys(i) = R*cos(s+phi);
    dxt(i) = dxs(i)*ds;
    dyt(i) = dys(i).*ds;
    X_d(i,:) = [x(i) y(i) z_0];
    dX_d(i,:) = [dxt(i) dyt(i) dzt_0];
end
%% Symbolic Jacobian
JCarsym1 = []; JdotCarsym1=[]; 
JCarsym2 = []; JdotCarsym2=[];
JCarsym3 = []; JdotCarsym3=[];
JCarsym4 = []; JdotCarsym4=[];
JCarsym5 = []; JdotCarsym5=[];
JCarsym6 = []; JdotCarsym6=[];
alpha = [pi/2,-pi/2,-pi/2,pi/2,pi/2,-pi/2,0];
a = [0,0,0,0,0,0,0];
d = [sym('d1'),0,sym('d3'),0,sym('d5'),0,sym('d7')];
qsym = [sym('q1'),sym('q2'),sym('q3'),sym('q4'),sym('q5'),sym('q6'),0];
dqsym = [sym('dq1');sym('dq2');sym('dq3');sym('dq4');sym('dq5');sym('dq6')]; %these are the joint velocities at evaluation point
[JCarsym1, JdotCarsym1] = computeJacobianSym(alpha(1:1),a(1:1),...
    d(1:1),qsym(1:1),dqsym(1:1));
[JCarsym2, JdotCarsym2] = computeJacobianSym(alpha(1:2),a(1:2),...
    d(1:2),qsym(1:2),dqsym(1:2));
[JCarsym3, JdotCarsym3] = computeJacobianSym(alpha(1:3),a(1:3),...
    d(1:3),qsym(1:3),dqsym(1:3));
[JCarsym4, JdotCarsym4] = computeJacobianSym(alpha(1:4),a(1:4),...
    d(1:4),qsym(1:4),dqsym(1:4));
[JCarsym5, JdotCarsym5] = computeJacobianSym(alpha(1:5),a(1:5),...
    d(1:5),qsym(1:5),dqsym(1:5));
[JCarsym6, JdotCarsym6] = computeJacobianSym(alpha(1:6),a(1:6),...
    d(1:6),qsym(1:6),dqsym(1:6));
%% Robot parameters for the KUKA-LWR
m = 3;
nJoints = 6;
d1=0.0;     % distance from base frame to second joint frame 
d3=0.4;     % distance from second joint frame to fourth joint frame
d5=0.39;    % distance from fourth joint frame to sixth joint frame 
d7=0.078;   % distance from sixth joint frame to EE frame; EE in the tip of KUKA without auxiliary addition

use_tau_f=false; % boolean telling whether tau_f should be used in the torque calculation or not

global robot
%Cartesian limits
X_max = 1*[1;1;1];
X_min = 1*[-1;-1;-1];
V_max = 1*[0.5;0.5;0.5];
V_min = 1*[-0.5;-0.5;-0.5];
A_max = 1*[0.5;0.5;0.5];
A_min = 1*[-0.5;-0.5;-0.5];
%% Control Loop
q = zeros(6,nSteps);
dq = zeros(6,nSteps);
ddq = zeros(6,nSteps);
Pee = zeros(4,nSteps);
dPee = zeros(3,nSteps);
T_stack_c = zeros(6,nSteps);
q(:,1)=deg2rad([0 -150 -60 -45 -15 0])';%for Pee=[0.5088,-0.2866,-0.5503] 6x1 vector with the current joint configuration (the 7th joint is assumed to be fixed q_7=0)  
%q(:,1)=deg2rad([0 30 -30 45 60 15])';
%q(:,1)=deg2rad([50 80 -90 -97 0 100])';
%q(:,1)=deg2rad([0 10 -60 45 30 0])';
dq(:,1)=zeros(6,1); % 6x1 vector with the current joint velocity (the 7th joint is assumed to be fixed q_7=0)
%%% no need to consider q7 if there is no EE desired orientation 
ddq(:,1) = zeros(6,1);

Kp = diag([80 80 80]); %position gain
Kv = 2*sqrt(Kp); %velocity gain
Kd = 10; %damping gain

tasks = [];
% if (0)
%     X_d = repmat(B,size(X_d,1),1);
%     dX_d = zeros(size(X_d));
% end
for i = 1:nSteps
    fprintf('-----------%d----------\n',i);
    robot=KUKA_LWR(d1,d3,d5,d7,use_tau_f,q(:,i),dq(:,i));
    [M,S,n,J,Jdot]=robot.computeMSnJJdot;
    [M,cc,g,J,Jdot]=robot.computeMccgJJdot;
    Pee(:,i) = robot.computeEnd_effectorPosition;
    dPee(:,i) = robot.computeEnd_effectorVelocity;
    ddPee(:,i) = Jdot*dq(:,i) + J*ddq(:,i);
    error_p(:,i) = (X_d(i,:)'-Pee(1:3,i));
    error_v(:,i) = (dX_d(i,:)'-dPee(:,i));
    f_pos_star = Kp*(X_d(i,:)'-Pee(1:3,i)) + Kv*(dX_d(i,:)'-dPee(:,i));
    f_damp_star = -Kd*dq(:,i);
    
    task1.f_star = f_pos_star;
    task1.J = J;
    %task1.Jsym = Jsym;
    task1.Jdot = Jdot;
    
    task2.f_star = f_damp_star;
    task2.J = eye(nJoints);
    %task2.Jsym = eye(n);
    task2.Jdot = zeros(nJoints);
    tasks = [task1;task2];
    %%%
    Parray(:,:,i) = robot.computeJointPositionArray;
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    q5 = q(5);
    q6 = q(6);
    
    dq1 = dq(1);
    dq2 = dq(2);
    dq3 = dq(3);
    dq4 = dq(4);
    dq5 = dq(5);
    dq6 = dq(6);
    symvalues = [d1 d3 d5 d7 q1 q2 q3 q4 q5 q6 dq1 dq2 dq3 dq4 dq5 dq6];
    JCar1 = double(subs(JCarsym1));
    JdotCar1 = double(subs(JdotCarsym1));
    JCar2 = double(subs(JCarsym2));
    JdotCar2 = double(subs(JdotCarsym2));
    JCar3 = double(subs(JCarsym3));
    JdotCar3 = double(subs(JdotCarsym3));
    JCar4 = double(subs(JCarsym4));
    JdotCar4 = double(subs(JdotCarsym4));
    JCar5 = double(subs(JCarsym5));
    JdotCar5 = double(subs(JdotCarsym5));
    JCar6 = double(subs(JCarsym6));
    JdotCar6 = double(subs(JdotCarsym6));
    dP1(:,i) = JCar1*dq(1:1,i);
    ddP1(:,i) = JdotCar1*dq(1:1,i) + JCar1*ddq(1:1,i);
    dP2(:,i) = JCar2*dq(1:2,i);
    ddP2(:,i) = JdotCar2*dq(1:2,i) + JCar2*ddq(1:2,i);
    dP3(:,i) = JCar3*dq(1:3,i);
    ddP3(:,i) = JdotCar3*dq(1:3,i) + JCar3*ddq(1:3,i);
    dP4(:,i) = JCar4*dq(1:4,i);
    ddP4(:,i) = JdotCar4*dq(1:4,i) + JCar4*ddq(1:4,i);
    dP5(:,i) = JCar5*dq(1:5,i);
    ddP5(:,i) = JdotCar5*dq(1:5,i) + JCar5*ddq(1:5,i);
    dP6(:,i) = JCar6*dq(1:6,i);
    ddP6(:,i) = JdotCar6*dq(1:6,i) + JCar6*ddq(1:6,i);
    %%%
    %T_stack_c(:,i) = getTorques(tasks, M, cc, g, dq(:,i), m, nJoints);%qsym, dqsym, symvalues
    %T_stack_c = J'*f_pos_star + n - S*dq(:,i);
    T_stack_c(:,i) = scs(M, cc, g, dq(:,i), J, Jdot, m, nJoints, Pee(1:3,i), dPee(:,i), dt, X_max, X_min, V_max, V_min, A_max, A_min, tasks, i);
    ddq(:,i+1) = M\(T_stack_c(:,i) - n);
    %%%%%
    %a = ddX_d(i,:)' + Kd*(dX_d(i,:)'-dPee(:,i)) + Kp*(X_d(i,:)'-Pee(1:3,i));
    %ddq(:,i+1) = pinv(J)*a - pinv(J)*Jdot*dq(:,i);
    dq(:,i+1) = dq(:,i) + ddq(:,i)*dt;
    q(:,i+1) = q(:,i) + dq(:,i)*dt + 0.5*ddq(:,i)*dt*dt;
end
P1 = squeeze(Parray(:,1,:));
P2 = squeeze(Parray(:,2,:));
P3 = squeeze(Parray(:,3,:));
P4 = squeeze(Parray(:,4,:));
P5 = squeeze(Parray(:,5,:));
P6 = squeeze(Parray(:,6,:));
size_Pee = size(Pee,2);
% figure;
% plot(tSteps(1:size_Pee),Pee(1,:));
% hold on
% plot(tSteps(1:size_Pee),X_d(1:size_Pee,1));
% xlabel('Time steps');
% ylabel('Traj x');
% legend('Output Traj x','Desired Traj x');
% 
% figure;
% plot(tSteps(1:size_Pee),Pee(2,:));
% hold on
% plot(tSteps(1:size_Pee),X_d(1:size_Pee,2));
% xlabel('Time steps');
% ylabel('Traj y');
% legend('Output Traj y','Desired Traj y');
% 
% figure;
% plot(tSteps(1:size_Pee),Pee(3,:));
% hold on
% plot(tSteps(1:size_Pee),X_d(1:size_Pee,3));
% xlabel('Time steps');
% ylabel('Traj z');
% legend('Output Traj z','Desired Traj z');

figure;
plot3(Pee(1,:)',Pee(2,:)',Pee(3,:)','linewidth',2);
hold on
plot3(X_d(:,1),X_d(:,2),X_d(:,3));
grid on
xlabel('Traj x');
ylabel('Traj y');
zlabel('Traj z');
legend('Output Traj','Desired Traj');
title('Cartesian Trajectory Tracking')

figure;
plot(tSteps(1:size_Pee),error_p(1,:),'linewidth',2);
hold on
plot(tSteps(1:size_Pee),error_p(2,:),'linewidth',2);
hold on
plot(tSteps(1:size_Pee),error_p(3,:),'linewidth',2);
grid on
xlabel('Time steps');
ylabel('error');
legend('in x','in y','in z');
title('Cartesian Error vs Time')

figure
plot(tSteps,T_stack_c(1,:),'linewidth',2);
hold on;
plot(tSteps,T_stack_c(2,:),'linewidth',2);
hold on;
plot(tSteps,T_stack_c(3,:),'linewidth',2);
hold on;
plot(tSteps,T_stack_c(4,:),'linewidth',2);
hold on;
plot(tSteps,T_stack_c(5,:),'linewidth',2);
hold on;
plot(tSteps,T_stack_c(6,:),'linewidth',2);
grid on
legend('j1','j2','j3','j4','j5','j6');
xlabel('T');
ylabel('Joint Torque');
title('Joint Torque vs Time')

figure
plot(tSteps,ddPee(1,:),'linewidth',2);
hold on;
plot(tSteps,ddPee(2,:),'linewidth',2);
hold on;
plot(tSteps,ddPee(3,:),'linewidth',2);
grid on
legend('x','y','z');
xlabel('T');
ylabel('Acceleration');
title('Cartesian Acceleration vs Time')

figure
plot(tSteps,dPee(1,:),'linewidth',2);
hold on;
plot(tSteps,dPee(2,:),'linewidth',2);
hold on;
plot(tSteps,dPee(3,:),'linewidth',2);
grid on
legend('x','y','z');
xlabel('T');
ylabel('Velocity');
title('Cartesian Velocity vs Time')

