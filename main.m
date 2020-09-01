clc
clear classes
clear all
close all
%%  Desired Circular Trajectory
dt = 0.001; %sampling time
T = 10; %Duration of the path
tSteps = 0:dt:T; %instant of time
nSteps = size(tSteps,2); %number of steps
tau_s = tSteps/T; %normalized time
% Define the circlular trajectory centered at (x0,y0), with radius 'R',
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
%% Robot parameters for the KUKA-LWR
m = 3;
nJoints = 6;
d1=0.0;     % distance from base frame to second joint frame 
d3=0.4;     % distance from second joint frame to fourth joint frame
d5=0.39;    % distance from fourth joint frame to sixth joint frame 
d7=0.078;   % distance from sixth joint frame to EE frame; EE in the tip of KUKA without auxiliary addition

use_tau_f=false; % boolean telling whether tau_f should be used in the torque calculation or not

global robot
%% Control Loop
q(:,1)=deg2rad([0 -150 -60 -45 -15 0])';%for Pee=[0.0360,-0.2866,0.7486] 6x1 vector with the current joint configuration (the 7th joint is assumed to be fixed q_7=0)  
%q(:,1)=deg2rad([50 80 -90 -97 0 100])';
dq(:,1)=zeros(6,1); % 6x1 vector with the current joint velocity (the 7th joint is assumed to be fixed q_7=0)
%%% no need to consider q7 if there is no EE desired orientation 
ddq(:,1) = zeros(6,1);

Kp = 20; %position gain
Kv = 10; %velocity gain
Kd = 10; %dampling gain

tasks = [];

for i = 1:nSteps
    %fprintf('-----------%d----------',i);
    robot=KUKA_LWR(d1,d3,d5,d7,use_tau_f,q(:,i),dq(:,i));
    [M,S,n,J,Jdot]=robot.computeMSnJJdot;
    [M,cc,g,J,Jdot]=robot.computeMccgJJdot;
    Pee(:,i) = robot.computeEnd_effectorPosition;
    dPee(:,i) = robot.computeEnd_effectorVelocity;
    error_p(:,i) = Kp*(X_d(i,:)'-Pee(1:3,i));
    error_v(:,i) = Kd*(dX_d(i,:)'-dPee(:,i));
    %%%
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
    [T_stack] = getTorques(tasks, M, cc, g, dq(:,i), m, nJoints);%qsym, dqsym, symvalues
    ddq(:,i+1) = inv(M)*(T_stack - n);
    %%%%%
    a = ddX_d(i,:)' + Kd*(dX_d(i,:)'-dPee(:,i)) + Kp*(X_d(i,:)'-Pee(1:3,i));
    ddq(:,i+1) = pinv(J)*a - pinv(J)*Jdot*dq(:,i);
    dq(:,i+1) = dq(:,i) + ddq(:,i)*dt;
    q(:,i+1) = q(:,i) + dq(:,i)*dt + 0.5*ddq(:,i)*dt*dt;
end

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
plot3(X_d(:,1),X_d(:,2),X_d(:,3),'linewidth',2);
xlabel('Traj x');
ylabel('Traj y');
zlabel('Traj z');
legend('Output Traj','Desired Traj');

figure;
plot(tSteps(1:size_Pee),error_p(1,:)/Kp);
hold on
plot(tSteps(1:size_Pee),error_p(2,:)/Kp);
hold on
plot(tSteps(1:size_Pee),error_p(3,:)/Kp);
xlabel('Time steps');
ylabel('error');
legend('in x','in y','in z');