clc
clear classes
clear all
%close all
%%  Define parameters
dt = 0.001; %sampling time
T = 10; %Duration of the path
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


%% Robot parameters for the KUKA-LWR
m = 3;
nJoints = 6;
d1=0.0;     % distance from base frame to second joint frame 
d3=0.4;     % distance from second joint frame to fourth joint frame
d5=0.39;    % distance from fourth joint frame to sixth joint frame 
d7=0.078;   % distance from sixth joint frame to EE frame; EE in the tip of KUKA without auxiliary addition

use_tau_f=false; % boolean telling whether tau_f should be used in the torque calculation or not

global robot
%% Limits
%Joint level limits
Q_max = deg2rad([165; 115; 165; 115; 165; 115]);
Q_min = -Q_max;
dQ_max = deg2rad([100; 10; 15; 130; 130; 180]);
dQ_min = -dQ_max;
ddQ_max = deg2rad([300; 300; 300; 300; 300; 300]);
ddQ_min = -ddQ_max;
% Q_max = pi*[0.7;0.25;0.8;0.7;0.9;0.85];
% Q_min = -Q_max;
% dQ_max = 5.5*ones(6,1);
% dQ_min = -dQ_max;
% ddQ_max = 5.5*ones(6,1);
% ddQ_min = -ddQ_max;
%Cartesian limits
X_max = [10;10;0.5];
X_min = [-10;-10;-0.5];
V_max = [10;10;0.2];
V_min = [-10;-10;-0.2];
A_max = [10;10;0.2];
A_min = [-10;-10;-0.2];
%% Symbolic Jacobian
nCar = 4;
JCarsym = []; JdotCarsym=[]; flagCar = 0;
if(nCar>0 && nCar<7)
    flagCar = 1;
    alpha = [pi/2,-pi/2,-pi/2,pi/2,pi/2,-pi/2,0];
    a = [0,0,0,0,0,0,0];
    d = [d1,0,d3,0,d5,0,d7];
    qsym = [sym('q1'),sym('q2'),sym('q3'),sym('q4'),sym('q5'),sym('q6'),0];
    dqsym = [sym('dq1');sym('dq2');sym('dq3');sym('dq4');sym('dq5');sym('dq6');0]; %these are the joint velocities at evaluation point
    [JCarsym, JdotCarsym] = computeJacobianSym(alpha(1:nCar),a(1:nCar),...
                        d(1:nCar),qsym(1:nCar),dqsym(1:nCar));
end
%% Control Loop
q = zeros(6,nSteps);
dq = zeros(6,nSteps);
ddq = zeros(6,nSteps);
Pee = zeros(4,nSteps);
dPee = zeros(3,nSteps);
T_stack_c = zeros(6,nSteps);
%q(:,1)=deg2rad([0 -150 -60 -45 -15 0])';%for Pee=[0.5088,-0.2866,-0.5503] 6x1 vector with the current joint configuration (the 7th joint is assumed to be fixed q_7=0)  
%q(:,1) = [0.1;-0.9;-1;0;-1;1];
%q(:,1) = [2.1;0.6;-2.4;-0.3;2.8;0.9];
%q(:,1)=deg2rad([0 30 -30 45 60 15])';
q(:,1)=deg2rad([50 80 -90 -97 0 100])';
%q(:,1)=deg2rad([0 10 -60 45 30 0])';
dq(:,1)=zeros(6,1); % 6x1 vector with the current joint velocity (the 7th joint is assumed to be fixed q_7=0)
%%% no need to consider q7 if there is no EE desired orientation 
ddq(:,1) = zeros(6,1);

% Kp = diag([800 800 800]); %position gain
% Kv = 2*sqrt(Kp); %velocity gain
% Kd = 10; %damping gain
Kp = diag([100,100,100]);
Kv = 2*sqrt(Kp);
Kd = norm(Kp)/80;
tasks = [];
% if (0)
%     X_d = repmat(B,size(X_d,1),1);
%     dX_d = zeros(size(X_d));
% end
robot=KUKA_LWR(d1,d3,d5,d7,use_tau_f,q(:,1),dq(:,1));
[M,S,n,J,Jdot]=robot.computeMSnJJdot;
Pee(:,1) = robot.computeEnd_effectorPosition;
%% Define the circlular trajectory centered at (x0,y0), with radius 'R',
% initial phase 'phi', and offset from the Z-axis 'z_0'.
% R=0.35; %raduis
% %x0 = Pee(1,1)-2*R;y0 = Pee(2,1);z_0 = Pee(3,1);
% x0 = 0.3;y0 = 0.3;z0 = 0.2;
% phi = 0; %initial phase
% dxt_0 = 0;
% z = zeros(nSteps,1);y = zeros(nSteps,1);
% dzs = zeros(nSteps,1);dys = zeros(nSteps,1);
% dzt = zeros(nSteps,1);dyt = zeros(nSteps,1);
% X_d = zeros(nSteps,3);dX_d = zeros(nSteps,3);ddX_d = zeros(nSteps,3);
% for i = 1:nSteps
%     s = 2*pi*(tau_s(i));
%     ds = 2*pi/T;
%     z(i) = z0 + R*cos(s+phi);
%     y(i) = y0 + R*sin(s+phi);
%     dzs(i) = -R*sin(s+phi);
%     dys(i) = R*cos(s+phi);
%     dzt(i) = dzs(i)*ds;
%     dyt(i) = dys(i).*ds;
%     X_d(i,:) = [x0 y(i) z(i)];
%     dX_d(i,:) = [dxt_0 dyt(i) dzt(i)];
% end
R= 0.35; %raduis
x0 = 0;y0 = -0.4;z0 = 0.4;
phi = deg2rad(20); %initial phase
dyt_0 = 0;
z = zeros(nSteps,1);x = zeros(nSteps,1);
dzs = zeros(nSteps,1);dxs = zeros(nSteps,1);
dzt = zeros(nSteps,1);dxt = zeros(nSteps,1);
X_d = zeros(nSteps,3);dX_d = zeros(nSteps,3);ddX_d = zeros(nSteps,3);
for k = 1:nSteps
    s = 2*pi*(tau_s(k));
    ds = 2*pi/T;
    z(k) = z0 + R*cos(s+phi);
    x(k) = x0 + R*sin(s+phi);
    dzs(k) = -R*sin(s+phi);
    dxs(k) = R*cos(s+phi);
    dzt(k) = dzs(k)*ds;
    dxt(k) = dxs(k).*ds;
    X_d(k,:) = [x(k) y0 z(k)];
    dX_d(k,:) = [dxt(k) dyt_0 dzt(k)];
end
%%
for i = 1:nSteps
    fprintf('\n-----------C%d----------\n',i);
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
    if(flagCar)
        PCar(:,i) = Parray(:,nCar,i);
        q1 = q(1,i);
        q2 = q(2,i);
        q3 = q(3,i);
        q4 = q(4,i);
        q5 = q(5,i);
        q6 = q(6,i);
        
        dq1 = dq(1,i);
        dq2 = dq(2,i);
        dq3 = dq(3,i);
        dq4 = dq(4,i);
        dq5 = dq(5,i);
        dq6 = dq(6,i);
        symvalues = [q1 q2 q3 q4 q5 q6 dq1 dq2 dq3 dq4 dq5 dq6];
        JCar = double(subs(JCarsym));
        JdotCar = double(subs(JdotCarsym));
        dPCar(:,i) = JCar*dq(1:nCar,i);
        ddPCar(:,i) = JdotCar*dq(1:nCar,i) + JCar*ddq(1:nCar,i);
        JCar = [JCar,zeros(3,(6-nCar))];
        JdotCar = [JdotCar,zeros(3,(6-nCar))];
    end

%     JCar = [J(:,1:nCar),zeros(3,(6-nCar))];
%     JdotCar = [Jdot(:,1:nCar),zeros(3,(6-nCar))];
%     PCar(:,i) = Parray(:,nCar,i);
%     dPCar(:,i) = JCar*dq(:,i);
%     ddPCar(:,i) = JCar*ddq(:,i)+JdotCar*dq(:,i);
    %%%
    %T_stack_c(:,i) = getTorques(tasks, M, cc, g, dq(:,i), m, nJoints);%qsym, dqsym, symvalues
    %T_stack_c = J'*f_pos_star + n - S*dq(:,i);
    [T_stack_cc(:,i),tasks_scs] = scs(M, cc, g, dq(:,i), J, Jdot, m, nJoints, Pee(1:3,i), dPee(:,i),...
                                             dt, X_max, X_min, V_max, V_min, A_max, A_min, tasks, i, JCar,...
                                             JdotCar, PCar(1:3,i), dPCar(:,i), nCar);
    fprintf('\n-----------J%d----------\n',i);
    [T_stack_c(:,i),tasks_sjs] = sjs(M, cc, g, (dq(:,i)), J, Jdot, m, nJoints, dt, (q(:,i)),...
                                 (Q_max), (Q_min), (dQ_max), (dQ_min), (ddQ_max), (ddQ_min), tasks_scs, i);

    ddq(:,i+1) = M\(T_stack_c(:,i) - n);
    %%%%%
    %a = ddX_d(i,:)' + Kd*(dX_d(i,:)'-dPee(:,i)) + Kp*(X_d(i,:)'-Pee(1:3,i));
    %ddq(:,i+1) = pinv(J)*a - pinv(J)*Jdot*dq(:,i);
    dq(:,i+1) = dq(:,i) + ddq(:,i+1)*dt;
    q(:,i+1) = q(:,i) + dq(:,i+1)*dt + 0.5*ddq(:,i+1)*dt*dt;
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
hold on
plot3(Pee(1,1)',Pee(2,1)',Pee(3,1)','-r*');
hold on
plot3(X_d(1,1),X_d(1,2),X_d(1,3),'-bo');
grid on
xlabel('Traj x');
ylabel('Traj y');
zlabel('Traj z');
legend('Output Traj','Desired Traj','Initial','Desired first');
title('Cartesian Trajectory Tracking')

figure
hold on
plot(tSteps,T_stack_c(:,:),'linewidth',2);
grid on
legend('j1','j2','j3','j4','j5','j6');
xlabel('T');
ylabel('Joint Torque');
title('Joint Torque vs Time')

%% End Effector Figures
figure;
plot(tSteps(1:size_Pee),error_p(:,:),'linewidth',2);
grid on
xlabel('Time steps');
ylabel('error');
legend('in x','in y','in z');
title('End-effector Error vs Time')

figure
plot(tSteps,dPee(:,:),'linewidth',2);
grid on
legend('x','y','z');
xlabel('T');
ylabel('Velocity');
title('End-effector Velocity vs Time')

figure
hold on
plot(tSteps,ddPee(:,:),'linewidth',2);
grid on
legend('x','y','z');
xlabel('T');
ylabel('Acceleration');
title('End-effector Acceleration vs Time')

%% Elbow Figures
size_t = size(tSteps,2);
figure
%plot(tSteps,dPCar(:,:),'linewidth',2);
plot(tSteps,bsxfun(@rdivide,PCar(1:3,1:size_t),X_max),'linewidth',2)
grid on
legend('x\_elbow','y\_elbow','z\_elbow');
xlabel('T');
ylabel('Position/X\_max');
title('Elbow Position/X\_max vs Time')

figure
%plot(tSteps,dPCar(:,:),'linewidth',2);
plot(tSteps,bsxfun(@rdivide,dPCar(:,1:size_t),V_max),'linewidth',2)
grid on
legend('x\_elbow','y\_elbow','z\_elbow');
xlabel('T');
ylabel('Velocity/V\_max');
title('Elbow Velocity/V\_max vs Time')

figure
hold on
%plot(tSteps,ddPCar(:,:),'linewidth',2);
plot(tSteps,bsxfun(@rdivide,ddPCar(:,1:size_t),A_max),'linewidth',2)
grid on
legend('x\_elbow','y\_elbow','z\_elbow');
xlabel('T');
ylabel('Elbow Acceleration/A\_max');
title('Elbow Acceleration/A\_max vs Time')

%% Joint Figures
size_t = size(tSteps,2);
figure
%plot(tSteps,q(:,1:size_t),'linewidth',2);
plot(tSteps,bsxfun(@rdivide,q(:,1:size_t),Q_max),'linewidth',2)
grid on
legend('j1','j2','j3','j4','j5','j6');
xlabel('T');
ylabel('Joint Position/Q\_max');
title('Joint Position/Q\_max vs Time')

figure
%plot(tSteps,dq(:,1:size_t),'linewidth',2);
plot(tSteps,bsxfun(@rdivide,dq(:,1:size_t),dQ_max),'linewidth',2);
grid on
legend('j1','j2','j3','j4','j5','j6');
xlabel('T');
ylabel('Joint Velocity/V\_max');
title('Joint Velocity/V\_max vs Time')

figure
%plot(tSteps,ddq(:,1:size_t),'linewidth',2);
plot(tSteps,bsxfun(@rdivide,ddq(:,1:size_t),ddQ_max),'linewidth',2);
grid on
legend('j1','j2','j3','j4','j5','j6');
xlabel('T');
ylabel('Joint Acceleration/A\_max');
title('Joint Acceleration/A\_max vs Time')

