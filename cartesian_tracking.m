clc
clear classes
clear all
close all
%% Desired Trajectory
%A = [5,1.5,2];
%A = [-0.3404,-0.0156,0.4884];
%A = [0.0360,-0.2866,0.7486]; %[0 10 -60 45 30 0]
A = [0.5088,0.2866,-0.5503]; %q=deg2rad([0 -150 -60 -45 -15 0])';
B = [-0.1440,-0.6546,0.5040]; %q=deg2rad([60 60 -60 30 15 0])';  
%B = [3,-0.0156,0.4884];
%B = [2,-2.5,2];
%B = [-0.1342,-0.1761,0.7726]; % for q=deg2rad([0 30 -30 45 60 15])';
%V = 2.5;
T = 5;
vin = [0,0,0];
vf = [0,0,0]; 
dt = 0.001; %sampling time

tSteps = 0:dt:T;
nSteps = size(tSteps,2);
tau = tSteps/T; %tau is between 0 and 1
X_d = zeros(nSteps,3);
dX_d = zeros(nSteps,3);
ddX_d = zeros(nSteps,3);
% for i = 1:nSteps
%     s_tau = -2*tau(i)*tau(i)*tau(i) + 3*tau(i)*tau(i);
%     ds_tau = -6*tau(i)*tau(i)/T + 6*tau(i)/T;
%     X_d(i,:) = A + (B-A)*s_tau;
%     dX_d(i,:) = (B-A)*ds_tau;
% end
dp = A; cp = vin; bp = (-3*(A-B)-(vin+vf)*T)/power(T,2); ap = (2*(A-B)+(vin+vf)*T)/power(T,3);
for i = 1:nSteps
    X_d(i,:) = ap*power(tSteps(i),3) + bp*power(tSteps(i),2) + cp*tSteps(i) + dp;
    dX_d(i,:) = 3*ap*power(tSteps(i),2) + 2*bp*tSteps(i) + cp;
    ddX_d(i,:) = 6*ap*tSteps(i) + 2*bp;
end
%% Robot parameters for the KUKA-LWR
d1=0.0;     % distance from base frame to second joint frame 
d3=0.4;     % distance from second joint frame to fourth joint frame
d5=0.39;    % distance from fourth joint frame to sixth joint frame 
d7=0.078;   % distance from sixth joint frame to EE frame; EE in the tip of KUKA without auxiliary addition

use_tau_f=false; % boolean telling whether tau_f should be used in the torque calculation or not

global robot
%% Control Loop
q(:,1)=deg2rad([0 -150 -60 -45 -15 0])';%for Pee=[0.0360,-0.2866,0.7486] 6x1 vector with the current joint configuration (the 7th joint is assumed to be fixed q_7=0)  
dq(:,1)=zeros(6,1); % 6x1 vector with the current joint velocity (the 7th joint is assumed to be fixed q_7=0)
%%% no need to consider q7 if there is no EE desired orientation 
ddq(:,1) = zeros(6,1);
Kp = diag(1*[75,20,20]);
%Ki = 1;
Kd = diag(1*[25,1,1]);

for i = 1:nSteps
    %fprintf('-----------%d----------',i);
    robot=KUKA_LWR(d1,d3,d5,d7,use_tau_f,q(:,i),dq(:,i));
    [M,S,n,J,Jdot]=robot.computeMSnJJdot;
    Pee(:,i) = robot.computeEnd_effectorPosition;
    dPee(:,i) = robot.computeEnd_effectorVelocity;
    error_p(:,i) = Kp*(X_d(i,:)'-Pee(1:3,i));
    error_v(:,i) = Kd*(dX_d(i,:)'-dPee(:,i));
    a = ddX_d(i,:)' + Kd*(dX_d(i,:)'-dPee(:,i)) + Kp*(X_d(i,:)'-Pee(1:3,i));
    ddq(:,i+1) = pinv(J)*a - pinv(J)*Jdot*dq(:,i);
    dq(:,i+1) = dq(:,i) + ddq(:,i)*dt;
    q(:,i+1) = q(:,i) + dq(:,i+1)*dt + 0.5*ddq(:,i)*dt*dt;
end

size_Pee = size(Pee,2);
figure;
plot(tSteps(1:size_Pee),Pee(1,:));
hold on
plot(tSteps(1:size_Pee),X_d(1:size_Pee,1));
xlabel('Time steps');
ylabel('Traj x');
legend('Output Traj x','Desired Traj x');

figure;
plot(tSteps(1:size_Pee),Pee(2,:));
hold on
plot(tSteps(1:size_Pee),X_d(1:size_Pee,2));
xlabel('Time steps');
ylabel('Traj y');
legend('Output Traj y','Desired Traj y');

figure;
plot(tSteps(1:size_Pee),Pee(3,:));
hold on
plot(tSteps(1:size_Pee),X_d(1:size_Pee,3));
xlabel('Time steps');
ylabel('Traj z');
legend('Output Traj z','Desired Traj z');

figure;
plot3(Pee(1,:)',Pee(2,:)',Pee(3,:)');;
hold on
plot3(X_d(:,1),X_d(:,2),X_d(:,3));
xlabel('Traj x');
ylabel('Traj y');
zlabel('Traj z');
legend('Output Traj','Desired Traj');