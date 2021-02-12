%% Some nice colors
Red       = [0.8500   0.3250   0.0980];
Yellow   = [0.929    0.694    0.125 ];
Blue     = [0        0.4470   0.7410];  
Violet   = [0.4940   0.1840   0.5560];
Green    = [0.4660   0.6740   0.1880];
Cyan     = [0.3010   0.7450   0.9330];
Bordeaux = [0.6350   0.0780   0.1840]; 
Silver = 1/255*[200,200,200];
%%
figure
plot3(Pee(1,:),Pee(2,:),Pee(3,:),'-','linewidth',1,'color',Blue),hold on
plot3(Pd(1,:),Pd(2,:),Pd(3,:),'--','linewidth',2,'color',Bordeaux),hold on
plot3(P_elbow(1,:),P_elbow(2,:),P_elbow(3,:),'-','linewidth',2,'color',Yellow),hold on
[X,Z] = meshgrid(0.06:0.06:0.42,-0.12:0.06:0.24);
Y = (XMin(axs)+0.003)*ones(size(X,1),size(X,1));    
s = surf(X,Y,Z);
s.EdgeColor = 'none';
s.FaceColor = Red;
s.FaceAlpha = 0.5;
% xlim([-0.5 0.5]);
% ylim([-0.7 0.3]);
% zlim([-0.2 0.8]);

xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
grid on;
legend('$\mathbf{X_{EE,d}}$','$\mathbf{X_{EE,cur}}$','$\mathbf{X_{elb,cur}}$','Obstacle','interpreter','latex','Location','southwest')
axis('normal')
% pause
% export_fig ('figures/exp2_fig1_origg','-eps','-pdf','-painters','-transparent')

figure
dT = (1:numSimSteps)*Ts;
yyaxis left
plot(dT,P_elbow(axs,1:numSimSteps),'-','linewidth',2)
% yline(XMax(axs),'--','linewidth',2,'color','#0072BD');
yline(XMin(axs),'--','linewidth',2,'color','#0072BD');
ylabel('Cartesian elbow position $[m]$','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
xlim([0 Tstop]);
ylim([min(P_elbow(axs,:))-0.1 max(P_elbow(axs,:))+0.1]);
grid on
yyaxis right
plot(dT,V_elbow(axs,1:numSimSteps),'linewidth',2)
% yline(VXMin(axs),'--','linewidth',2,'color','#D95319')
yline(VXMax(axs),'--','linewidth',2,'color','#D95319')
xlim([0 Tstop]);
ylim([min(V_elbow(axs,:))-0.1 max(V_elbow(axs,:))+0.1]);
ylabel('Cartesian elbow velocity $[m/s]$','Interpreter','latex');
legend('Current elbow position in y','Minimal elbow position in y','Current elbow velicity in y','Maximal elbow velocity in y','Interpreter','latex')

% export_fig ('figures/exp2_fig2_orig','-eps','-pdf','-painters','-transparent')

figure
dT = (1:numSimSteps)*Ts;
plot(dT,A_elbow(axs,1:numSimSteps)./AXMax(axs),'linewidth',2,'color',Blue),grid;
yline(1,'--','linewidth',2,'color','#0072BD')
yline(-1,'--','linewidth',2,'color','#0072BD')
xlim([0 Tstop]);
ylim([min(A_elbow(axs,:))/AXMax(axs)-1 max(A_elbow(axs,:))/AXMax(axs)+1]);
ylabel('Cartesian elbow acceleration $[\ddot{X}_{elb}/{A_{c,max}}]$','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
% export_fig ('figures/exp2_fig3_orig','-eps','-pdf','-painters','-transparent')
%%
figure
dT = (1:numSimSteps)*Ts;
plot(dT,q(:,1:numSimSteps)./QMax,'linewidth',2),grid on;
legend('joint 1','joint 2','joint 3','joint 4','joint 5','joint 6','Orientation','horizontal','Location','northoutside'	,'Interpreter','latex');
xlim([0 Tstop]);
ylim([-1.1 1.1]);
ylabel('Joint position $[q/Q_{max}]$','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
% export_fig ('figures/exp2_fig4_orig','-eps','-pdf','-painters','-transparent')

figure
dT = (1:numSimSteps)*Ts;
plot(dT,dq(:,1:numSimSteps)./VJntMax,'linewidth',2),grid;
legend('joint 1','joint 2','joint 3','joint 4','joint 5','joint 6','Orientation','horizontal','Location','northoutside'	,'Interpreter','latex');
xlim([0 Tstop]);
ylim([-1.1 1.1]);
ylabel('Joint velocity $[\dot{q}/V_{j,max}]$','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
% export_fig ('figures/exp2_fig5_orig','-eps','-pdf','-painters','-transparent')

figure
dT = (1:numSimSteps)*Ts;
plot(dT,ddq(:,1:numSimSteps)./AJntMax,'linewidth',2),grid;
legend('joint 1','joint 2','joint 3','joint 4','joint 5','joint 6','Orientation','horizontal','Location','northoutside'	,'Interpreter','latex');
xlim([0 Tstop]);
ylim([-1.1 1.1]);
ylabel('Joint acceleration $[\ddot{q}/A_{j,max}]$','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
% export_fig ('figures/exp2_fig6_orig','-eps','-pdf','-painters','-transparent')


figure
dT = (1:numSimSteps)*Ts;
plot(dT,tau_stack_SJS(:,1:numSimSteps),'linewidth',2),grid on;
legend('joint 1','joint 2','joint 3','joint 4','joint 5','joint 6','Orientation','horizontal','Location','northoutside'	,'Interpreter','latex');
xlim([0 Tstop]);
ylabel('Joint torque $[Nm]$','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
% export_fig ('figures/exp2_fig7_orig','-eps','-pdf','-painters','-transparent')


figure
dT = (1:numSimSteps)*Ts;
eP = Pd-Pee(1:3,1:size(Pee,2)-1);
plot(dT,eP(1,:),dT,eP(2,:),dT,eP(3,:),'linewidth',2),grid;
legend('$\mathbf{e_x}$','$\mathbf{e_y}$','$\mathbf{e_z}$','Interpreter','latex');
xlim([0 Tstop]);
ylim([-3.5e-1 3.5e-1]);
ylabel('End-effector position error [m]','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
% export_fig ('figures/exp2_fig8_orig','-eps','-pdf','-painters','-transparent')
