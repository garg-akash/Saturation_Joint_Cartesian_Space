function [tau_stack,tasks_sjs] = sjs(M, c, g, dq, J, Jdot, m, nJoints, dt, q, Q_max, Q_min, dQ_max, dQ_min, ddQ_max, ddQ_min, tasks, num_iter)
disp('11')
n = c+g;
% I = eye(6);
% Z = zeros(6);
iter = 1;
tol = 1e-3;
J_lim = [];
Jdot_lim = [];
ddq_sat = [];
len_sat_prev = 0;
flag_overwrite = false;
flag_sat = 1;
% if num_iter==474
%     keyboard;
% end
while (flag_sat)
    if iter>5
        %keyboard;
        tau_stack = tau_stackTmp(:,iter-1);
        tasks_sjs = tasks;
        break;
    end
    %tau_stackTmp(:,iter) = getTorques_revision(tasks,M,n,dq); 
    tau_stackTmp(:,iter) = getTorques(tasks,M,c,g,dq,nJoints); 
    ddq(:,iter) = M^(-1)*(tau_stackTmp(:,iter) - n);
    %dq_ = dq + ddq(:,iter)*dt;
    %q_ = q + dq*dt + 0.5*ddq(:,iter)*dt.^2;
    
    V_min = max([(Q_min-q)/dt,dQ_min,-sqrt(2*ddQ_max.*abs(q-Q_min))],[],2);
    V_max = min([(Q_max-q)/dt,dQ_max,sqrt(2*ddQ_max.*abs(Q_max-q))],[],2);
    
    A_min = max([2*(Q_min-q-dq*dt)/dt^2, (V_min-dq)/dt, ddQ_min],[],2);
    A_max = min([2*(Q_max-q-dq*dt)/dt^2, (V_max-dq)/dt, ddQ_max],[],2);    
    
    idx0 = find(A_max<A_min & A_max<ddQ_min);
    if idx0>0
        for x=1:size(idx0,1)
            if (A_max(idx0(x))<A_min(idx0(x)))
                A_max(idx0(x))=A_min(idx0(x));
            end
        end
    end
    idx1 = find((A_min>A_max) & (A_min>ddQ_max));
    if idx1>0
        for x=1:size(idx1,1)
            if (A_min(idx1(x))>A_max(idx1(x)))
                A_min(idx1(x))=A_max(idx1(x));
            end
        end
    end
    
    if all( ddq(:,iter)>(A_min-tol) ) && all( ddq(:,iter)<(A_max+tol) )
        flag_sat = false;
    else
        flag_sat = true;
        for i = 1:nJoints
            if (ddq(i,iter) < (A_min(i) - tol))
%                 if (A_min(i) < A_max(i))
%                     ddq_sat = [ddq_sat;A_min(i)];
%                 else
%                     if (A_max(i)>ddQ_min(i))
%                         ddq_sat = [ddq_sat;A_max(i)];
%                     else
%                         ddq_sat = [ddq_sat;ddQ_min(i)];
%                     end
%                 end
                ddq_sat = [ddq_sat;A_min(i)];
                J_lim_temp = zeros(1,nJoints);
                J_lim_temp(i) = 1;
                J_lim = [J_lim;J_lim_temp];
%                 J_lim = [J_lim;I(i,:)];
%                 Jdot_lim = [Jdot_lim;Z(i,:)];
                Jdot_lim = zeros(size(J_lim));
            elseif (ddq(i,iter) > (A_max(i) + tol))
%                 if (A_max(i) > A_min(i))
%                     ddq_sat = [ddq_sat;A_max(i)];
%                 else
%                     if (A_min(i)<ddQ_max(i))
%                         ddq_sat = [ddq_sat;A_min(i)];
%                     else
%                         ddq_sat = [ddq_sat;ddQ_max(i)];
%                     end
%                 end
                ddq_sat = [ddq_sat;A_max(i)];
                J_lim_temp = zeros(1,nJoints);
                J_lim_temp(i) = 1;
                J_lim = [J_lim;J_lim_temp];
%                 J_lim = [J_lim;I(i,:)];
%                 Jdot_lim = [Jdot_lim;Z(i,:)];
                Jdot_lim = zeros(size(J_lim));
            end
        end
%     if iter > 1
%         flag_overwrite = true;  
%         [J_lim,ia,~]=unique(J_lim,'rows','last');
%          ddq_sat = ddq_sat(ia);
%          Jdot_lim = Jdot_lim(ia,:);
%     end
    end
%     if (flag_sat)
%         tasks = reOrderTasks_revision(tasks, J_lim, Jdot_lim, ddq_sat, flag_overwrite);
%     end
    if (length(ddq_sat)>len_sat_prev)
        if (len_sat_prev>0)
            tasks(1) = []; %the new task contains the previous saturations so remove the previous task
        end
        task(1).J = J_lim;
        task(1).f_star = ddq_sat;
        task(1).Jdot = Jdot_lim;
        tasks = [task;tasks];
        len_sat_prev = length(ddq_sat);
    end
    iter = iter + 1;
end

if (flag_sat == false)
    tau_stack = tau_stackTmp(:,iter-1);
    tasks_sjs = tasks;
end

end