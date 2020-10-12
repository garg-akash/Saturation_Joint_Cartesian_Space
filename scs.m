function [tau_stack,tasks_scs] = scs(M, c, g, dq, J, Jdot, m, nJoints, X, dX, dt, X_max, X_min, dX_max,...
                 dX_min, ddX_max, ddX_min, tasks, num_iter, JCar, JdotCar, XCar, dXCar, nCar)
disp('11')
n = c+g;
dq_=dq;
% I = eye(6);
% Z = zeros(6);
iter = 1;
tol = 1e-3;
%J_lim = [];
%Jdot_lim = [];
J_carLim = [];
J_carLim_dot = [];
%ddq_sat = [];
ddX_sat = [];
len_sat_prev = 0;
%flag_overwrite = false;
flag_sat = 1;
while (flag_sat)
    if iter>5
        keyboard;
        tau_stack = tau_stackTmp(:,iter-1);
        tasks_scs = tasks;
        break;
    end
    %tau_stackTmp(:,iter) = getTorques_revision(tasks,M,n,dq_); %%  dq_
    tau_stackTmp(:,iter) = getTorques(tasks,M,c,g,dq,nJoints); 
    ddq(:,iter) = M^(-1)*(tau_stackTmp(:,iter) - n);
    %dq_ = dq + ddq(:,iter)*dt;
    %q_ = q + dq*dt + 0.5*ddq(:,iter)*dt.^2;
    ddX(:,iter) = JCar*ddq(:,iter) + JdotCar*dq;
    %dXCar_ = JCar*dq_;
    %XCar_ = XCar + dXCar_*dt + 0.5*ddX(:,iter)*dt.^2;
    
    V_min = max([(X_min-XCar)/dt,dX_min,-sqrt(2*ddX_max.*abs(XCar-X_min))],[],2);
    V_max = min([(X_max-XCar)/dt,dX_max,sqrt(2*ddX_max.*abs(X_max-XCar))],[],2);
    
    A_min = max([2*(X_min-XCar-dXCar*dt)/dt^2, (V_min-dXCar)/dt, ddX_min],[],2);
    A_max = min([2*(X_max-XCar-dXCar*dt)/dt^2, (V_max-dXCar)/dt, ddX_max],[],2);    
    
    idx0 = find((A_max<A_min) & (A_max<ddX_min));
    if idx0>0
        for x=1:size(idx0,1)
            if (A_max(idx0(x))<A_min(idx0(x)))
                A_max(idx0(x))=A_min(idx0(x));
            end
        end
    end
    idx1 = find((A_min>A_max) & (A_min>ddX_max));
    if idx1>0
        for x=1:size(idx1,1)
            if (A_min(idx1(x))>A_max(idx1(x)))
                A_min(idx1(x))=A_max(idx1(x));
            end
        end
    end
    
    %totNumSat = nnz((ddX(:,iter) < (ddXMin - tol)))+nnz((ddX(:,iter) > (ddXMax + tol)));

    
    if all( ddX(3,iter)>(A_min-tol) ) && all( ddX(3,iter)<(A_max+tol) )
        flag_sat = false;
    else
        flag_sat = true;
        for i = 3:m
            if (ddX(i,iter) < (A_min(i) - tol))
%                 if (A_min(i) < A_max(i))
%                     ddX_sat = [ddX_sat;A_min(i)];
%                 else
%                     if (A_max(i)>ddX_min(i))
%                         ddX_sat = [ddX_sat;A_max(i)];
%                     else
%                         ddX_sat = [ddX_sat;ddX_min(i)];
%                     end
%                 end
                ddX_sat = [ddX_sat;A_min(i)];
                J_carLim = [J_carLim;JCar(i,:)];
%                 J_lim = [J_lim;I(i,:)];
%                 Jdot_lim = [Jdot_lim;Z(i,:)];
                J_carLim_dot = [J_carLim_dot;JdotCar(i,:)];
            elseif (ddX(i,iter) > (A_max(i) + tol))
%                 if (A_max(i) > A_min(i))
%                     ddX_sat = [ddX_sat;A_max(i)];
%                 else
%                     if (A_min(i)<ddX_max(i))
%                         ddX_sat = [ddX_sat;A_min(i)];
%                     else
%                         ddX_sat = [ddX_sat;ddX_max(i)];
%                     end
%                 end
                ddX_sat = [ddX_sat;A_max(i)];
                J_carLim = [J_carLim;JCar(i,:)];
%                 J_lim = [J_lim;I(i,:)];
%                 Jdot_lim = [Jdot_lim;Z(i,:)];
                J_carLim_dot = [J_carLim_dot;JdotCar(i,:)];
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
    if (length(ddX_sat)>len_sat_prev)
        if (len_sat_prev>0)
            tasks(1) = []; %the new task contains the previous saturations so remove the previous task
        end
        task(1).J = J_carLim;
        task(1).f_star = ddX_sat;
        task(1).Jdot = J_carLim_dot;
        tasks = [task;tasks];
        len_sat_prev = length(ddX_sat);
    end
    iter = iter + 1;
end

if (flag_sat == false)
    tau_stack = tau_stackTmp(:,iter-1);
    tasks_scs = tasks;
end

end