function [tau_stack] = sjs(M, n, J, q, dq, QMin, QMax, VJntMin, VJntMax, AJntMin, AJntMax, tasks)
global tol Ts
[~, nJnt] = size(J);
ddqSat = [];iter = 1;
I = eye(6);Z = zeros(6);
J_LIM = [];Jdot_LIM = [];
LimitExceed = true;
flag_overwrite = false;

while (LimitExceed)
    if iter>6 
        tau_stack = tau_stackTmp(:,iter-1); %return the torque
        break;
    end
    tau_stackTmp(:,iter) = getTorques(tasks,M,n,dq); %compute torques with Algorithm1
    ddq(:,iter) = M^(-1)*(tau_stackTmp(:,iter) - n); %compute the resulting joint acceleration
    
    %compute joint acceleration bounds
    dQMin = max([(QMin-q)/Ts,VJntMin,-sqrt(2*AJntMax.*abs(q-QMin))],[],2);
    dQMax = min([(QMax-q)/Ts,VJntMax,sqrt(2*AJntMax.*abs(QMax-q))],[],2);
     
    ddQMin = max([2*(QMin-q-dq*Ts)/Ts^2, (dQMin-dq)/Ts, AJntMin],[],2);
    ddQMax = min([2*(QMax-q-dq*Ts)/Ts^2, (dQMax-dq)/Ts, AJntMax],[],2);
    
    idx0 = find(ddQMax<ddQMin & ddQMax<AJntMin);
    if idx0>0
        for x=1:size(idx0,1)
            if (ddQMax(idx0(x))<ddQMin(idx0(x)))
                ddQMax(idx0(x))=ddQMin(idx0(x));
            end
        end
    end
    idx1 = find((ddQMin>ddQMax) & (ddQMin>AJntMax));
    if idx1>0
        for x=1:size(idx1,1)
            if (ddQMin(idx1(x))>ddQMax(idx1(x)))
                ddQMin(idx1(x))=ddQMax(idx1(x));
            end
        end
    end
    %check if the accelerations are within their bounds  
    if all( ddq(:,iter)>(ddQMin-tol) ) && all( ddq(:,iter)<(ddQMax+tol) )
        LimitExceed = false;
    else
        LimitExceed = true;
        for i = 1:nJnt
            if (ddq(i,iter) < (ddQMin(i) - tol))
                ddqSat = [ddqSat;ddQMin(i)];
                J_LIM = [J_LIM;I(i,:)];
                Jdot_LIM = [Jdot_LIM;Z(i,:)];
            elseif (ddq(i,iter) > (ddQMax(i) + tol))
                ddqSat = [ddqSat;ddQMax(i)];
                J_LIM = [J_LIM;I(i,:)];
                Jdot_LIM = [Jdot_LIM;Z(i,:)];
            end
        end
        if iter > 1
            flag_overwrite = true;
        end     
    % reorder the tasks and set the joint limiting task as the highest priority
    tasks = ReOrderTasks(tasks, J_LIM, Jdot_LIM, ddqSat, flag_overwrite);
    end
    iter = iter + 1;
end
% return the torque
tau_stack = tau_stackTmp(:,iter-1);

end