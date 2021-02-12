function [tasksOutputSCS] = scs(M, n, J, Jdot , dq, X, dX, XMin, XMax, VXMin, VXMax, AXMin, AXMax, tasks)
global tol Ts
[mTsk, ~] = size(J);
ddXSat = [];iter = 1;
J_CarLIM=[];Jdot_CarLIM=[];
LimitExceed = true;
flag_overwrite = false;
while (LimitExceed)
    if iter>3
        tasksOutputSCS = tasks; %return the cartesian limiting task
        break;
    end
    tau_stackTmp(:,iter) = getTorques(tasks,M,n,dq); %compute torques with Algorithm1
    ddq(:,iter) = M^(-1)*(tau_stackTmp(:,iter) - n); %compute the resulting joint acceleration
    ddX(:,iter) = J*ddq(:,iter) + Jdot*dq; %compute the resulting elbow acceleration
    
    %compute cartesian acceleration bounds
    dXMin = max([(XMin-X)/Ts,VXMin,-sqrt(2*AXMax.*abs(X-XMin))],[],2);
    dXMax = min([(XMax-X)/Ts,VXMax,sqrt(2*AXMax.*abs(XMax-X))],[],2);
    
    ddXMin = max([2*(XMin-X-dX.*Ts)./Ts.^2, (dXMin-dX)./Ts, AXMin],[],2);
    ddXMax = min([2*(XMax-X-dX.*Ts)./Ts.^2, (dXMax-dX)./Ts, AXMax],[],2);
    
    idx0 = find((ddXMax<ddXMin) & (ddXMax<AXMin));
    if idx0>0
        for x=1:size(idx0,1)
            if (ddXMax(idx0(x))<ddXMin(idx0(x)))
                ddXMax(idx0(x))=ddXMin(idx0(x));
            end
        end
    end
    idx1 = find((ddXMin>ddXMax) & (ddXMin>AXMax));
    if idx1>0
        for x=1:size(idx1,1)
            if (ddXMin(idx1(x))>ddXMax(idx1(x)))
                ddXMin(idx1(x))=ddXMax(idx1(x));
            end
        end
    end
    
    %check if the accelerations are within their bounds    
    if all( ddX(:,iter)>(ddXMin-tol) ) && all( ddX(:,iter)<(ddXMax+tol) )
        LimitExceed = false;
    else
        LimitExceed = true;
        for i = 1:mTsk 
            if (ddX(i,iter) < (ddXMin(i) - tol))
                ddXSat = [ddXSat;ddXMin(i)];
                J_CarLIM = [J_CarLIM;J(i,:)];
                Jdot_CarLIM = [Jdot_CarLIM;Jdot(i,:)];
            elseif (ddX(i,iter) > (ddXMax(i) + tol))
                ddXSat = [ddXSat;ddXMax(i)];
                J_CarLIM = [J_CarLIM;J(i,:)];
                Jdot_CarLIM = [Jdot_CarLIM;Jdot(i,:)];
            end
        end
        if iter > 1
            flag_overwrite = true;
        end
    % reorder the tasks and set the cartesian limiting task as the highest priority
    tasks = ReOrderTasks(tasks, J_CarLIM, Jdot_CarLIM, ddXSat, flag_overwrite);
    end
    iter = iter + 1;
end
    tasksOutputSCS = tasks; %return the tasks containing the cartesian limiting task
    % this will become the input of the SJS algorithm
end