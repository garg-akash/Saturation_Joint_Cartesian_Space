function [tau_stack] = SJCS(M, n, J, J_ev, dJ_ev, q, dq, QMin, QMax, VJntMin, VJntMax, AJntMin, AJntMax, X, dX, XMin, XMax, VXMin, VXMax, AXMin, AXMax, tasks)
global tol
iter = 1;
[mTsk, nJnt] = size(J);
I = eye(6);Z = zeros(6);
ddXSat = [];ddqSat = [];
J_LIM = [];Jdot_LIM = [];
J_CarLIM = [];Jdot_CarLIM = [];
WJnt = eye(nJnt);WCar = eye(mTsk);
LimitExceed=true;flag_overwrite=false;

while (LimitExceed)
    if iter>9 
        tau_stack = tau_stackTmp(:,iter-1); %return the torque
        break;
    end
    tau_stackTmp(:,iter) = getTorques(tasks,M,n,dq); %compute torques with Algorithm1
    ddq(:,iter) = M^(-1)*(tau_stackTmp(:,iter) - n); %compute the resulting joint acceleration
    ddX(:,iter) = J_ev*ddq(:,iter) + dJ_ev*dq; %compute the resulting elbow acceleration
    
    %Compute the joint and cartesian bounds
    [ddXMin,ddXMax] = cartesianAccelerationBounds(X,dX,XMin,XMax,VXMin,VXMax,AXMin,AXMax);
    [ddQMin,ddQMax] = jointAccelerationBounds(q,dq,QMin,QMax,VJntMin,VJntMax,AJntMin,AJntMax);
    
    %check if the accelerations are within their bounds
    if all( ddq(:,iter)>(ddQMin-tol) ) && all( ddq(:,iter)<(ddQMax+tol) )...
            && all( ddX(:,iter)>(ddXMin-tol) ) && all( ddX(:,iter)<(ddXMax+tol) )
        LimitExceed = false;
    else 
        LimitExceed = true;
        % get the acceleration ratios
        [sJnt,sCar] = getAccelerationRatios(WJnt,WCar,ddq(:,iter),ddQMin,ddQMax,ddX(:,iter),ddXMin,ddXMax);
        
        if min(sJnt) <= min(sCar) %the critical case corresponds to a joint acc.
            [~, jntIdx] = min(sJnt);
            J_ast = [J_LIM;I(jntIdx,:)];
            WJnt(jntIdx, jntIdx) = 0;
            if rank(J_ast)>rank(J_LIM)
                J_LIM = J_ast;
                Jdot_LIM = [Jdot_LIM;Z(jntIdx,:)];
                ddqSat = [ddqSat;min(max(ddQMin(jntIdx), ddq(jntIdx,iter)), ddQMax(jntIdx))];
            end
        else %the critical case corresponds to a cartesian acc.
            [~, carIdx] = min(sCar);
            J_ast = [J_LIM;J_ev(carIdx,:)];
            WCar(carIdx, carIdx) = 0;
            if rank(J_ast)>rank(J_LIM)
                J_LIM = J_ast;
                Jdot_LIM = [Jdot_LIM;dJ_ev(carIdx,:)];
                ddqSat = [ddqSat;min(max(ddXMin(carIdx), ddX(carIdx,iter)), ddXMax(carIdx))];
            end
        end
        if iter > 1
            flag_overwrite = true;
        end
    % reorder the tasks and set the limiting task as the highest priority
    tasks = ReOrderTasks(tasks, J_LIM, Jdot_LIM, ddqSat, flag_overwrite);
    end
    iter = iter + 1;
end
% return the torque
tau_stack = tau_stackTmp(:,iter-1);

end