function [T_c] = scs(M, c, g, dq, J, Jdot, m, n, X, dX, dt, X_max, X_min, V_max,...
                 V_min, A_max, A_min, tasks, iter, JCar, JdotCar, XCar, dXCar, nCar)
    disp('11')
    for i=1:m
        dX_max(i) = min(V_max(i),sqrt(2*A_max(i)*(X_max(i)-XCar(i))));
        dX_min(i) = max(V_min(i),-sqrt(2*A_max(i)*(X_max(i)-XCar(i))));
    end
    flag_sat = 1;
    J_carLim = [];
    J_carLim_dot = [];
    ddX_sat = [];
    len_sat_prev = 0;
%     if (iter==1473)
%         keyboard;
%     end
    %%%%Start loop%%%%
    while flag_sat
        T_stack = getTorques(tasks, M, c, g, dq, m, n);
        ddq = M\(T_stack - g - c);
        ddX = JCar*ddq + JdotCar*dq;
        ddX_max = min(min(2*(X_max-XCar-dXCar*dt)/power(dt,2), (dX_max'-dXCar)/dt), A_max);
        ddX_min = max(max(2*(X_min-XCar-dXCar*dt)/power(dt,2), (dX_min'-dXCar)/dt), A_min);
        
%         J_carLim = [];
%         J_carLim_dot = [];
%         ddX_sat = [];
%         len_sat_prev = 0;
        for i = 1:m
            if( (ddX(i)>(ddX_max(i)+0.005)) || (ddX(i)<(ddX_min(i)-0.005)) )
                disp('saturation found!');
                J_carLim = [J_carLim;JCar(i,:)];
                J_carLim_dot = [J_carLim_dot;JdotCar(i,:)];
                if(ddX(i)>ddX_max(i))
                    ddX_sat = [ddX_sat;ddX_max(i)];
                else
                    ddX_sat = [ddX_sat;ddX_min(i)];
                end
            end
        end
        if (length(ddX_sat)>len_sat_prev)
            if (len_sat_prev>0)
                tasks(1) = []; %the new task contains the previous saturations so remove the previous task
            end
            task.J = J_carLim;
            task.Jdot = J_carLim_dot;
            task.f_star = ddX_sat;  
            tasks = [task;tasks];
            len_sat_prev = length(ddX_sat);
            %%%%back to start of loop%%%%
        else
            flag_sat = 0;
            %%%%end loop%%%%
        end
    end
    T_c = T_stack;
end
