function task = ReOrderTasks(task,J_LIM,Jdot_LIM,ddqSat, flag_overwrite)
if (flag_overwrite==false) %shift the tasks in the first iteration
    nTasks = size(task,2);
    for i = nTasks:-1:1
        task(i+1) = task(i);
    end
end
task(1).J = J_LIM; %overwrite the joint limiting task (highest priority) in the next iterations
task(1).F_ast = ddqSat;
task(1).Jdot = Jdot_LIM;
end