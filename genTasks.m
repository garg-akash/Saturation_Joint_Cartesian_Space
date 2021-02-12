function task = genTasks(J,Jdot,Pd,Pee,dPd,Vee,dq) %two primary tasks
global Kp Kv kd
task(1).J = J;
task(1).Jdot = Jdot;
task(1).F_ast = Kp*(Pd-Pee)+Kv*(dPd-Vee);
task(2).J = eye(6,6);
task(2).Jdot = zeros(6,6);
task(2).F_ast = -kd*dq;
end