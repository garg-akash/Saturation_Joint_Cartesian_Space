function [J, Jdot] = computeJacobianSym(alpha,a,d,theta,dq)
    M_net = 1;
    for i = 1:length(alpha)
        M_ind(:,:,i) = DHmatrix(alpha(i),a(i),d(i),theta(i));
        M_net = M_net*M_ind(:,:,i);
    end
    p_all = M_net*[0;0;0;1];
    p = p_all(1:3);
    
    DIFF_LEN = length(alpha);
    J_dot = 0;
    if i==7
        DIFF_LEN = 6; %We have fixed joint 7 and don't want to differentiate
    end
    for i = 1:DIFF_LEN
        J(:,i) = diff(p,theta(i));
    end
    
    for i = 1:DIFF_LEN
        Jd = diff(J,theta(i))*dq(i);
        J_dot = J_dot + Jd;
    end
    Jdot = simplify(J_dot);
end