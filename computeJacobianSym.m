function [J, Jdot] = computeJacobianSym(alpha,a,d,theta,dq)
    for i = 1:length(alpha)
        M_ind(:,:,i) = DHmatrix(alpha(i),a(i),d(i),theta(i));
    end
    M_net = (M_ind(:,:,1)*M_ind(:,:,2)*M_ind(:,:,3)*M_ind(:,:,4)*M_ind(:,:,5)*M_ind(:,:,6)*M_ind(:,:,7));

    p_all = M_net*[0;0;0;1];
    p = p_all(1:3);
    J1 = diff(p,theta(1));
    J2 = diff(p,theta(2));
    J3 = diff(p,theta(3));
    J4 = diff(p,theta(4));
    J5 = diff(p,theta(5));
    J6 = diff(p,theta(6));
    %J7 = diff(p,theta(7));
    J = [J1 J2 J3 J4 J5 J6];

    Jd_1 = diff(J,theta(1))*dq(1);
    Jd_2 = diff(J,theta(2))*dq(2);
    Jd_3 = diff(J,theta(3))*dq(3);
    Jd_4 = diff(J,theta(4))*dq(4);
    Jd_5 = diff(J,theta(5))*dq(5);
    Jd_6 = diff(J,theta(6))*dq(6);
    J_dot = Jd_1 + Jd_2 + Jd_3 + Jd_4 + Jd_5 + Jd_6; %Manually observe to take common and form brackets
    Jdot = simplify(J_dot);
end