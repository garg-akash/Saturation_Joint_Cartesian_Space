%n is the no of joints in the robot
%J_h = J hat.......Jacobian i projected onto the null space of all higher
%proiority tasks
%J_h = J hat
%J_hc = J hat cap
function [T_stack] = getTorques(tasks, M, c, g, dq, nJoint)%qsym, dqsym, symvalues,
    disp('22')
    n = c+g;
    k = length(tasks);
    %fprintf('Number of tasks given : %d',k);
    
    T_proj = zeros(nJoint,k);
    %T_proj = [];
    i = 1;
    while(i<=length(tasks))
        J = tasks(i).J;
        f_star = tasks(i).f_star;
        %Jsym = tasks(i).Jsym;
        Jdot = tasks(i).Jdot;
        if i==1
            N_aug = eye(nJoint);
        else
            N_aug = N_aug_prev*(eye(nJoint) - J_h_prev'*J_hc_prev');
        end
%         Jsym_h = Jsym*N_aug';
%         %%%%%%%%%%%%%%
%         J_h_dot1 = diff(Jsym_h,qsym(1))*dqsym(1);
%         J_h_dot2 = diff(Jsym_h,qsym(2))*dqsym(2);
%         J_h_dot3 = diff(Jsym_h,qsym(3))*dqsym(3);
%         J_h_dot4 = diff(Jsym_h,qsym(4))*dqsym(4);
%         J_h_dot5 = diff(Jsym_h,qsym(5))*dqsym(5);
%         J_h_dot6 = diff(Jsym_h,qsym(6))*dqsym(6);
%         J_h_dot_sym = J_h_dot1 + J_h_dot2 + J_h_dot3 + J_h_dot4 + J_h_dot5 + J_h_dot6;
%         %%%%%%%%%%%%%%       
%         J_h = double(subs(Jsym_h))
%         J_h_dot = double(subs(J_h_dot_sym))
        J_h = J*N_aug';
        J_h_dot = Jdot*N_aug';
        [U,S,V] = svd(J_h*(M\J_h')); %svd of lambda_inverse
%         if (i==1 && k>2)
%             for j = size(S,2):-1:1
%                 if S(j,j)<0.001
%                     U(:,j) = [];
%                     S(:,j) = [];
%                     S(j,:) = [];
%                 end
%             end
%         else
%         for j = size(S,2):-1:1
%             if S(j,j)<0.1
%                 U(:,j) = [];
%                 S(:,j) = [];
%                 S(j,:) = [];
%             end
%         end
%         end
        if size(S,2)>1
            flagOk = true;
            for j = size(S,2):-1:1
                if S(j,j)<0.001
                    U(:,j) = []; S(:,j) = []; S(j,:) = [];
                end
            end
        else
            if S<0.001
                flagOk = false;
                tasks(i) = [];
            else
                flagOk = true;
            end
        end
        if(flagOk)
            lambda_h = U*S^(-1)*U';
            J_hc = M\(J_h'*lambda_h);
            mu_h = J_hc'*c - lambda_h*J_h_dot*dq;
            p_h = J_hc'*g;
            f_h = lambda_h*f_star + mu_h + p_h;
            T = J'*f_h;
            T_proj(:,i) = N_aug*T;
            i = i+1;
            N_aug_prev = N_aug;
            J_h_prev = J_h;
            J_hc_prev = J_hc;
        end
    end
    T_stack = sum(T_proj,2);
end