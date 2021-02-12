function tau_stack = getTorques(task,M,n,dq)
global tol
tau_stack = zeros(6,1);
N_aug = cell(1,size(task,2));
J_hat = cell(1,size(task,2));
J_hat_inv = cell(1,size(task,2));
tau_p = zeros(6,size(task,2));
i=1;
while(i<=size(task,2))
    J = task(i).J;
    Jdot = task(i).Jdot;
    F_ast = task(i).F_ast;
    if i==1
        N_aug{i} = eye(6);
    else
        N_aug{i} = N_aug{i-1}*(eye(6)-J_hat{i-1}'*J_hat_inv{i-1}');
    end
    J_hat{i} = J*N_aug{i}';
    Jdot_proj = Jdot*N_aug{i}';
    [U,S,~] = svd(J_hat{i}*M^(-1)*J_hat{i}');
    if size(S,2)>1
        flagOk = true;
        for j = size(S,2):-1:1
            if S(j,j)<tol
                U(:,j) = []; S(:,j) = []; S(j,:) = [];
            end
        end
    else
        if S<tol
            flagOk = false;
            task(i) = [];
        else
            flagOk = true;
        end
    end
    if (flagOk)
        A_hat = U*S^(-1)*U';
        J_hat_inv{i} = M^(-1)*J_hat{i}'*A_hat;
        n_hat = J_hat_inv{i}'*n-A_hat*Jdot_proj*dq;
        F_hat = A_hat*F_ast+n_hat;
        tau = J'*F_hat;
        tau_p(:,i) = N_aug{i}*tau;
        tau_stack = tau_p(:,i) + tau_stack;
        i = i + 1;
    end
end
end