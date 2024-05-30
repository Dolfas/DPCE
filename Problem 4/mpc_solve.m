function u0 = mpc_solve(x0, H, R, A, B, C, ref,alpha)

    maxTemp = 55;

    n=size(A,1); % Dimension of the state
    
    % Compute the Augmented Cost matrices
    Q = C' * C;

    for k = 1:H
        if(k==1)
            Q_aug=Q;
        else
            Q_aug=blkdiag(Q_aug,Q);
        end
    end

    Q_aug = [zeros(size(Q_aug,1),n),Q_aug];
    Q_aug = [zeros(n,size(Q_aug,2));Q_aug];

    R_aug = R*eye(H,H); 

    
    % Sparse Formulation functions 
    F = 2 .* blkdiag(Q_aug,R_aug,alpha*eye(H,H)); 
    f = zeros(size(F,2),1);

    % Augmented Dynamics matrices
    for k = 1:H
        if(k==1)
            A_aug=A;
        else
            A_aug=blkdiag(A_aug,A);
        end
    end

    A_aug = [A_aug, zeros(size(A_aug,1),n)];
    A_aug = [zeros(n,size(A_aug,2)); A_aug];

    for k = 1:H
        if(k==1)
            B_aug=B;
        else
            B_aug=blkdiag(B_aug,B);
        end
    end
    B_aug = [zeros(n,size(B_aug,2));B_aug];
    
    for k = 1:H
        if(k==1)
            C_aug=C;
        else
            C_aug=blkdiag(C_aug,C);
        end
    end
    C_aug = [C_aug, zeros(size(C_aug,1),n)];

    E = [eye(n,n);zeros(H*n,n)];
    
    % Equality constraint matrices 
    Aeq = [A_aug - eye(size(A_aug,1),size(A_aug,2)), B_aug];
    Aeq = [Aeq, zeros(size(A_aug,1),H)]; 
    beq = -E*x0;
    
    % Inequality constraint matrices 
    g_aug = repmat(min([maxTemp-ref,0]),H,1);
    G_aug = eye(H,H);

    Aineq = [G_aug*C_aug, zeros(H,H), -eye(H,H)];
    bineq = g_aug; 

    % Lower and Higher bounds for the different variables 
    lb = -inf*ones(H*(n+2)+n,1); ub = inf*ones(H*(n+2)+n,1); % Bounds for the state
    lb(end-2*H+1:end-H,1)= -30; ub(end-2*H+1:end-H,1) = 68;  % Bounds for the control action
    lb(end-H+1:end,1)= 0; ub(end-H+1:end,1) = inf;           % Bounds for the eta

    % Set options to suppress quadprog output
    options = optimoptions('quadprog', 'Display', 'off');
    [z, ~, exitflag] = quadprog(F, f, Aineq, bineq, Aeq, beq, lb, ub, [], options);

    if exitflag == 1
        disp('Constrained optimization successful.');
    else
        disp('Constrained optimization failed.');
    end
    
    % Extract the optimal control action 
    u0=z(n*(H+1)+1);
end    