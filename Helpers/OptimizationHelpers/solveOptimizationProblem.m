function deta_dt = solveOptimizationProblem(H,Aeq,beq)
    options = optimoptions('quadprog', ...
        'Algorithm', 'trust-region-reflective');
    
    n = size(H,1) / 2;
    deta_dt_stack = quadprog(H, zeros(2*n,1), [], [], Aeq, beq, [], [], [], options);
    deta_dt_x = deta_dt_stack(1:n);
    deta_dt_y = deta_dt_stack(n+1:end);
    deta_dt = complex(deta_dt_x, deta_dt_y);
end