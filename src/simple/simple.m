% =========================================================================
% SIMPLE Algorithm Iterative Solver
% 
% This section executes the primary iterative routine of the Semi-Implicit 
% Method for Pressure Linked Equations (SIMPLE). Within each iteration, it 
% successively resolves the momentum equations to yield intermediate velocity 
% fields, evaluates the pressure-correction equation, and subsequently 
% applies these corrections to update the pressure and velocity fields. 
% Furthermore, a reference pressure is established to prevent matrix 
% singularity. The iterative cycle persists until the maximum computational 
% residual across the system descends below a predefined tolerance, thereby 
% confirming mathematical convergence.
% =========================================================================

n = 1; % Initialise N Iterrations
maxRes = 1000; % Initial high residual

while (n<=max_iter) && (maxRes>=tol)
    
    % Solve U-Momentum
    [U_star, dU, uRes] = solve_UMomentum(U, V, P, dx, dy, Re, alpha_u, Nx, Ny);

    % Solve V-Momentum
    [V_star, dV, vRes] = solve_VMomentum(U, V, P, dx, dy, Re, alpha_v, Nx, Ny);

    % Solve Pressure-Correction Equation
    [P_prime, pRes] = solve_PCorrection(U_star, V_star, dU, dV, dx, dy, Nx, Ny);

    % Correct Velocities and Pressure
    [U_new, V_new, P_new] = correct_Fields(U_star, V_star, P, P_prime, dU, dV, alpha_p, Nx, Ny);

    % Normalize pressure (pin one node to 0 to set a reference)
    P_new = P_new - P_new(2,2); 
    
    % Update fields for next iteration
    U = U_new;
    V = V_new;
    P = P_new;
    
    % Check maximum residual
    maxRes = max([uRes, vRes, pRes]);
    fprintf('Iter: %d | Res: %.4e \n', n, maxRes);

    % Update iteration
    n = n + 1;

end % END LOOP

% Check convergence
if maxRes < tol
    fprintf('Convergence achieved at iteration %d. Final residual: %.4e\n', n-1, maxRes);
end