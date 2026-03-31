% =========================================================================
% Field Correction: correct_Fields
% 
% This routine executes the final phase of the SIMPLE algorithm iteration 
% by updating the kinematic and thermodynamic fields. It  applies the
% derived pressure correction gradients to the intermediate velocity
% components, thereby enforcing strict mass conservation (continuity)
% across the computational domain. Furthermore, it computes the newly
% corrected pressure field by incorporating an explicit under-relaxation
% factor, ensuring numerical stability for subsequent iterative cycles.
% =========================================================================

function [U_new, V_new, P_new] = correct_Fields(U_star, V_star, P, P_prime, dU, dV, alpha_p, Nx, Ny)
    % correct_Fields Updates the velocity and pressure fields
    % This function applies the calculated pressure correction (P_prime) to enforce 
    % mass conservation (continuity) across the computational domain.
    
    % Initialize outputs with the guessed fields
    U_new = U_star;
    V_new = V_star;
    P_new = P;
    
    % Grid dimensions
    NxU = Nx + 1; NyU = Ny + 2;
    NxV = Nx + 2; NyV = Ny + 1;
    
    % Correct U velocities
    for i = 2:(NxU-1)
        for j = 2:(NyU-1)
            U_new(i,j) = U_star(i,j) - dU(i,j) * (P_prime(i+1, j) - P_prime(i, j));
        end
    end
    
    % Correct V velocities
    for i = 2:(NxV-1)
        for j = 2:(NyV-1)
            V_new(i,j) = V_star(i,j) - dV(i,j) * (P_prime(i, j+1) - P_prime(i, j));
        end
    end
    
    % Correct Pressure 
    P_new = P + alpha_p * P_prime;
end