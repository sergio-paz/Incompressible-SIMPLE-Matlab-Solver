% =========================================================================
% Line-by-Line Tri-Diagonal Matrix Algorithm (TDMA) Solvers
% 
% The following functions provide a highly efficient computational suite for 
% solving line-implicit systems in two-dimensional fluid dynamics grids. 
% Utilising the Tri-Diagonal Matrix Algorithm (Thomas algorithm) via 
% forward elimination and backward substitution, these complementary routines 
% evaluate the discretised linear equations along continuous grid lines. 
% They are explicitly designed to perform alternating directional sweeps 
% (across both the i-axis and j-axis), a methodology that is instrumental 
% in expediting the convergence of the overarching numerical simulation.
% =========================================================================

% TDMA i
function u = tridag_i(d, a, b, c, e, r, u, n, l)
    % tridag_i Tri-diagonal matrix solver (derived from Numerical Recipes)
    % Solves the system along i-lines by iterating across the domain in the j-direction.
    
    % Preallocate the gam vector to optimise computational performance
    gam = zeros(n, 1); 
    
    % Iterate through the j-direction
    for j = 2:(l-1)
        
        % Forward recurrence phase
        bet = b(1, j);
        u(1, j) = r(1, j) / bet;
        
        for m = 2:n
            gam(m) = -c(m-1, j) / bet;
            bet = b(m, j) + a(m, j) * gam(m);
            rhs = r(m, j) + d(m, j) * u(m, j-1) + e(m, j) * u(m, j+1);
            u(m, j) = (rhs + a(m, j) * u(m-1, j)) / bet;
        end
        
        % Backward recurrence phase
        for m = (n-1):-1:1
            u(m, j) = u(m, j) - gam(m+1) * u(m+1, j);
        end
        
    end
end

% TDMA j
function u = tridag_j(d, a, b, c, e, r, u, n, l)
    % tridag_j Tri-diagonal matrix solver (derived from Numerical Recipes)
    % Solves the system along j-lines by iterating across the domain in the i-direction.

    % Preallocate the gam vector to optimise computational performance
    gam = zeros(n, 1); 
    
    % Iterate through the i-direction
    for i = 2:(l-1)
        
        % Forward recurrence phase
        bet = b(i, 1);
        u(i, 1) = r(i, 1) / bet;
        
        for m = 2:n
            gam(m) = -c(i, m-1) / bet;
            bet = b(i, m) + a(i, m) * gam(m);
            rhs = r(i, m) + d(i, m) * u(i-1, m) + e(i, m) * u(i+1, m);
            u(i, m) = (rhs + a(i, m) * u(i, m-1)) / bet;
        end
        
        % Backward recurrence phase
        for m = (n-1):-1:1
            u(i, m) = u(i, m) - gam(m+1) * u(i, m+1);
        end
        
    end
end