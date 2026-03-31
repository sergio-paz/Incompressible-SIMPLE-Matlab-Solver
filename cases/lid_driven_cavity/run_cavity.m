% =========================================================================
% Incompressible Flow Solver for lid-driven cavity flow problem
% FVM SIMPLE-TDMA
% Structured Uniform Staggered Grid
% =========================================================================
clc; clear all; close all;

%% Read parameters
fileID = fopen('params.txt', 'r');
if fileID == -1
    error('Failed to open the file.');
end

data = textscan(fileID, '%f %*[^\n]'); % Reading only the number (skip everything else)
fclose(fileID);
params = data{1}; % Create a vector with the parameters


%% Problem Setup
H              = params(1);    % Domain Height
W              = params(2);    % Domain Width
Re             = params(3);    % Reynolds Number
Nx             = params(4);    % Number of cells in x-direction
Ny             = params(5);    % Number of cells in y-direction
max_iter       = params(6);    % Maximum Iterations (Corregido: Maximum)
tol            = params(7);    % Convergence criteria
alpha_p        = params(8);    % P under-relaxation
alpha_u        = params(9);    % U under-relaxation
alpha_v        = params(10);   % V under-relaxation

% Print Parameters
fprintf('--- Parameters loaded successfully ---\n'); % Corregido: loaded / successfully
fprintf('Domain: %.2f x %.2f | Re: %.1f\n', W, H, Re);
fprintf('Grid: %d x %d\n', Nx, Ny);
fprintf('Convergence: %.e\n', tol);
fprintf('---------------------------------------\n');


%% Spatial Discretisation (Structured Uniform Grid)
% Set cell size
dx = W / Nx; % x-direction
dy = H / Ny; % y-direction

% Set nodes
NxP = Nx + 2; % P nodes in x-direction
NyP = Ny + 2; % P nodes in y-direction
NxU = Nx + 1; % U nodes in x-direction
NyU = Ny + 2; % U nodes in y-direction
NxV = Nx + 2; % V nodes in x-direction
NyV = Ny + 1; % V nodes in y-direction

disp('Spatial discretisation completed successfully.');


%% Initialisation of Variables
% Variables
P = zeros(NxP, NyP);
U = zeros(NxU, NyU);
V = zeros(NxV, NyV);

% Correction Variables
dU = zeros(NxU, NyU);
dV = zeros(NxV, NyV);
dP = zeros(NxP, NyP);

disp('Initialisation completed successfully.');


%% SIMPLE Algortihm
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


%% Export Results
export_to_nodes(U, V, P, dx, dy, Nx, Ny, 'results.txt');


%% Functions
% U-Momentum
function [U, dU, uRes] = solve_UMomentum(U_old, V_old, P, dx, dy, Re, alpha_u, Nx, Ny)
    % solve_UMomentum Computes the intermediate u-velocity field (u_star)
    % This function assembles the convective and diffusive transport coefficients
    % and incorporates the pressure gradient source term before invoking the TDMA

    % Initialise coefficient matrices
    NxU = Nx + 1; NyU = Ny + 2;
    aE = zeros(NxU, NyU);
    aW = zeros(NxU, NyU);
    aN = zeros(NxU, NyU);
    aS = zeros(NxU, NyU);
    aP = ones(NxU, NyU);
    s = zeros(NxU, NyU);
    
    dU = zeros(NxU, NyU);
    
    % Loop over interior U-nodes
    for i = 2:(NxU-1)
        for j = 2:(NyU-1)
            % Diffusive terms & Convective fluxes
            de = dy / (dx * Re);
            ge = dy * (U_old(i, j) + U_old(i+1, j)) / 2.0;
            dw = dy / (dx * Re);
            gw = dy * (U_old(i-1, j) + U_old(i, j)) / 2.0;

           if j == 2
               ds = 2.0 * dx / (dy * Re);
               gs = 0.0;
           else
               ds = dx / (dy * Re);
               gs = dx * (V_old(i, j-1) + V_old(i+1, j-1)) / 2.0;
           end
           if j == NyU-1
               dn = 2.0 * dx / (dy * Re);
               gn = 0.0;
           else
               dn = dx / (dy * Re);
               gn = dx * (V_old(i, j) + V_old(i+1, j)) / 2.0;
           end     
            
            % Coefficients (Central Differencing Scheme CDS)
            aE(i, j) = de - ge/2.0;
            aW(i, j) = dw + gw/2.0;
            aN(i, j) = dn - gn/2.0;
            aS(i, j) = ds + gs/2.0;
            
            aP(i, j) = de + dw + dn + ds;
            
            % Source term: Pressure gradient
            s(i, j) = (P(i, j) - P(i+1, j)) * dy; 
        end
    end

    % Apply Boundary Conditions for U
    % WEST WALL
    aE(1, :) = 0.0;
    aW(1, :) = 0.0;
    aN(1, :) = 0.0;
    aS(1, :) = 0.0;
    aP(1, :) = 1.0;
    s(1, :) = 0.0;
    % EAST WALL
    aE(NxU, :) = 0.0;
    aW(NxU, :) = 0.0;
    aN(NxU, :) = 0.0;
    aS(NxU, :) = 0.0;
    aP(NxU, :) = 1.0;
    s(NxU, :) = 0.0;
    % NORTH WALL
    aE(:, NyU) = 0.0;
    aW(:, NyU) = 0.0;
    aN(:, NyU) = 0.0;
    aS(:, NyU) = 0.0;
    aP(:, NyU) = 1.0;
    s(:, NyU) = 1.0;
    % SOUTH WALL
    aE(:, 1) = 0.0;
    aW(:, 1) = 0.0;
    aN(:, 1) = 0.0;
    aS(:, 1) = 0.0;
    aP(:, 1) = 1.0;
    s(:, 1) = 0.0;

    % Residual Calculation
    uRes = 0.0;
    for i = 2:(NxU-1)
        for j = 2:(NyU-1)
            tempRes = abs(aP(i,j)*U_old(i,j) ...
                        - aE(i,j)*U_old(i+1,j) - aW(i,j)*U_old(i-1,j) ...
                        - aN(i,j)*U_old(i,j+1) - aS(i,j)*U_old(i,j-1) ...
                        - s(i,j));
            if tempRes > uRes
                uRes = tempRes;
            end
        end
    end

    % Apply Underrelaxation & Store dU for Pressure Correction
    for i = 2:(NxU-1)
        for j = 2:(NyU-1)
            aP(i,j) = aP(i,j) / alpha_u; 
            s(i,j) = s(i,j) + ...
                (1.0 - alpha_u) * aP(i,j) * U_old(i,j);
            
            dU(i,j) = dy / aP(i,j);
        end
    end

    % Solve the Linear System using TDMA
    U = tridag_i(aS, aW, aP, aE, aN, s, U_old, NxU, NyU);
    U = tridag_j(aW, aS, aP, aN, aE, s, U, NyU, NxU);
        
end

% V-Momentum
function [V, dV, vRes] = solve_VMomentum(U_old, V_old, P, dx, dy, Re, alpha_v, Nx, Ny)
    % Solve_VMomentum Computes the intermediate v-velocity field (v_star)
    % This function constructs the discretised momentum balance in the y-direction,
    % ensuring the coupling of advection and diffusion terms is maintained.

    % Initialise coefficient matrices
    NxV = Nx + 2; NyV = Ny + 1;
    aE = zeros(NxV, NyV);
    aW = zeros(NxV, NyV);
    aN = zeros(NxV, NyV);
    aS = zeros(NxV, NyV);
    aP = ones(NxV, NyV);
    s = zeros(NxV, NyV);
    
    dV = zeros(NxV, NyV);
    
    % Loop over interior V-nodes
    for i = 2:(NxV-1)
        for j = 2:(NyV-1)
            % Diffusive terms & Convective fluxes
            dn = dx / (dy * Re);
            gn = dx * (V_old(i, j) + V_old(i, j+1)) / 2.0;
            ds = dx / (dy * Re);
            gs = dx * (V_old(i, j-1) + V_old(i, j)) / 2.0;
            if i == 2
                dw = 2.0 * dy / (dx * Re);
                gw = 0.0;
            else
                dw = dy / (dx * Re);
                gw = dy * (U_old(i-1, j) + U_old(i-1, j+1)) / 2.0;
            end
            if i == NxV-1
                de = 2.0 * dy / (dx * Re);
                ge = 0.0;
            else
                de = dy / (dx * Re);
                ge = dy * (U_old(i, j) + U_old(i, j+1)) / 2.0;
            end
          
            % Coefficients (Central Differencing Scheme CDS)
            aE(i, j) = de - ge/2.0;
            aW(i, j) = dw + gw/2.0;
            aN(i, j) = dn - gn/2.0;
            aS(i, j) = ds + gs/2.0;
            
            aP(i, j) = de + dw + dn + ds;
            
            % Source term: Pressure gradient
            s(i, j) = (P(i, j) - P(i, j+1)) * dy; 
        end
    end

    % Apply Boundary Conditions for V
    % WEST WALL
    aE(1, :) = 0.0;
    aW(1, :) = 0.0;
    aN(1, :) = 0.0;
    aS(1, :) = 0.0;
    aP(1, :) = 1.0;
    s(1, :) = 0.0;
    % EAST WALL
    aE(NxV, :) = 0.0;
    aW(NxV, :) = 0.0;
    aN(NxV, :) = 0.0;
    aS(NxV, :) = 0.0;
    aP(NxV, :) = 1.0;
    s(NxV, :) = 0.0;
    % NORTH WALL
    aE(:, NyV) = 0.0;
    aW(:, NyV) = 0.0;
    aN(:, NyV) = 0.0;
    aS(:, NyV) = 0.0;
    aP(:, NyV) = 1.0;
    s(:, NyV) = 0.0;
    % SOUTH WALL
    aE(:, 1) = 0.0;
    aW(:, 1) = 0.0;
    aN(:, 1) = 0.0;
    aS(:, 1) = 0.0;
    aP(:, 1) = 1.0;
    s(:, 1) = 0.0;

    % Residual Calculation
    vRes = 0.0;
    for i = 2:(NxV-1)
        for j = 2:(NyV-1)
            tempRes = abs(aP(i,j)*V_old(i,j) ...
                        - aE(i,j)*V_old(i+1,j) - aW(i,j)*V_old(i-1,j) ...
                        - aN(i,j)*V_old(i,j+1) - aS(i,j)*V_old(i,j-1) ...
                        - s(i,j));
            if tempRes > vRes
                vRes = tempRes;
            end
        end
    end

    % Apply Underrelaxation & Store dU for Pressure Correction   
    for i = 2:(NxV-1)
        for j = 2:(NyV-1)
            aP(i,j) = aP(i,j) / alpha_v; 
            s(i,j) = s(i,j) + ...
                (1.0 - alpha_v) * aP(i,j) * V_old(i,j);
            
            dV(i,j) = dx / aP(i,j);
        end
    end

    % Solve the Linear System using TDMA
    V = tridag_i(aS, aW, aP, aE, aN, s, V_old, NxV, NyV);
    V = tridag_j(aW, aS, aP, aN, aE, s, V, NyV, NxV);

end

% P correction
function [P_prime, pRes] = solve_PCorrection(U, V, dU, dV, dx, dy, Nx, Ny)
    % solve_PCorrection Solves the Poisson equation for pressure correction
    % This function enforces mass conservation by translating the velocity 
    % divergence (continuity error) into a pressure correction field.

    % Initialize coefficient matrices
    NxP = Nx + 2; NyP = Ny + 2;
    aE = zeros(NxP, NyP);
    aW = zeros(NxP, NyP);
    aN = zeros(NxP, NyP);
    aS = zeros(NxP, NyP);
    aP = ones(NxP, NyP);
    s = zeros(NxP, NyP);
    
    % Loop over interior Pressure nodes
    for i = 2:(NxP-1)
        for j = 2:(NyP-1)
            % Coefficients based on the dU and dV saved earlier
            aE(i,j) = dy * dU(i, j);
            aW(i,j) = dy * dU(i-1, j);
            aN(i,j) = dx * dV(i, j);
            aS(i,j) = dx * dV(i, j-1);
            
            aP(i,j) = aW(i,j) + aE(i,j) + aS(i,j) + aN(i,j);
            
            % Source term: Mass imbalance (Continuity Error) 
            % Mass In - Mass Out
            s(i,j) = dy * (U(i-1, j) - U(i, j)) + ...
                      dx * (V(i, j-1) - V(i, j));
        end
    end

    % Apply Boundary Conditions for U
    % WEST WALL
    aE(1, :) = 1.0;
    aW(1, :) = 0.0;
    aN(1, :) = 0.0;
    aS(1, :) = 0.0;
    aP(1, :) = 1.0;
    s(1, :) = 0.0;
    % EAST WALL
    aE(NxP, :) = 0.0;
    aW(NxP, :) = 1.0;
    aN(NxP, :) = 0.0;
    aS(NxP, :) = 0.0;
    aP(NxP, :) = 1.0;
    s(NxP, :) = 0.0;
    % NORTH WALL
    aE(:, NyP) = 0.0;
    aW(:, NyP) = 0.0;
    aN(:, NyP) = 0.0;
    aS(:, NyP) = 1.0;
    aP(:, NyP) = 1.0;
    s(:, NyP) = 0.0;
    % SOUTH WALL
    aE(:, 1) = 0.0;
    aW(:, 1) = 0.0;
    aN(:, 1) = 1.0;
    aS(:, 1) = 0.0;
    aP(:, 1) = 1.0;
    s(:, 1) = 0.0;
    
    % Residual Calculation
    pRes = 0.0;
    for i = 2:(NxP-1)
        for j = 2:(NyP-1)
            tempRes = abs(s(i,j));
            if tempRes > pRes
                pRes = tempRes;
            end
        end
    end

    % Solve the Linear System using TDMA
    P_prime = zeros(NxP, NyP);
    for i=1:10
        P_prime = tridag_i(aS, aW, aP, aE, aN, s, P_prime, NxP, NyP);
        P_prime = tridag_j(aW, aS, aP, aN, aE, s, P_prime, NyP, NxP);
    end

end

% Fields correction
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

function export_to_nodes(U, V, P, dx, dy, Nx, Ny, filename)
    % export_to_nodes Exports converged data to an external text file
    % This function formats the spatial coordinates and respective scalar/vector 
    % fields into a structured ASCII file for subsequent analysis.

    % Open the file for writing
    fileID = fopen(filename, 'w');
    
    % Write a clean header
    % There are Nx+1 and Ny+1 physical nodes (corners) in the mesh
    fprintf(fileID, 'Nx_nodes= %d, Ny_nodes= %d\n', Nx+1, Ny+1);
    fprintf(fileID, 'X\tY\tU\tV\tP\n');
    
    % Loop over the physical nodes
    for j = 1:(Ny+1)
        for i = 1:(Nx+1)
            % Physical coordinates of the node (origin at 0,0)
            x_coord = (i - 1) * dx;
            y_coord = (j - 1) * dy;
            
            % Interpolation to the nodes (corners)
            % P: Average of the four surrounding cells
            p_node = 0.25 * (P(i, j) + P(i+1, j) + P(i, j+1) + P(i+1, j+1));
            
            % U: Vertical average
            u_node = 0.5 * (U(i, j) + U(i, j+1));
            
            % V: Horizontal average
            v_node = 0.5 * (V(i, j) + V(i+1, j));
            
            % Write data line
            fprintf(fileID, '%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n', ...
                    x_coord, y_coord, u_node, v_node, p_node);
        end
    end
    
    fclose(fileID);
    fprintf('Results successfully exported to: %s\n', filename);
end


% END