% =========================================================================
% Horizontal Momentum Solver: solve_UMomentum
% 
% This routine resolves the discrete horizontal (u) momentum equation to 
% formulate an intermediate velocity field. It meticulously assembles the 
% convective and diffusive transport coefficients using a Central 
% Differencing Scheme, incorporates the pressure gradient source terms, 
% applies requisite boundary conditions, and integrates under-relaxation 
% before solving the linearised system via the Tri-Diagonal Matrix 
% Algorithm (TDMA).
% =========================================================================

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