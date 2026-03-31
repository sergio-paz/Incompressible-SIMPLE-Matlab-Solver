% =========================================================================
% Vertical Momentum Solver: solve_VMomentum
% 
% This routine resolves the discrete vertical (v) momentum equation to 
% formulate an intermediate velocity field. It meticulously assembles the 
% convective and diffusive transport coefficients using a Central 
% Differencing Scheme, incorporates the pressure gradient source terms, 
% applies requisite boundary conditions, and integrates under-relaxation 
% before solving the linearised system via the Tri-Diagonal Matrix 
% Algorithm (TDMA).
% =========================================================================

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