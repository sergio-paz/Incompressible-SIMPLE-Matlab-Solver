% =========================================================================
% Pressure Correction Solver: solve_PCorrection
% 
% This function resolves the Poisson equation for pressure correction, an 
% essential step in enforcing the strict conservation of mass. By 
% calculating the local velocity divergence (continuity error) from the 
% newly computed intermediate momentum fields, it establishes a source term 
% to determine the requisite pressure adjustments. The resulting field is 
% subsequently approximated across multiple TDMA iterations.
% =========================================================================

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