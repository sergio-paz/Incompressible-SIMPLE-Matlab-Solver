% =========================================================================
% Data Export: export_to_nodes
% 
% This function facilitates the exportation of converged computational data 
% to an external text file. It systematically maps and interpolates the 
% staggered velocity and pressure fields onto the physical grid nodes 
% (corners), subsequently formatting the spatial coordinates alongside 
% their respective scalar and vector variables into a structured ASCII 
% format suitable for post-processing analysis.
% =========================================================================

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