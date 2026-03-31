% =========================================================================
% CFD Post-Processing and Visualisation Routine
% Extracts nodal data from text output and generates flow contours.
% =========================================================================
clc; clear; close all;

% Set default interpreter to LaTeX for publication-quality rendering
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% 1. Data Extraction
filename = 'results.txt';

% Open file to read the header for grid dimensions
fid = fopen(filename, 'r');
if fid == -1
    error('File %s could not be opened. Please ensure it exists.', filename);
end

% Read the first line and extract numbers
header_line = fgetl(fid);
tokens = regexp(header_line, '\d+', 'match');
Nx_nodes = str2double(tokens{1});
Ny_nodes = str2double(tokens{2});
fclose(fid);

fprintf('Successfully detected grid size: %d x %d nodes.\n', Nx_nodes, Ny_nodes);

% Read the numerical data block (skipping the 2 header lines)
data = readmatrix(filename, 'NumHeaderLines', 2);

% Reshape data arrays. MATLAB fills column-wise, whilst the output loop 
% iterated x first, then y. Therefore, reshaping requires a transpose.
X_mat = reshape(data(:, 1), [Nx_nodes, Ny_nodes])';
Y_mat = reshape(data(:, 2), [Nx_nodes, Ny_nodes])';
U_mat = reshape(data(:, 3), [Nx_nodes, Ny_nodes])';
V_mat = reshape(data(:, 4), [Nx_nodes, Ny_nodes])';
P_mat = reshape(data(:, 5), [Nx_nodes, Ny_nodes])';

%% 2. Colormap Definition (Blue - White - Red)
n_colors = 256; 
half_n   = n_colors / 2;

c_min = [0.0, 0.0, 1.0]; % Pure Blue 
c_mid = [1.0, 1.0, 1.0]; % Pure White 
c_max = [1.0, 0.0, 0.0]; % Pure Red 

R = [linspace(c_min(1), c_mid(1), half_n), linspace(c_mid(1), c_max(1), half_n)];
G = [linspace(c_min(2), c_mid(2), half_n), linspace(c_mid(2), c_max(2), half_n)];
B = [linspace(c_min(3), c_mid(3), half_n), linspace(c_mid(3), c_max(3), half_n)];

cmap_cfd = [R', G', B'];

%% 3. Centreline Extraction
% Determine the indices nearest to the geometric centre
mid_x_idx = round(Nx_nodes / 2);
mid_y_idx = round(Ny_nodes / 2);

U_centreline = U_mat(:, mid_x_idx);
Y_centreline = Y_mat(:, mid_x_idx);

V_centreline = V_mat(mid_y_idx, :);
X_centreline = X_mat(mid_y_idx, :);

%% 4. Graphical Visualisation

% Figure 1: U-Velocity Contour
figure('Position', [100, 100, 600, 500]);
contourf(X_mat, Y_mat, U_mat, 30, 'LineColor', 'none');
colormap(cmap_cfd);
cbar = colorbar;
ylabel(cbar, 'Horizontal Velocity, $u$');
xlabel('Domain Width, $x$');
ylabel('Domain Height, $y$');
title('Contours of Horizontal Velocity ($u$)');
axis equal tight;
box on;

% Figure 2: V-Velocity Contour
figure('Position', [750, 100, 600, 500]);
contourf(X_mat, Y_mat, V_mat, 30, 'LineColor', 'none');
colormap(cmap_cfd);
cbar = colorbar;
ylabel(cbar, 'Vertical Velocity, $v$');
xlabel('Domain Width, $x$');
ylabel('Domain Height, $y$');
title('Contours of Vertical Velocity ($v$)');
axis equal tight;
box on;

% Figure 3: Centreline Profiles
figure('Position', [400, 300, 800, 400]);

% Subplot A: U velocity along vertical centreline
subplot(1, 2, 1);
plot(U_centreline, Y_centreline, 'k-', 'LineWidth', 1.5);
xlabel('Horizontal Velocity, $u$');
ylabel('Vertical Coordinate, $y$');
title(sprintf('$u$-velocity along $x = L/2$'));
grid on;

% Subplot B: V velocity along horizontal centreline
subplot(1, 2, 2);
plot(X_centreline, V_centreline, 'k-', 'LineWidth', 1.5);
xlabel('Horizontal Coordinate, $x$');
ylabel('Vertical Velocity, $v$');
title(sprintf('$v$-velocity along $y = H/2$'));
grid on;

fprintf('Visualisation complete. Figures generated successfully.\n');