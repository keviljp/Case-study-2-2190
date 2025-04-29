%% Task 2, Step 1
% Analytical eigenfunctions:
% psi_{n_x,n_y}(x,y) = (2/L)*sin(n_x*pi*x/L)*sin(n_y*pi*y/L)
% The probability density is |psi(x,y)|^2.

% Clear workspace and figures
clear; clc; close all;

%%

% Parameters
L = 1;            % Length of the square well
N = 200;          % Number of grid points in each direction for smooth plots
x = linspace(0, L, N);
y = linspace(0, L, N);
[X, Y] = meshgrid(x, y);

% Define the normalized eigenfunction as an anonymous function
psi = @(nx, ny, X, Y, L) (2/L) * sin(nx * pi * X / L) .* sin(ny * pi * Y / L);

% Define the quantum numbers for the lowest five eigenvalues
states = [1, 1;
          1, 2;
          2, 1;
          2, 2;
          1, 3;
          3, 1;
          2, 3;
          3, 2];

% Create a figure for the contour plots
figure('Name', 'Analytical Probability Densities |\psi(x,y)|^2', 'NumberTitle', 'off');
% Loop over the states and plot the probability densities
for k = 1:size(states, 1)/2
    nx = states(k, 1);
    ny = states(k, 2);
    
    % Compute the eigenfunction and its probability density
    psi_val = psi(nx, ny, X, Y, L);
    prob_density = abs(psi_val).^2;
    
    % Create subplot for each eigenstate
    subplot(2, 2, k);
    contour(X, Y, prob_density, 20); 
    set(gca, 'FontSize', 15);
    title(['\psi_{', num2str(nx), ',', num2str(ny), '}']);
    xlabel('x');
    ylabel('y');
    colorbar;
    axis equal;
end

figure('Name', 'Analytical Probability Densities |\psi(x,y)|^2', 'NumberTitle', 'off');
for k = 1:size(states, 1)/2
    nx = states(k+4, 1);
    ny = states(k+4, 2);
    
    % Compute the eigenfunction and its probability density
    psi_val = psi(nx, ny, X, Y, L);
    prob_density = abs(psi_val).^2;
    
    % Create subplot for each eigenstate
    subplot(2, 2, k);
    contour(X, Y, prob_density, 20); 
    set(gca, 'FontSize', 15);
    title(['\psi_{', num2str(nx), ',', num2str(ny), '}']);
    xlabel('x');
    ylabel('y');
    colorbar;
    axis equal;
end

