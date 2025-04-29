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

%% Task 2, Step 2 & 3
% MATLAB Code to Solve the 2D Time-Independent Schrödinger Equation 

%% Step 2 Numerical solution
% Numerical Solution for 2D Infinite Square Well using Finite Difference Method
% Assumption: Lx = Ly = L = 1
% Constants (set hbar = 1, 2m = 1 for simplicity)
L = 1;            % Length of the well
N = 150;           % Number of interior points in each dimension (excluding boundaries)
h = L / (N+1);    % Grid spacing

% Create 1D Laplacian matrix
e = ones(N,1);
D = spdiags([e -2*e e], -1:1, N, N);
I = speye(N);

% 2D Laplacian operator using Kronecker sum
H = -(1/h^2) * (kron(D, I) + kron(I, D));

% Solve for the lowest 5 eigenvalues and eigenvectors
num_eigenstates = 8;
[Psi, E_matrix] = eigs(H, num_eigenstates, 'smallestreal');
E_numerical = diag(E_matrix); % Eigenvalues

% Set up grid for plotting
x = linspace(0, L, N+2); % include boundary points
y = linspace(0, L, N+2);
[X, Y] = meshgrid(x, y);

% Plotting numerical probability densities
figure;
for k = 1:num_eigenstates-4
    % Reshape eigenvector to 2D function
    psi_vec = Psi(:,k);
    psi_2D = reshape(psi_vec, [N,N]);
    
    % Add boundary zeros (since wavefunction is zero at the walls)
    psi_full = zeros(N+2,N+2);
    psi_full(2:end-1,2:end-1) = psi_2D;
    
    % Probability density
    prob_density = abs(psi_full).^2;
    
    % Plot
    subplot(2,2,k);
    contour(X, Y, prob_density, 20);
    colorbar;
    axis equal;
    title(['Numerical |\psi_{', num2str(k), '}|^2']);
end

figure;
for k = 5:num_eigenstates
    % Reshape eigenvector to 2D function
    psi_vec = Psi(:,k);
    psi_2D = reshape(psi_vec, [N,N]);
    
    % Add boundary zeros (since wavefunction is zero at the walls)
    psi_full = zeros(N+2,N+2);
    psi_full(2:end-1,2:end-1) = psi_2D;
    
    % Probability density
    prob_density = abs(psi_full).^2;
    
    % Plot
    subplot(2,2,k-4);
    contour(X, Y, prob_density, 20);
    colorbar;
    axis equal;
    title(['Numerical |\psi_{', num2str(k), '}|^2']);
end

% Print eigenvalues
disp('Numerical Eigenvalues (in units where \hbar^2/2m = 1):');
disp(E_numerical);

%% Step 3 Comparison
% Solve eigenvalue problem using eigs
[Psi, E_matrix] = eigs(H, num_eigenstates, 'smallestreal');
E_numerical = diag(E_matrix);

% Reshape numerical eigenfunctions
psi_numerical = cell(num_eigenstates,1);
max_prob_numerical = zeros(num_eigenstates,1);
for k = 1:num_eigenstates
    psi_vec = Psi(:,k);
    % Normalize properly for quantum mechanics
    psi_vec = psi_vec / sqrt(h^2 * sum(abs(psi_vec).^2)); 

    psi_2D = reshape(psi_vec, [N,N]);
    psi_full = zeros(N+2,N+2);
    psi_full(2:end-1,2:end-1) = psi_2D;
    psi_numerical{k} = abs(psi_full).^2;
    max_prob_numerical(k) = max(psi_numerical{k}(:));
end


%% Analytical Solution
% Quantum numbers corresponding to lowest eight states
quantum_numbers = [1 1; 1 2; 2 1; 2 2; 3 1; 1 3; 2 3; 3 2];

psi_analytical = cell(num_eigenstates,1);
E_analytical = zeros(num_eigenstates,1);
max_prob_analytical = zeros(num_eigenstates,1);

for k = 1:num_eigenstates
    nx = quantum_numbers(k,1);
    ny = quantum_numbers(k,2);
    psi = (2/L)*sin(nx*pi*X/L).*sin(ny*pi*Y/L);
    prob_density = abs(psi).^2;
    psi_analytical{k} = prob_density;
    max_prob_analytical(k) = max(prob_density(:));
    E_analytical(k) = pi^2*(nx^2 + ny^2); % Since \hbar^2/2m = 1
end

%% Plotting Analytical vs Numerical side by side
figure;
for k = 1:2
    % Analytical
    subplot(2,2,2*k-1);
    contour(X, Y, psi_analytical{k}, 20);
    colorbar;
    axis equal;
    title(['Analytical |\psi_{', num2str(quantum_numbers(k,1)), ',', num2str(quantum_numbers(k,2)), '}(x,y)|^2']);
    
    % Numerical
    subplot(2,2,2*k);
    contour(X, Y, psi_numerical{k}, 20);
    colorbar;
    axis equal;
    title(['Numerical |\psi_{', num2str(quantum_numbers(k,1)), ',', num2str(quantum_numbers(k,2)), '}(x,y)|^2']);
end

figure;
for k = 3:4
    % Analytical
    subplot(2,2,2*k-5);
    contour(X, Y, psi_analytical{k}, 20);
    colorbar;
    axis equal;
    title(['Analytical |\psi_{', num2str(quantum_numbers(k,1)), ',', num2str(quantum_numbers(k,2)), '}(x,y)|^2']);
    
    % Numerical
    subplot(2,2,2*k-4);
    contour(X, Y, psi_numerical{k}, 20);
    colorbar;
    axis equal;
    title(['Numerical |\psi_{', num2str(quantum_numbers(k,1)), ',', num2str(quantum_numbers(k,2)), '}(x,y)|^2']);
end

figure;
for k = 5:6
    % Analytical
    subplot(2,2,2*k-9);
    contour(X, Y, psi_analytical{k}, 20);
    colorbar;
    axis equal;
    title(['Analytical |\psi_{', num2str(quantum_numbers(k,1)), ',', num2str(quantum_numbers(k,2)), '}(x,y)|^2']);
    
    % Numerical
    subplot(2,2,2*k-8);
    contour(X, Y, psi_numerical{k}, 20);
    colorbar;
    axis equal;
    title(['Numerical |\psi_{', num2str(quantum_numbers(k,1)), ',', num2str(quantum_numbers(k,2)), '}(x,y)|^2']);
end

figure;
for k = 7:8
    % Analytical
    subplot(2,2,2*k-13);
    contour(X, Y, psi_analytical{k}, 20);
    colorbar;
    axis equal;
    title(['Analytical |\psi_{', num2str(quantum_numbers(k,1)), ',', num2str(quantum_numbers(k,2)), '}(x,y)|^2']);
    
    % Numerical
    subplot(2,2,2*k-12);
    contour(X, Y, psi_numerical{k}, 20);
    colorbar;
    axis equal;
    title(['Numerical |\psi_{', num2str(quantum_numbers(k,1)), ',', num2str(quantum_numbers(k,2)), '}(x,y)|^2']);
end

%% Create comparison table
ComparisonTable = table((1:num_eigenstates)', quantum_numbers(:,1), quantum_numbers(:,2), ...
    E_analytical, E_numerical, max_prob_analytical, max_prob_numerical, ...
    'VariableNames', {'State', 'n_x', 'n_y', 'E_analytical', 'E_numerical', 'MaxProb_Analytical', 'MaxProb_Numerical'});

disp('Comparison Table:');
disp(ComparisonTable);
