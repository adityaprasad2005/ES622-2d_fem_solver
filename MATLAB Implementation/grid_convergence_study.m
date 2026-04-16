% === grid_convergence_study.m ===
clear; clc; close all;

disp('Starting Grid Convergence Study...');

% 1. Define the subfolder and the mesh files to test
mesh_folder = 'grid convergence'; 

mesh_files = {
    fullfile(mesh_folder, 'mesh_coarse'), ...
    fullfile(mesh_folder, 'mesh_medium'), ...
    fullfile(mesh_folder, 'mesh_fine'), ...
    fullfile(mesh_folder, 'mesh_ultra')
};

num_studies = length(mesh_files);

% Pre-allocate arrays to store results
results_max_stress = zeros(num_studies, 1);
results_elements = zeros(num_studies, 1);
results_dofs = zeros(num_studies, 1);
results_max_u = zeros(num_studies, 1);
results_energy = zeros(num_studies, 1);
results_time = zeros(num_studies, 1);

% 2. Loop through each mesh and run the solver
for i = 1:num_studies
    current_mesh = mesh_files{i};
    fprintf('\nAnalyzing %s...\n', current_mesh);
    
    try
        % Call the updated solver (false = do not pop up individual plots)
        [max_stress, numElems, totalDOFs, max_u, strain_energy, solve_time] = run_fem_solver(current_mesh, false);
        
        % Store the outputs
        results_max_stress(i) = max_stress;
        results_elements(i) = numElems;
        results_dofs(i) = totalDOFs;
        results_max_u(i) = max_u;
        results_energy(i) = strain_energy;
        results_time(i) = solve_time;
        
        fprintf('  Elements: %d | DOFs: %d\n', numElems, totalDOFs);
        fprintf('  Max Stress: %.2f MPa | Max Disp: %.5f mm\n', max_stress, max_u);
        fprintf('  Strain Energy: %.2f N-mm | Time: %.3f s\n', strain_energy, solve_time);
        
    catch ME
        fprintf('  ERROR reading %s: %s\n', current_mesh, ME.message);
    end
end

% 3. Plot the Comprehensive Convergence Curves
disp('Generating Convergence Curves...');
fig_dashboard = figure('Name', 'Comprehensive Grid Convergence', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 800]);

% Subplot 1: Maximum Stress (Local Convergence)
subplot(2, 2, 1);
plot(results_elements, results_max_stress, '-o', 'LineWidth', 2, 'Color', '#D95319', 'MarkerFaceColor', '#D95319');
title('Max Von Mises Stress (Local)');
xlabel('Number of Elements'); ylabel('Stress (MPa)'); grid on;

% Subplot 2: Maximum Displacement (Kinematic Convergence)
subplot(2, 2, 2);
plot(results_elements, results_max_u, '-s', 'LineWidth', 2, 'Color', '#0072BD', 'MarkerFaceColor', '#0072BD');
title('Max X-Displacement');
xlabel('Number of Elements'); ylabel('Displacement (mm)'); grid on;

% Subplot 3: Total Strain Energy (Global Convergence)
subplot(2, 2, 3);
plot(results_elements, results_energy, '-^', 'LineWidth', 2, 'Color', '#77AC30', 'MarkerFaceColor', '#77AC30');
title('Total Strain Energy (Global)');
xlabel('Number of Elements'); ylabel('Energy (N-mm)'); grid on;

% Subplot 4: Computational Cost
subplot(2, 2, 4);
plot(results_elements, results_time, '-d', 'LineWidth', 2, 'Color', '#7E2F8E', 'MarkerFaceColor', '#7E2F8E');
title('Computational Cost');
xlabel('Number of Elements'); ylabel('Solve Time (Seconds)'); grid on;

% Save the dashboard figure automatically to the grid convergence folder
dashboard_path = fullfile(mesh_folder, 'convergence_dashboard.png');
exportgraphics(fig_dashboard, dashboard_path, 'Resolution', 300);

disp('Convergence Study Complete! Dashboard saved.');