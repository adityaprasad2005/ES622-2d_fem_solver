% === analytical_validation.m ===
clear; clc; close all;

disp('Starting Analytical Validation & Convergence Study...');

% 1. Analytical Baseline Setup
sigma0 = 100;          % Applied far-field tension (MPa)
Kt_infinite = 3.0;     % Theoretical Stress Concentration Factor (Infinite Plate)
analytical_limit = Kt_infinite * sigma0; 

% 2. Define Meshes to Evaluate
mesh_folder = 'grid convergence'; 
mesh_files = {
    fullfile(mesh_folder, 'mesh_coarse'), ...
    fullfile(mesh_folder, 'mesh_medium'), ...
    fullfile(mesh_folder, 'mesh_fine'), ...
    fullfile(mesh_folder, 'mesh_ultra')
};

num_studies = length(mesh_files);
fem_max_stress = zeros(num_studies, 1);
element_counts = zeros(num_studies, 1);

% 3. Run the FEM Solver for Each Mesh
for i = 1:num_studies
    current_mesh = mesh_files{i};
    fprintf('\nExtracting data from %s...\n', current_mesh);
    
    try
        % Call your existing solver (show_plot = false)
        [max_stress, numElems, ~, ~, ~, ~] = run_fem_solver(current_mesh, false);
        
        % Store the outputs
        fem_max_stress(i) = max_stress;
        element_counts(i) = numElems;
        
        fprintf('  Elements: %d | FEM Peak Stress: %.2f MPa\n', numElems, max_stress);
    catch ME
        fprintf('  ERROR reading %s: %s\n', current_mesh, ME.message);
    end
end

% 4. Plotting the Validation Curve
disp('Generating Validation Plot...');
fig_val = figure('Name', 'Analytical Validation', 'NumberTitle', 'off', 'Position', [200, 200, 800, 500]);
hold on; grid on;

% Plot FEM Convergence
plot(element_counts, fem_max_stress, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', '#0072BD', 'Color', '#0072BD', 'DisplayName', 'Custom FEM Solver');

% Plot Analytical Baseline (Infinite Plate)
yline(analytical_limit, '--r', 'LineWidth', 2, 'DisplayName', 'Kirsch Analytical (Infinite Plate, K_t = 3)');

% Formatting
title('Convergence of Peak Von Mises Stress vs. Analytical Baseline');
xlabel('Number of Elements (Mesh Density)');
ylabel('Maximum Von Mises Stress (MPa)');
legend('Location', 'southeast', 'FontSize', 11);

% Annotate the final converged FEM value
final_elems = element_counts(end);
final_stress = fem_max_stress(end);
text(final_elems, final_stress, sprintf('  Converged: %.1f MPa', final_stress), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'Color', '#0072BD');

% Save the plot automatically
val_save_path = fullfile(mesh_folder, 'analytical_validation.png');
exportgraphics(fig_val, val_save_path, 'Resolution', 300);

disp('Validation complete! Plot saved to grid convergence folder.');