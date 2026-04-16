function [max_stress, numElems, totalDOFs, max_u, strain_energy, solve_time] = run_fem_solver(mesh_filename, show_plot)
    % run_fem_solver: Calculates metrics for a 2D plate with a hole.
    % Inputs: 
    %   mesh_filename - string name of the .msh file (can include folder path)
    %   show_plot     - boolean (true/false) to show the plots on screen
    % Outputs:
    %   max_stress    - Maximum Von Mises stress (MPa)
    %   numElems      - Total number of elements
    %   totalDOFs     - Total degrees of freedom
    %   max_u         - Maximum X-displacement (mm)
    %   strain_energy - Total global strain energy (N-mm or Joules)
    %   solve_time    - Time to solve the linear system (seconds)
    
    % 1. Load the Mesh
    [NodeCoords, Connectivity] = read_gmsh_v2(mesh_filename);
    
    % 2. Physics Initialization
    E = 210000;         
    nu = 0.3;           
    t = 1.0;            
    sigma0 = 100;       
    
    D = (E / (1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    
    numNodes = size(NodeCoords, 1);
    numElems = size(Connectivity, 1);
    totalDOFs = 2 * numNodes; 
    
    K = sparse(totalDOFs, totalDOFs); 
    F = zeros(totalDOFs, 1);          
    
    gaussPts = [-1/sqrt(3), 1/sqrt(3)];
    weights = [1, 1];
    
    % 3. Global Stiffness Matrix Assembly
    for e = 1:numElems
        elemNodes = Connectivity(e, :);
        x = NodeCoords(elemNodes, 1);
        y = NodeCoords(elemNodes, 2);
        
        sctr = [2*elemNodes(1)-1, 2*elemNodes(1), ...
                2*elemNodes(2)-1, 2*elemNodes(2), ...
                2*elemNodes(3)-1, 2*elemNodes(3), ...
                2*elemNodes(4)-1, 2*elemNodes(4)];
                
        Ke = zeros(8, 8);
        
        for i = 1:2
            for j = 1:2
                xi = gaussPts(i);
                eta = gaussPts(j);
                w_i = weights(i);
                w_j = weights(j);
                
                dN_dxi = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta)];
                dN_deta = 0.25 * [-(1-xi), -(1+xi),  (1+xi),  (1-xi)];
                
                J = [dN_dxi * x,  dN_dxi * y; 
                     dN_deta * x, dN_deta * y];
                 
                detJ = det(J);
                if detJ <= 0
                    error('Element %d has a negative or zero Jacobian!', e);
                end
                
                dN_dx_dy = J \ [dN_dxi; dN_deta]; 
                dN_dx = dN_dx_dy(1, :);
                dN_dy = dN_dx_dy(2, :);
                
                B = zeros(3, 8);
                B(1, 1:2:7) = dN_dx;
                B(2, 2:2:8) = dN_dy;
                B(3, 1:2:7) = dN_dy;
                B(3, 2:2:8) = dN_dx;
                
                Ke = Ke + (B' * D * B) * detJ * w_i * w_j * t;
            end
        end
        K(sctr, sctr) = K(sctr, sctr) + Ke;
    end
    
    % --- CRITICAL: Save unpenalized K for true energy calculation ---
    K_orig = K; 
    
    % 4. Apply Boundary Conditions & Loads
    tolerance = 1e-6;
    left_nodes = find(abs(NodeCoords(:,1) - 0) < tolerance);
    bottom_nodes = find(abs(NodeCoords(:,2) - 0) < tolerance);
    right_nodes = find(abs(NodeCoords(:,1) - max(NodeCoords(:,1))) < tolerance);
    
    fixed_DOFs = unique([2 * left_nodes - 1; 2 * bottom_nodes]);
    
    [~, sort_idx] = sort(NodeCoords(right_nodes, 2));
    sorted_right_nodes = right_nodes(sort_idx);
    
    for i = 1:(length(sorted_right_nodes) - 1)
        n1 = sorted_right_nodes(i);
        n2 = sorted_right_nodes(i+1);
        dy = NodeCoords(n2, 2) - NodeCoords(n1, 2);
        edge_force = sigma0 * dy * t;
        F(2*n1 - 1) = F(2*n1 - 1) + edge_force / 2;
        F(2*n2 - 1) = F(2*n2 - 1) + edge_force / 2;
    end
    
    % 5. Solve Linear System
    penalty_value = 1e15; 
    K(fixed_DOFs, fixed_DOFs) = K(fixed_DOFs, fixed_DOFs) + penalty_value * speye(length(fixed_DOFs));
    
    tic; % Start timer
    U = K \ F;
    solve_time = toc; % Stop timer
    
    % --- Calculate Global Metrics ---
    U_x = U(1:2:end);
    U_y = U(2:2:end);
    max_u = max(U_x); % Max displacement in X (mm)
    
    % Strain Energy = 1/2 * U^T * K_original * U
    strain_energy = 0.5 * dot(U, K_orig * U); 
    
    % 6. Post-Processing (Stress & Energy Recovery)
    sigma_vm_elements = zeros(numElems, 1);
    strain_energy_elements = zeros(numElems, 1); % NEW: Elemental strain energy density
    
    for e = 1:numElems
        elemNodes = Connectivity(e, :);
        x = NodeCoords(elemNodes, 1);
        y = NodeCoords(elemNodes, 2);
        
        sctr = [2*elemNodes(1)-1, 2*elemNodes(1), ...
                2*elemNodes(2)-1, 2*elemNodes(2), ...
                2*elemNodes(3)-1, 2*elemNodes(3), ...
                2*elemNodes(4)-1, 2*elemNodes(4)];
        d_e = U(sctr);
        
        xi = 0; eta = 0;
        dN_dxi = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta)];
        dN_deta = 0.25 * [-(1-xi), -(1+xi),  (1+xi),  (1-xi)];
        J = [dN_dxi * x,  dN_dxi * y; dN_deta * x, dN_deta * y];
        dN_dx_dy = J \ [dN_dxi; dN_deta];
        dN_dx = dN_dx_dy(1, :); dN_dy = dN_dx_dy(2, :);
        
        B = zeros(3, 8);
        B(1, 1:2:7) = dN_dx; B(2, 2:2:8) = dN_dy; B(3, 1:2:7) = dN_dy; B(3, 2:2:8) = dN_dx;
        
        strain = B * d_e;
        stress = D * strain;
        sig_x = stress(1); sig_y = stress(2); tau_xy = stress(3);
        
        % Von Mises Stress
        sigma_vm_elements(e) = sqrt(sig_x^2 + sig_y^2 - sig_x*sig_y + 3*tau_xy^2);
        
        % Elemental Strain Energy Density (1/2 * stress * strain)
        strain_energy_elements(e) = 0.5 * dot(stress, strain); 
    end
    
    max_stress = max(sigma_vm_elements);
    max_energy_density = max(strain_energy_elements); % Extract max density for plotting
    
    % === 7. Visualization & Automated Saving ===
    [filepath, name, ~] = fileparts(mesh_filename);
    if isempty(filepath)
        filepath = pwd; 
    end
    if show_plot
        fig_vis = 'on';
    else
        fig_vis = 'off'; 
    end
    
    % Calculate scaled deformed coordinates once for all plots
    max_disp = max(abs(U));
    scale_factor = (max(NodeCoords(:,1)) * 0.05) / max_disp; 
    DefNodeCoords = NodeCoords;
    DefNodeCoords(:,1) = NodeCoords(:,1) + scale_factor * U_x;
    DefNodeCoords(:,2) = NodeCoords(:,2) + scale_factor * U_y;

    % 7.1 Undeformed Mesh Plot
    fig_mesh = figure('Name', sprintf('Mesh - %s', name), 'NumberTitle', 'off', 'Visible', fig_vis, 'Position', [100, 100, 800, 500]);
    patch('Faces', Connectivity, 'Vertices', NodeCoords, 'FaceColor', 'c', 'EdgeColor', 'k');
    axis equal;
    title(sprintf('Undeformed Mesh: %s', name), 'Interpreter', 'none');
    xlabel('X (mm)'); ylabel('Y (mm)');
    mesh_save_path = fullfile(filepath, sprintf('%s_mesh.png', name));
    exportgraphics(fig_mesh, mesh_save_path, 'Resolution', 300);

    % 7.2 Von Mises Stress Contour Plot
    fig_stress = figure('Name', sprintf('Stress - %s', name), 'NumberTitle', 'off', 'Visible', fig_vis, 'Position', [150, 150, 800, 500]);
    patch('Faces', Connectivity, 'Vertices', DefNodeCoords, ...
          'FaceVertexCData', sigma_vm_elements, 'FaceColor', 'flat', 'EdgeColor', '#444444');
    colormap(fig_stress, 'jet'); c_stress = colorbar; c_stress.Label.String = 'Von Mises Stress (MPa)';
    axis equal; 
    title(sprintf('Von Mises Stress: %s | Max: %.2f MPa', name, max_stress), 'Interpreter', 'none');
    xlabel('X (mm)'); ylabel('Y (mm)');
    stress_save_path = fullfile(filepath, sprintf('%s_stress.png', name));
    exportgraphics(fig_stress, stress_save_path, 'Resolution', 300);

    % 7.3 X-Displacement Contour Plot (Nodal Data -> Interpolated color)
    fig_disp = figure('Name', sprintf('X-Displacement - %s', name), 'NumberTitle', 'off', 'Visible', fig_vis, 'Position', [200, 200, 800, 500]);
    patch('Faces', Connectivity, 'Vertices', DefNodeCoords, ...
          'FaceVertexCData', U_x, 'FaceColor', 'interp', 'EdgeColor', '#444444');
    colormap(fig_disp, 'jet'); c_disp = colorbar; c_disp.Label.String = 'X-Displacement (mm)';
    axis equal; 
    title(sprintf('X-Displacement: %s | Max: %.5f mm', name, max_u), 'Interpreter', 'none');
    xlabel('X (mm)'); ylabel('Y (mm)');
    disp_save_path = fullfile(filepath, sprintf('%s_dispX.png', name));
    exportgraphics(fig_disp, disp_save_path, 'Resolution', 300);

    % 7.4 Strain Energy Density Contour Plot
    fig_energy = figure('Name', sprintf('Strain Energy - %s', name), 'NumberTitle', 'off', 'Visible', fig_vis, 'Position', [250, 250, 800, 500]);
    patch('Faces', Connectivity, 'Vertices', DefNodeCoords, ...
          'FaceVertexCData', strain_energy_elements, 'FaceColor', 'flat', 'EdgeColor', '#444444');
    colormap(fig_energy, 'jet'); c_energy = colorbar; c_energy.Label.String = 'Strain Energy Density (N-mm/mm^3)';
    axis equal; 
    title(sprintf('Strain Energy Density: %s | Max Density: %.2f | Total Energy: %.2f N-mm', name, max_energy_density, strain_energy), 'Interpreter', 'none');
    xlabel('X (mm)'); ylabel('Y (mm)');
    energy_save_path = fullfile(filepath, sprintf('%s_energy.png', name));
    exportgraphics(fig_energy, energy_save_path, 'Resolution', 300);

    % Close figures to free up memory if they are meant to run in the background
    if ~show_plot
        close(fig_mesh);
        close(fig_stress);
        close(fig_disp);
        close(fig_energy);
    end

    % === 8. Full Plate Symmetry Reflection & Visualization ===
    disp('Generating full-plate symmetry reflections...');
    
    % 8.1 Replicate Nodes for all 4 quadrants
    % Q1: Original (Top-Right), Q2: Top-Left, Q3: Bottom-Left, Q4: Bottom-Right
    Nodes_Q1 = NodeCoords;
    Nodes_Q2 = [-NodeCoords(:,1), NodeCoords(:,2)];
    Nodes_Q3 = [-NodeCoords(:,1), -NodeCoords(:,2)];
    Nodes_Q4 = [NodeCoords(:,1), -NodeCoords(:,2)];
    AllNodes = [Nodes_Q1; Nodes_Q2; Nodes_Q3; Nodes_Q4];
    
    % 8.2 Replicate Connectivity
    % Shift the node IDs for each quadrant block
    Conn_Q1 = Connectivity;
    Conn_Q2 = Connectivity + numNodes;
    Conn_Q3 = Connectivity + 2 * numNodes;
    Conn_Q4 = Connectivity + 3 * numNodes;
    AllConn = [Conn_Q1; Conn_Q2; Conn_Q3; Conn_Q4];
    
    % 8.3 Replicate Displacements (Vector reflection)
    AllUx = [U_x; -U_x; -U_x; U_x];
    AllUy = [U_y; U_y; -U_y; -U_y];
    
    % 8.4 Replicate Scalar Fields
    AllStress = [sigma_vm_elements; sigma_vm_elements; sigma_vm_elements; sigma_vm_elements];
    AllEnergy = [strain_energy_elements; strain_energy_elements; strain_energy_elements; strain_energy_elements];
    
    % 8.5 Calculate Full Plate Deformed Coordinates
    DefAllNodes = AllNodes;
    DefAllNodes(:,1) = AllNodes(:,1) + scale_factor * AllUx;
    DefAllNodes(:,2) = AllNodes(:,2) + scale_factor * AllUy;
    
    % --- Generate and Save Full Plate Plots ---
    
    % Full Undeformed Mesh
    fig_full_mesh = figure('Name', sprintf('Full Mesh - %s', name), 'NumberTitle', 'off', 'Visible', fig_vis, 'Position', [100, 100, 900, 500]);
    patch('Faces', AllConn, 'Vertices', AllNodes, 'FaceColor', 'c', 'EdgeColor', 'k', 'LineWidth', 0.1);
    axis equal;
    title(sprintf('Full Undeformed Mesh: %s', name), 'Interpreter', 'none');
    xlabel('X (mm)'); ylabel('Y (mm)');
    full_mesh_save_path = fullfile(filepath, sprintf('%s_full_mesh.png', name));
    exportgraphics(fig_full_mesh, full_mesh_save_path, 'Resolution', 300);

    % Full Von Mises Stress
    fig_full_stress = figure('Name', sprintf('Full Stress - %s', name), 'NumberTitle', 'off', 'Visible', fig_vis, 'Position', [150, 150, 900, 500]);
    patch('Faces', AllConn, 'Vertices', DefAllNodes, ...
          'FaceVertexCData', AllStress, 'FaceColor', 'flat', 'EdgeColor', 'none'); % EdgeColor none for smoother look on full plate
    colormap(fig_full_stress, 'jet'); c_full_stress = colorbar; c_full_stress.Label.String = 'Von Mises Stress (MPa)';
    axis equal; 
    title(sprintf('Full Plate Von Mises Stress: %s | Max: %.2f MPa', name, max_stress), 'Interpreter', 'none');
    xlabel('X (mm)'); ylabel('Y (mm)');
    full_stress_save_path = fullfile(filepath, sprintf('%s_full_stress.png', name));
    exportgraphics(fig_full_stress, full_stress_save_path, 'Resolution', 300);

    % Full X-Displacement
    fig_full_disp = figure('Name', sprintf('Full X-Disp - %s', name), 'NumberTitle', 'off', 'Visible', fig_vis, 'Position', [200, 200, 900, 500]);
    patch('Faces', AllConn, 'Vertices', DefAllNodes, ...
          'FaceVertexCData', AllUx, 'FaceColor', 'interp', 'EdgeColor', 'none');
    colormap(fig_full_disp, 'jet'); c_full_disp = colorbar; c_full_disp.Label.String = 'X-Displacement (mm)';
    axis equal; 
    title(sprintf('Full Plate X-Displacement: %s | Max: %.5f mm', name, max_u), 'Interpreter', 'none');
    xlabel('X (mm)'); ylabel('Y (mm)');
    full_disp_save_path = fullfile(filepath, sprintf('%s_full_dispX.png', name));
    exportgraphics(fig_full_disp, full_disp_save_path, 'Resolution', 300);

    % Full Strain Energy Density
    fig_full_energy = figure('Name', sprintf('Full Energy - %s', name), 'NumberTitle', 'off', 'Visible', fig_vis, 'Position', [250, 250, 900, 500]);
    patch('Faces', AllConn, 'Vertices', DefAllNodes, ...
          'FaceVertexCData', AllEnergy, 'FaceColor', 'flat', 'EdgeColor', 'none');
    colormap(fig_full_energy, 'jet'); c_full_energy = colorbar; c_full_energy.Label.String = 'Strain Energy Density (N-mm/mm^3)';
    axis equal; 
    title(sprintf('Full Plate Strain Energy Density: %s | Max Density: %.2f | Total Energy: %.2f N-mm', name, max_energy_density, strain_energy), 'Interpreter', 'none');
    xlabel('X (mm)'); ylabel('Y (mm)');
    full_energy_save_path = fullfile(filepath, sprintf('%s_full_energy.png', name));
    exportgraphics(fig_full_energy, full_energy_save_path, 'Resolution', 300);

    if ~show_plot
        close(fig_full_mesh);
        close(fig_full_stress);
        close(fig_full_disp);
        close(fig_full_energy);
    end
end