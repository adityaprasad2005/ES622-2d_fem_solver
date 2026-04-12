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
    max_u = max(U_x); % Max displacement in X (mm)
    
    % Strain Energy = 1/2 * U^T * K_original * U
    strain_energy = 0.5 * dot(U, K_orig * U); 
    
    % 6. Post-Processing (Stress Recovery)
    sigma_vm_elements = zeros(numElems, 1);
    
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
        
        sigma_vm_elements(e) = sqrt(sig_x^2 + sig_y^2 - sig_x*sig_y + 3*tau_xy^2);
    end
    
    max_stress = max(sigma_vm_elements);
    
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

    % Undeformed Mesh Plot
    fig_mesh = figure('Name', sprintf('Mesh - %s', name), 'NumberTitle', 'off', 'Visible', fig_vis, 'Position', [100, 100, 800, 500]);
    patch('Faces', Connectivity, 'Vertices', NodeCoords, 'FaceColor', 'c', 'EdgeColor', 'k');
    axis equal;
    title(sprintf('Undeformed Mesh: %s', name), 'Interpreter', 'none');
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    mesh_save_path = fullfile(filepath, sprintf('%s_mesh.png', name));
    exportgraphics(fig_mesh, mesh_save_path, 'Resolution', 300);

    % Deformed Contour Plot
    U_y = U(2:2:end);
    max_disp = max(abs(U));
    scale_factor = (max(NodeCoords(:,1)) * 0.05) / max_disp; 
    
    DefNodeCoords = NodeCoords;
    DefNodeCoords(:,1) = NodeCoords(:,1) + scale_factor * U_x;
    DefNodeCoords(:,2) = NodeCoords(:,2) + scale_factor * U_y;
    
    fig_contour = figure('Name', sprintf('Stress - %s', name), 'NumberTitle', 'off', 'Visible', fig_vis, 'Position', [150, 150, 800, 500]);
    patch('Faces', Connectivity, 'Vertices', DefNodeCoords, ...
          'FaceVertexCData', sigma_vm_elements, 'FaceColor', 'flat', 'EdgeColor', '#444444');
    colormap('jet'); c = colorbar; c.Label.String = 'Von Mises Stress (MPa)';
    axis equal; 
    title(sprintf('Mesh: %s | Max Stress: %.2f MPa', name, max_stress), 'Interpreter', 'none');
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    contour_save_path = fullfile(filepath, sprintf('%s_contour.png', name));
    exportgraphics(fig_contour, contour_save_path, 'Resolution', 300);

    if ~show_plot
        close(fig_mesh);
        close(fig_contour);
    end
end