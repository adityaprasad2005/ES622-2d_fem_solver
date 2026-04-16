function [NodeCoords, Connectivity] = read_gmsh_v2(filename)
    % Opens and reads a Gmsh Version 2 ASCII (.msh) file
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    
    NodeCoords = [];
    Connectivity = [];
    
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        
        % --- Parse Nodes ---
        if strcmp(line, '$Nodes')
            numNodes = str2double(fgetl(fid));
            % Preallocate for speed
            NodeCoords = zeros(numNodes, 2); 
            
            for i = 1:numNodes
                data = sscanf(fgetl(fid), '%f');
                % data(1) is Node ID, data(2) is X, data(3) is Y, data(4) is Z
                NodeCoords(data(1), :) = [data(2), data(3)]; 
            end
        end
        
        % --- Parse Elements ---
        if strcmp(line, '$Elements')
            numElems = str2double(fgetl(fid));
            % We don't know exactly how many are Q4s vs lines/points yet, 
            % so we will dynamically grow this array (or preallocate and trim)
            tempConn = zeros(numElems, 4);
            quadCount = 0;
            
            for i = 1:numElems
                data = sscanf(fgetl(fid), '%f');
                elemType = data(2);
                
                % In Gmsh V2, Element Type 3 represents a 4-node quadrangle
                if elemType == 3
                    quadCount = quadCount + 1;
                    % The last 4 numbers in the row are the node IDs
                    tempConn(quadCount, :) = data(end-3:end);
                end
            end
            % Trim the array to only include the quad elements
            Connectivity = tempConn(1:quadCount, :);
        end
    end
    
    fclose(fid);
    fprintf('Successfully loaded %d nodes and %d Q4 elements.\n', size(NodeCoords,1), size(Connectivity,1));
end