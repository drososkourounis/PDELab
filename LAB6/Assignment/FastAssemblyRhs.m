function b = assembleRhs(b, Mesh)

    assembly_time_start = tic;

    X = Mesh.Points(:,1);
    Y = Mesh.Points(:,2);
    f = source(X, Y);
    I = Mesh.Elements'; 
    F = f(I);

    Me = makeLocalMass(1, Mesh);
    b = Me*F;
    I = I(:);
    J = ones(Mesh.NumberOfVertices*Mesh.NumberOfElements, 1);

    % here we create a sparse matrix of size n x 1, which is practically
    % a sparse vector. We do it this way because the scattered FE vector
    % contributed by each element will be assembled on a global vector
    % once it becomes a sparse matrix. Then using the method 'find' we 
    % extract back our vector held in obj.b

    B = sparse(I, J, b(:), Mesh.NumberOfNodes, 1);
    b=full(B);

    assembly_time = toc(assembly_time_start);
    fprintf('ASSEMBLY RHS TIME: %10.2f\n', assembly_time);
end
