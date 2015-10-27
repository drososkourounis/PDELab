function [A, M, L] = assembleLhs(A, M, L, Mesh)

    assembly_time_start = tic;

    % Local element mass matrix 
    Me= makeLocalMass(1, Mesh);

    % Local element diffusion matrix 
    Le= makeLocalLaplacian(1, Mesh);

    % Element matrices unfolded as a column vector column by column
    me    = Me(:);
    le    = Le(:);

    n     = Mesh.NumberOfNodes;
    nel   = Mesh.NumberOfElements;
    nv    = Mesh.NumberOfVertices;

    % We repeat this vectors as many times as the number of elements 
    Ms     = repmat(me, nel, 1);
    Ls     = repmat(le, nel, 1);

    i=1:nv;
    j=ones(1, nv);

    % ig = [1..nv, 1..nv, ... 1..nv] nv times
    ig=repmat(i, 1, nv);

    % jg = [1 1 ... 1, 2 2 ... 2, ... , nv nv ... nv] each number
    % repeated nv times
    jg=repmat(1:nv, nv, 1);

    %  (:) produces a column and we need a row of indices
    jg=jg(:)';


    % [1 1 ... 1, 2 2 ... 2, ... , nv nv ... nv ] with each number
    % repeated nv^2 times, these will be the indices to choose the
    % weights for each element matrix (diffusion coefficients), place
    % them in an array and multiply the arrays with the corresponding
    % matrix
    cg  = repmat(1:nel, nv^2, 1);
    
    Is  = Mesh.Elements(:, ig)';
    Js  = Mesh.Elements(:, jg)';

    Operators.L = sparse(Is(:), Js(:), Ls, n, n);
    Operators.M = sparse(Is(:), Js(:), Ms, n, n);

    assembly_time = toc(assembly_time_start);
    fprintf('ASSEMBLY LHS TIME: %10.2f\n', assembly_time);    

end



