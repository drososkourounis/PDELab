classdef Simulator < handle
    
    properties
        
        mesh
        Identity3D  % for local mass
        Laplace3D   % for local Laplace
        M           % the global mass matrix
        L           % the global Laplacian
        P           % the preconditioner (Cholesky)
        A           % the lhs (Laplace + boundary conditions)
        b           % the rhs
        x           % the solution vector

    end
    
    
    methods
        

        function obj = Simulator(N, dimens)
            
            obj.init();
            % create mesh object
            obj.mesh       = Mesh;
            
            % and differential operator objects
            obj.Identity3D = Identity3DQ1;
            obj.Laplace3D  = Laplace3DQ1;
            
            % generate the box-type mesh

            mesh_time_start = tic;
            
            obj.mesh.makeRectilinearGrid(N, dimens);

            mesh_time = toc(mesh_time_start);
            
            fprintf('MESH TIME: %10.2f\n', mesh_time);
            
            
        end
        
        
        
        function init(obj)
            
            obj.A   = [];
            obj.L   = [];
            obj.M   = [];
            obj.P   = [];
            obj.b   = [];
            obj.x   = [];

        end


        function assemble(obj)

            obj.assembleV();
            obj.applyBoundaryConditions();

        end



        function applyBoundaryConditions(obj)

            assembly_time_start = tic;
            I          = obj.mesh.DirichletPoints;
            u0         = obj.u0(obj.mesh.Points(I,:)); 
            A_I        = obj.A(:, I);
            obj.b      = obj.b - A_I*u0;
            obj.b(I)   = u0;

            obj.A(:,I) = 0; 
            obj.A(I,:) = 0;
            obj.A(I,I) = speye(length(I));


            assembly_time = toc(assembly_time_start);
            fprintf('BOUNDARY CONDITIONS TIME: %10.2f\n', assembly_time);
            fprintf('We have in total %d boundary points\n', length(I));
        end


        
        function assembleV(obj)
            
            assembly_time_start = tic;
            
            
            % first make sure that we clean everything related to the 
            % linear system to avoid unpleasant surprices.
            
            clear obj.A;
            clear obj.M;
            clear obj.L;
            clear obj.P;
            clear obj.b;
            clear obj.x;

            % %%%%%%%%%%%%%%%% %
            % Assemble the LHS %
            % %%%%%%%%%%%%%%%% %

            % Local element mass matrix 
            Me    = obj.Identity3D.makeForBrickElement(obj.mesh.dimens);

            % Local element diffusion matrix 
            Le    = obj.Laplace3D.makeForBrickElement(obj.mesh.dimens);
            
            % Element matrices unfolded as a column vector column by column
            me    = Me(:);
            le    = Le(:);
            
            n     = obj.mesh.Np;
            nel   = obj.mesh.Ne;
            nv    = obj.mesh.Nv;
            
            % We repeat this vectors as many times as the number of elements 
            M     = repmat(me, nel, 1);
            L     = repmat(le, nel, 1);
            
            % ig = [1..nv, 1..nv, ... 1..nv] nv times
            ig=repmat(1:nv,1,nv);

            % jg = [1 1 ... 1, 2 2 ... 2, ... , nv nv ... nv] each number
            % repeated nv times
            jg=repmat(1:nv,nv,1);

            %  (:) produces a column and we need a row of indices
            jg=jg(:)';


            % [1 1 ... 1, 2 2 ... 2, ... , nv nv ... nv ] with each number
            % repeated nv^2 times, these will be the indices to choose the
            % weights for each element matrix (diffusion coefficients), place
            % them in an array and multiply the arrays with the corresponding
            % matrix
             
            I    = obj.mesh.Elements(:,ig)';
            J    = obj.mesh.Elements(:,jg)';
            I    = I(:); J = J(:);
                
            obj.A = sparse(I, J, L, n, n);
            %obj.M = sparse(I, J, M, n, n);

            %K     = sparse(I, J, L+i*M, n, n);
            %obj.A = real(K);
            %obj.L = obj.A;
            %obj.M = imag(K);

            

            assembly_time = toc(assembly_time_start);
            fprintf('ASSEMBLY LHS TIME: %10.2f\n', assembly_time);
            

            % %%%%%%%%%%%%%%%% %
            % Assemble the RHS %
            % %%%%%%%%%%%%%%%% %


            f = obj.makeSource(obj.mesh.Points);
            I = obj.mesh.Elements'; 
            F = f(I);
            b = Me*F;
            I = I(:);
            J = ones(nv*nel, 1);

            % here we create a sparse matrix of size n x 1, which is practically
            % a sparse vector. We do it this way because the scattered FE vector
            % contributed by each element will be assembled on a global vector
            % once it becomes a sparse matrix. Then using the method 'find' we 
            % extract back our vector held in obj.b

            assembly_time_start = tic;
            B = sparse(I, J, b(:), n, 1);

            obj.b = full(B);

            assembly_time = toc(assembly_time_start);
            fprintf('ASSEMBLY RHS TIME: %10.2f\n', assembly_time);
            
            
            clear ig, jg;
            clear I, J;
            clear L, M;
        end
        

        
        
        function solve(obj)
            
            solution_time_start = tic;
            
            %obj.x = obj.A \ obj.b;

            obj.P         = ichol(obj.A);
            solution_time = toc(solution_time_start);
            fprintf('ichol TIME: %10.2f\n', solution_time);
            solution_time_start = tic;

            %[obj.x,flag,relres,iter,resvec] = gmres(obj.A,obj.b,100,1e-8,1,@(x)precondition(obj, x));
            [obj.x,flag,relres,iter,resvec] = pcg(obj.A,obj.b,1e-8,100,@(x)precondition(obj, x));

            fprintf('iter %4d  %25.16e\n', [1:length(resvec) ; resvec' ]); 
            solution_time = toc(solution_time_start);
            fprintf('SOLUTION TIME: %10.2f\n', solution_time);
            
            r     = obj.b - obj.A*obj.x;
            fprintf('|| A x - b ||_2 = %25.16e\n', norm(r));
                        
            obj.mesh.Parameters(3,:) = obj.x;
        end
        
        
        function y = precondition(obj, x)
            y = obj.P'\(obj.P\x);
        end
        
        
        
        function f = makeSource(obj, p)
          
          x = p(:,1); 
          y = p(:,2); 
          z = p(:,3); 

          exponential = -x .* exp(-( y - 1 ).^2 .* ( z - 1 ).^2 );
          f = exponential.*( 4*(y - 1).^2 .* (z - 1).^4 - 2*(z - 1).^2 - 2*(y - 1).^2 + 4*(y - 1).^4 .* (z - 1).^2);
        end
        
        
         function u = u0(obj, p)
          
          x = p(:,1); 
          y = p(:,2); 
          z = p(:,3); 
          u = x .* exp(-( y - 1 ).^2 .* ( z - 1 ).^2 );
          
        end
        




       function errors(obj)

           U      = obj.u0(obj.mesh.Points);
           U_Uh   = U - obj.x;
           obj.mesh.Parameters(4,:) = abs(U_Uh);

           normL2 = U_Uh' * (obj.M * U_Uh);    

           normH1 = normL2 + U_Uh' * (obj.L * U_Uh);    

           fprintf('|| u - u_h ||_inf= %25.16e\n', max(abs(U_Uh)));
           fprintf('|| u - u_h ||_L2 = %25.16e\n', sqrt(normL2));
           fprintf('|| u - u_h ||_H1 = %25.16e\n', sqrt(normH1));

        end


        function writeSolutionToFile(obj, filename)
        
            writeMeshToVTKFile(filename,             ...
                               obj.mesh.Elements, ...
                               obj.mesh.Points,   ...
                               obj.mesh.Parameters,  ...
                               obj.mesh.type);

        end

       
    end   % methods
end
