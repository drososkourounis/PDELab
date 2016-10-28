classdef Simulator < handle
    
    properties
        
        mesh
        Identity2D
        Identity3D
        Laplace2D
        Laplace3D
        Ig
        Jg
        Kg 
        iKg
        A
        b
        x
        x0
        omega
    end
    
    
    methods
        
        
        
        function obj = Simulator(N, dimens)
            
            obj.x0 = [0.5 0.5 0.1];

            % create mesh object
            obj.mesh       = Mesh;
            
            % and differential operator objects
            obj.Identity3D = Identity3DQ1;
            obj.Identity2D = Identity2DQ1;
            
            obj.Laplace3D  = Laplace3DQ1;
            obj.Laplace2D  = Laplace2DQ1;
            
            % generate the box-type mesh

            mesh_time_start = tic;
            
            obj.mesh.makeRectilinearGrid(N, dimens);
            % obj.mesh.neighborhood();
            obj.mesh.makeParameters();

            mesh_time = toc(mesh_time_start);
            
            fprintf('MESH TIME: %10.2f\n', mesh_time);
            
            
        end
        
        
        
        function init(obj)
            
            obj.A   = [];
            obj.b   = [];
            obj.x   = [];
            obj.Ig  = [];
            obj.Jg  = [];
            obj.Kg  = [];
            obj.iKg = [];

        end


        function assemble(obj, w)

            obj.init();

            obj.omega = w;
            obj.assembleV();
            obj.assembleS();

            n     = obj.mesh.numberOfPoints;
          
            %obj.A = sparse(obj.Ig, obj.Jg, obj.Kg, n, n);
            obj.A = sparse(obj.Ig, obj.Jg, obj.Kg + i*obj.iKg, n, n);
   

        end



        function assembleV(obj)
            
            Me   = obj.Identity3D.makeForBrickElement(obj.mesh.dimens);
            Le   = obj.Laplace3D.makeForBrickElement(obj.mesh.dimens);
            
            me   = Me(:);
            le   = Le(:);
                        
            n     = obj.mesh.numberOfPoints;
            nel   = obj.mesh.numberOfElements;
            nv    = obj.mesh.numberOfVertices;
            
            assembly_time_start = tic;


            nx = obj.mesh.N(1);
            ny = obj.mesh.N(2);
            nz = obj.mesh.N(3);

            % volume contributions will be nel*nv^2
            % surface contributions will be nx*ny*2 + ny*nz*2 + nx*nz*2*4^2
            vol_size  = nel*nv^2;
            surf_size = (nx*ny*2 + ny*nz*2 + nx*nz*2)*16;

            obj.Kg=zeros(1, vol_size + surf_size);
            obj.iKg=zeros(1, vol_size + surf_size);

            i=1:nv;
            j=ones(1,nv);

            % ig = [1..nv, 1..nv, ... 1..nv] nv times
            ig=repmat(i,1,nv);

            % jg = [1 1 ... 1, 2 2 ... 2, ... , nv nv ... nv] each number
            % repeated nv times
            jg=repmat(1:nv,nv,1);

            %  (:) produces a column and we need a row of indices
            jg=jg(:)';

            obj.b = zeros(n,1);
                         
            k=1:nv^2;
            
            for e = 1:nel
                
                I          = obj.mesh.ElementList(e,:);
                x_e        = obj.mesh.PointList(I,:);
                c_e        = obj.mesh.Parameters(1, e);
                f_e        = obj.makeSource(x_e);
                
                obj.Ig(k)  = I(ig);
                obj.Jg(k)  = I(jg);
                obj.Kg(k)  = c_e^2*le - obj.omega^2*me;
                
                k          = k + nv^2;
                
                obj.b(I)   = obj.b(I) + Me*f_e;
            end
            
            assembly_time = toc(assembly_time_start);
            fprintf('OPERATOR ASSEMBLY TIME: %10.2f\n', assembly_time);
            
            
        end
        
        
        
        
        function assembleS(obj)
            Me   = obj.Identity2D.makeForBrickElement(obj.mesh.dimens);
            
            me   = Me(:);
                        
            nv    = 4;

            assembly_time_start = tic;


            % ig = [1..nv, 1..nv, ... 1..nv] nv times
            ig=repmat(1:nv,1,nv);

            % jg = [1 1 ... 1, 2 2 ... 2, ... , nv nv ... nv] each number
            % repeated nv times
            jg=repmat(1:nv,nv,1);

            %  (:) produces a column and we need a row of indices
            jg=jg(:)';

            k      = 1:nv^2;
            offset = obj.mesh.numberOfElements*64; %length(obj.Kg);
            k      = k + offset;
            
            surface_elements = size(obj.mesh.BoundaryFacet, 1);
            fprintf('surface elements %d\n', surface_elements) ;

            Facets= [ 1 4 3 2; 5 6 7 8; 1 5 8 4; 2 3 7 6; 1 2 6 5; 3 4 8 7];


            for se = 1:surface_elements
                
                e          = obj.mesh.BoundaryFacet(se, 1);
                f          = obj.mesh.BoundaryFacet(se, 2);
                F          = obj.mesh.ElementList( e, Facets(f, :) );
                area       = obj.mesh.Areas(f);
                c_e        = obj.mesh.Parameters(1, e);
                
                obj.Ig(k)  = F(ig);
                obj.Jg(k)  = F(jg);

                obj.iKg(k) = -c_e*me*area;
                
                k          = k + nv^2;
            end
            
            
            assembly_time = toc(assembly_time_start);
            fprintf('BOUNDARY ASSEMBLY TIME: %10.2f\n', assembly_time);
            
            
        end
 



        
        
        function assembleFast(obj, w)
            
            assembly_time_start = tic;
            
            
            % first make sure that we clean everything related to the 
            % linear system to avoid unpleasant surprices.
            
            clear obj.A;
            clear obj.b;
            clear obj.x;

            % %%%%%%%%%%%%%%%% %
            % Assemble the LHS %
            % %%%%%%%%%%%%%%%% %
            obj.omega      = w;
            

            % Local element mass matrix 
            Me    = obj.Identity3D.makeForBrickElement(obj.mesh.dimens);

            % Local element diffusion matrix 
            Le    = obj.Laplace3D.makeForBrickElement(obj.mesh.dimens);
            
            % Element matrices unfolded as a column vector column by column
            me    = Me(:);
            le    = Le(:);
            
            n     = obj.mesh.numberOfPoints;
            nel   = obj.mesh.numberOfElements;
            nv    = obj.mesh.numberOfVertices;
            
            % We repeat this vectors as many times as the number of elements 
            M     = repmat(me, nel, 1);
            L     = repmat(le, nel, 1);
            
            i=1:nv;
            j=ones(1,nv);

            % ig = [1..nv, 1..nv, ... 1..nv] nv times
            ig=repmat(i,1,nv);

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
            cg=repmat(1:nel, nv^2, 1);
             
            Ig    = obj.mesh.ElementList(:,ig)';
            Jg    = obj.mesh.ElementList(:,jg)';
            C     = obj.mesh.Parameters(1, cg(:))'.^2;
            
            Kg    = C.*L - obj.omega^2 * M;
            

            
                
            obj.A = sparse(Ig(:), Jg(:), Kg, n, n);
            
            clear i; 
            clear j; 

            clear ig; 
            clear jg; 
            clear cg;

            clear C;
            

            

            assembly_time = toc(assembly_time_start);
            fprintf('ASSEMBLY LHS TIME: %10.2f\n', assembly_time);
            


            % %%%%%%%%%%%%%%%% %
            % Assemble the RHS %
            % %%%%%%%%%%%%%%%% %


            f = obj.makeSource(obj.mesh.PointList);
            I = obj.mesh.ElementList'; 
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

            [I, J, obj.b] = find(B);

            assembly_time = toc(assembly_time_start);
            fprintf('ASSEMBLY RHS TIME: %10.2f\n', assembly_time);
            
            
        end
        

        
        
        
        
        
        
        
        function solve(obj)
            
            
            solution_time_start = tic;
            
            obj.x = obj.A \ obj.b;
            
            solution_time = toc(solution_time_start);
            fprintf('SOLUTION TIME: %10.2f\n', solution_time);
            
            r     = obj.b - obj.A*obj.x;
            fprintf('|| A x - b ||_2 = %25.16e\n', norm(r));
            
                        
            obj.mesh.Parameters(2,:) = real(obj.x)';
            obj.mesh.Parameters(3,:) = imag(obj.x)';
            
            
        end
        
        
        
        
        
        function f = makeSource(obj, p)
            
          % a single shot centered at x0, f0(x) = n*e^(-10*n*||x-x0||^2)
  
          n   = size(p, 1);

          dx  = p-repmat(obj.x0, n, 1);
          x   = sum(dx.^2', 1)'; 
          f   = -10*exp(-10*x);
        end
        
        
        


        function writeSolutionToFile(obj, filename)
        
            writeMeshToVTKFile(filename,             ...
                               obj.mesh.ElementList, ...
                               obj.mesh.PointList,   ...
                               obj.mesh.Parameters,  ...
                               obj.mesh.type);

        end
        
    end
end
