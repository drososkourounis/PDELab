classdef Mesh < handle
    
    % Documentation
    % here we read a mesh written in the Tetgen or Triangle input
    % format; the user needs to provide the full path to the prefix of
    % the mesh files (filename without .ele or .node)
    
    properties %(GetAccess='private', SetAccess='private')
        dimension
        Ne
        Nv
        Np
        dimens
        N
        Elements
        Points
        PointMarkers
        DirichletPoints
        Parameters
        type
    end
    


    methods
        
        
        function obj = Mesh()
            
            obj.init();
            
        end
        
        
        function init(obj)
            obj.dimension          = 0;
            obj.Ne                 = 0;
            obj.Nv                 = 0;
            obj.Np                 = 0;
            obj.N                  = [];
            obj.dimens             = [];
            obj.Elements           = [];
            obj.Points             = [];
            obj.PointMarkers       = [];
            obj.DirichletPoints    = [];
            obj.Parameters         = [];
            obj.type               = 0;
        end
        
        


        function makeParameters(obj)
            fprintf('number of points is: %d\n', obj.Np);
            
            obj.Parameters(1,:) = ones(1,obj.Np);
            
            for e=1:obj.Ne
                I   = obj.Elements(e,:);
                x   = obj.Points(I,:);
                x_c = sum(x)/8;
                x_c2= x_c(2);
                Z   = x_c(3);
                
                if (0.4+0.1*x_c2 >= Z)
                    obj.Parameters(1,e) = 2.0;
                elseif (0.8 - 0.2*x_c2 >= Z)
                    obj.Parameters(1,e) = 1.5;
                else
                    obj.Parameters(1,e) = 3.0;
                end
            end
           
        end
        
   
        
        
        
        function makeRectilinearGrid(obj, N, D)
            
            obj.init();
            
            % N = [nx; ny; nz]
            % D = [dx; dy; dz]
            
            Nx            = N(1);
            Ny            = N(2);
            Nz            = N(3);
            m             = N(4);

            nx            = Nx + 2*m;
            ny            = Ny + 2*m;
            nz            = Nz + 2*m;
            
            obj.Ne        = nx*ny*nz;
            obj.Np        = (nx+1)*(ny+1)*(nz+1);
            obj.Nv        = 8;
            obj.dimension = 3;
            obj.type      = 12;

            dx            = D(1);
            dy            = D(2);
            dz            = D(3);
            
            obj.N         = N;
            obj.dimens    = D;
            
            x             = 0:dx:nx*dx;
            y             = 0:dy:ny*dy;
            z             = 0:dz:nz*dz;
            [X, Y, Z]     = meshgrid(x, y, z); 
            obj.Points    = [ X(:), Y(:), Z(:) ];
            fprintf('mesh::makeRectilinearGrid constructing mesh::Points ...\n');
            obj.PointMarkers = zeros(obj.Np, 1);
            I = find(obj.Points(:,1) == 0 | obj.Points(:,1) == 1 | obj.Points(:,2) == 0 | obj.Points(:,3) == 0);
            obj.PointMarkers(I) = 1;
            obj.DirichletPoints = I;

            
            fprintf('mesh::makeRectilinearGrid constructing mesh::Elements ...\n');
            
            nel             = nx*ny*nz;
            obj.Parameters  = zeros(2,(nx+1)*(ny+1)*(nz+1));

            k=1:nz;
            j=1:ny;
            i=1:nx;

            [I, J, K] = meshgrid(i, j, k);
            I = I(:); J = J(:); K = K(:);
            
            v1               = (K-1)*(nx+1)*(ny+1) + (J-1)*(nx+1) + I;
            v2               = (K-1)*(nx+1)*(ny+1) + (J-1)*(nx+1) + I+1;
            v3               = (K-1)*(nx+1)*(ny+1) +  J   *(nx+1) + I+1;
            v4               = (K-1)*(nx+1)*(ny+1) +  J   *(nx+1) + I;

            obj.Elements     = [v1 v2 v3 v4  [v1 v2 v3 v4]+(nx+1)*(ny+1) ];

        end
        

    end % methods
    
end





