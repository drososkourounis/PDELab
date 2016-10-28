type Mesh
  Ne::Int
  Np::Int
  Nv::Int
  Points
  Elements
  PointMarkers



    function Mesh(nx::Int, ny::Int, nz::Int)
        
        this = new();
        print("- Creating mesh ", nx, "x", ny, "x", nz, "\n");
    
        tic()
        dx = 1.0 / nx;
        dy = 1.0 / ny;
        dz = 1.0 / nz;
    
        this.Nv = 8;
        this.Ne = nx*ny*nz;
        this.Np = (nx+1)*(ny+1)*(nz+1);
    
        this.Points       = Array{Float64}(this.Np, 3);   
        this.Elements     = Array{Int}(this.Ne, this.Nv);   
        this.PointMarkers = Array{Int}(this.Np,1);
        for k in 0:nz
          z = k*dz;   
          for j in 0:ny
            y = j*dy;   
            for i in 0:nx
              x = i*dx;
              p = k*(nx+1)*(ny+1) + j*(nx+1) + i+1;
              this.Points[p, :] = [ x, y, z ];
              if  (i==1) || (i==nx+1) || (j==1) || (k==1) 
                this.PointMarkers[p] = 1;
              else
                this.PointMarkers[p] = 0;
              end
            end
          end
        end
    
        for k in 1:nz
          for j in 1:ny
            for i in 1:nx
               e                = (k-1)*nx*ny + (j-1)*nx + i;
               v1               = (k-1)*(nx+1)*(ny+1) + (j-1)*(nx+1) + i;
               v2               = (k-1)*(nx+1)*(ny+1) + (j-1)*(nx+1) + i+1;
               v3               = (k-1)*(nx+1)*(ny+1) +  j   *(nx+1) + i+1;
               v4               = (k-1)*(nx+1)*(ny+1) +  j   *(nx+1) + i;
               v5               =  k   *(nx+1)*(ny+1) + (j-1)*(nx+1) + i;
               v6               =  k   *(nx+1)*(ny+1) + (j-1)*(nx+1) + i+1;
               v7               =  k   *(nx+1)*(ny+1) +  j   *(nx+1) + i+1;
               v8               =  k   *(nx+1)*(ny+1) +  j   *(nx+1) + i;
               this.Elements[e, :] = [ v1,v2,v3,v4,v5,v6,v7,v8 ];
            end
          end
        end
              
        toc()
        return this; 
    
    end

end
