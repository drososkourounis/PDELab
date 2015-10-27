function main(L_x, L_y, Ns)
   number_of_discretizations = length(Ns);
   fprintf('%5s  %15s %15s %12s %12s\n', 'N', '||u-u_h||_L2', '||u-u_h||_H1', 'assembly (s)', 'solve (s)');
   fprintf('%5s  %15s %15s %12s %12s\n', '-', '------------', '------------', '------------', '---------');
   
  
   for n=1:number_of_discretizations
       N_x = Ns(n);
       N_y = Ns(n);
       
       [L2, H1] = solvePDE(L_x, L_y, N_x, N_y);
       Y(n,:) = [L2, H1];
       H(n,1) = max(L_x / N_x, L_y / N_y);
   end

   loglog(H,Y);
   
end






function [L2,H1] = solvePDE(L_x, L_y, N_x, N_y)

% step 1: make grid
 [mesh] = makeGrid(L_x, L_y, N_x, N_y, 'triangles');
 
 assembly_time = tic;
% step 2: assemble discrete operators (matrices)
[M,K,b] = assembleDiscreteOperators(mesh);
 assembly_time = toc(assembly_time);

%step 2+1/2: enforce Dirichlet
[K, b] = EnforceDirichletBC(K, b);

solve_time = tic;
% step 3: solve linear system
u_h = solve(K,b);
solve_time = toc(solve_time);

% step 4: compute L^2 and H^1 norms
[L2, H1] = conv(u_h, mesh, M, K);

% the values of L^2 and H^1 norms are respectively
fprintf('%5d %15.6e %15.6e  %12.2f  %12.2f\n', N_x, L2, H1, assembly_time, solve_time);


% step 5: print vtk solution file

end





function [mesh] = makeGrid(L_x, L_y, N_x, N_y, grid_type)

if grid_type == 'triangles'
    N_e = 2 * N_x * N_y;
elseif grid_type == 'quadrilaterals'
    N_e = N_x * N_y;
end
    

N   = (N_x+1)*(N_y+1); %number of points 
N_v = 3;               %number of vertices

delta_x  = L_x / N_x;
delta_y  = L_y / N_y;

Points       = zeros(N, 2);
Elements     = zeros(N_e, N_v);
PointMarkers = zeros(N,1);

for j=1:(N_y+1)
    for i=1:(N_x+1)
        Pointmarkers(N,1)=1;
        for j=2:(N_y-1)
            for i=2:(N_x-1)
                PointMarkers(N,1)=0;
            end
        end
    end
end


for j=1:(N_y+1)
    for i=1:(N_x+1)
        x = (i-1)*delta_x;
        y = (j-1)*delta_y;
        n = i + (N_x+1)*(j-1);
        Points(n,1)=x;
        Points(n,2)=y;
    end
end


        

for j=1:N_y
    for i=1:N_x
        e   = i + (j-1)*N_x;
        n_1 = i     + (j-1)*(N_x+1);
        n_2 = (i+1) + (j-1)*(N_x+1);
        n_3 = (i+1) + j*(N_x+1);
        n_4 = i     + j*(N_x+1);
        if     grid_type == 'triangles'
            Elements(2*e-1,:)=[n_1, n_2, n_4];
            Elements(2*e,  :)=[n_2, n_3, n_4];
        elseif grid_type == 'quadrilaterals'
            Elements(e,   :)=[n_1, n_2, n_3, n_4];
        end
    end
end

mesh = struct('grid_type', 'triangles', 'L_x', L_x, 'L_y', L_y, ...
    'N_x', N_x, 'N_y', N_y, 'N_e', N_e, 'N', N, 'N_v', N_v, ...
    'delta_x', delta_x, 'delta_y', delta_y, ...
    'Points', Points, 'Elements', Elements, 'PointMarkers', PointMarkers);


%writeMeshToVTKFile(mesh, 'mesh');

end % function makeGrid(L_x, L_y, N_x, N_y, grid_type)




 
function [M,K,b] = assembleDiscreteOperators(mesh)
  N   = size(mesh.Points, 1);
  Ne  = size(mesh.Elements,1);
  Nv  = size(mesh.Elements,2);
  M   = sparse(N,N); K = sparse(N,N);
  b   = zeros(N,1);

  for e=1:Ne
    Me      = makeMe(e, mesh);
    Ke      = makeKe(e, mesh);
    fe      = makefe(e, mesh);
    I       = mesh.Elements(e, :);
    M(I, I) = M(I, I) + Me;
    K(I, I) = K(I, I) + Ke;
    b(I)    = b(I) + fe;
  end % e loop
end





function Me = makeMe(e, mesh)


V = (mesh.delta_x*mesh.delta_y)/2;

Me = (V/12) * [2 1 1; ...
               1 2 1; ...
               1 1 2];


end





function Ke = makeKe(e, mesh)

n      = mesh.Elements(e,:);
P      = [ mesh.Points(n, :) [1 1 1]'];

Ve     = (1/2)  * abs(det(P));

P      = inv(P);
P      = P(1:2,:);

Ke     = P' * P *  Ve;

end





function z = f(x, y)

    z = (-2*x.^2).*(2*x - 3).*(2*y - 3) - (2*y.^2).*(2*x - 3).*(2*y - 3) ...
      - (8*x.^2.*y).*(2*x - 3) - (8*x.*y.^2).*(2*y - 3);
end



function fe = makefe(e, mesh)

  n1 = mesh.Elements(e,1);
  n2 = mesh.Elements(e,2);
  n3 = mesh.Elements(e,3);
  
  x1 = mesh.Points(n1,1);
  y1 = mesh.Points(n1,2);
  x2 = mesh.Points(n2,1);
  y2 = mesh.Points(n2,2);
  x3 = mesh.Points(n3,1);
  y3 = mesh.Points(n3,2);

  z11 = 7.9742698535308731e-01;
  z12 = 1.0128650732345633e-01;
  z21 = 4.7014206410511505e-01;
  z22 = 5.9715871789769892e-02;
  z3  = 3.3333333333333331e-01;
  
  V = 0.5*abs(det([x1 y1 1; x2 y2 1; x3 y3 1]));
  T = inv([x1 y1 1; x2 y2 1; x3 y3 1]);
  
  c(1, 1) = z11*x1 + z12*x2 + z12*x3;
  c(1, 2) = z11*y1 + z12*y2 + z12*y3;
  c(1, 3) = 1.0;
  
  c(2, 1) = z12*x1 + z11*x2 + z12*x3;
  c(2, 2) = z12*y1 + z11*y2 + z12*y3;
  c(2, 3) = 1.0;
  
  c(3, 1) = z12*x1 + z12*x2 + z11*x3;
  c(3, 2) = z12*y1 + z12*y2 + z11*y3;
  c(3, 3) = 1.0;
  
  c(4, 1) = z22*x1 + z21*x2 + z21*x3;
  c(4, 2) = z22*y1 + z21*y2 + z21*y3;
  c(4, 3) = 1.0;
  
  c(5, 1) = z21*x1 + z22*x2 + z21*x3;
  c(5, 2) = z21*y1 + z22*y2 + z21*y3;
  c(5, 3) = 1.0;
  
  c(6, 1) = z21*x1 + z21*x2 + z22*x3;
  c(6, 2) = z21*y1 + z21*y2 + z22*y3;
  c(6, 3) = 1.0;

  c(7, 1) = z3*(x1 + x2 + x3);
  c(7, 2) = z3*(y1 + y2 + y3);    
  c(7, 3) = 1.0;
  
  weights(1,1) = 1.2593918054482717e-01;
  weights(2,1) = 1.2593918054482717e-01;
  weights(3,1) = 1.2593918054482717e-01;
  
  weights(4,1) = 1.3239415278850616e-01;
  weights(5,1) = 1.3239415278850616e-01;
  weights(6,1) = 1.3239415278850616e-01;
  
  weights(7) = 2.25e-01;


  F = f(c(:,1), c(:,2));
  N = c*T;

  fe = N' * (F.*weights) * V;

end





function w = u_0(x,y)

w = (x.^2.*y.^2).*(2*x-3).*(2*y-3);
end




function [K, b] = EnforceDirichletBC(K, b)
K(1,:)  = 0;
K(1, 1) = 1;
b(1)    = 0;
end





function u_h = solve(K,b)
  u_h = K \ b;
end



function [normL2, normH1] = conv(u_h, mesh, M, K)


x = mesh.Points(:,1);
y = mesh.Points(:,2);

u = u_0(x,y);


normL2 = sqrt((u-u_h)' * M * (u-u_h));

normH1 = sqrt((u-u_h)' * M *(u - u_h) + (u - u_h)' * K * (u - u_h));

end



