classdef Laplace3DQ1 < handle
  
  % Documentation
  % this method reads the coordinates for the current element and calculates the
  % barycentric coordinates and the volume for the specific element
  
  properties
    c        % the coordinates of the trilinear reference element
    C        % the coefficients of the map from reference to current
    T        % C = T \ X
    D        % The local gradient for all cubature points 
    Points   % cubature points
    Weights  % cubature weights
  end
  
  methods
    
     


    
    function L = makeUsingCubature(obj, number, mesh)
      
      I = mesh.ElementList(number, :)';
      X = mesh.PointList(I, :);
      C = obj.T\X;
      

      %     [ d x       d y       d z   ]             [ d xi     d eta    d zeta  ]
      %     [ ----     -----    ------- ]             [ ----     -----    ------- ]
      %     [ d xi     d xi      d xi   ]             [ d x      d x       d x    ]
      %     [                           ]             [                           ]
      %     [ d x       d y       d z   ]         -1  [ d xi     d eta    d zeta  ]
      % G = [ ----     -----    ------- ]        G  = [ ----     -----    ------- ]
      %     [ d eta    d eta     d eta  ]             [ d y      d y       d y    ]
      %     [                           ]             [                           ]
      %     [ d x       d y       d z   ]             [ d xi     d eta    d zeta  ]
      %     [ ----     -----    ------- ]             [ ----     -----    ------- ]
      %     [ d zeta   d zeta    d zeta ]             [ d z      d z       d z    ]
 

      %     [  d   ]     [ d xi     d eta    d zeta  ]   [  d   ] 
      %     [ ---- ]     [ ----     -----    ------- ]   [ ---- ] 
      %     [  dx  ]     [ d x      d x       d x    ]   [ dxi  ] 
      %     [      ]     [                           ]   [      ] 
      %     [  d   ]     [ d xi     d eta    d zeta  ]   [  d   ] 
      % D = [ ---- ]  =  [ ----     -----    ------- ]   [ ---- ] = grad
      %     [  dy  ]     [ d y      d y       d y    ]   [ deta ] 
      %     [      ]     [                           ]   [      ] 
      %     [  d   ]     [ d xi     d eta    d zeta  ]   [  d   ] 
      %     [ ---- ]     [ ----     -----    ------- ]   [ ---- ] 
      %     [  dz  ]     [ d z      d z       d z    ]   [ dzeta] 


      
      %  |)
      %  |
      %  |  grad(u) * grad(v) dVe  =   
      %  |      x         x
      % (|
      %   Ve
      
      %                                       Np                                           
      %  |)                                  ----                                            
      %  |   -1           -1                 \      -1                -1                   
      %  |  G   grad(u) * G grad(v) d V  =    )     G (Xi) grad(u) * G (Xi) grad(v) Wi * Ve
      %  |          X           X      0     /                 Xi               Xi         
      % (|                                   ----                                          
      %   V                                  i=1                                           
      %    0
      
      
      L       = zeros(8,8);
      npoints = length(obj.Weights);
      
      fprintf('npoints= %d\n', npoints);
      for   k = 1:npoints
       
       xi     = obj.Points(k, 1);
       eta    = obj.Points(k, 2);
       zeta   = obj.Points(k, 3);

      %     [ d x       d y       d z   ]             [ d xi     d eta    d zeta  ]
      %     [ ----     -----    ------- ]             [ ----     -----    ------- ]
      %     [ d xi     d xi      d xi   ]             [ d x      d x       d x    ]
      %     [                           ]             [                           ]
      %     [ d x       d y       d z   ]         -1  [ d xi     d eta    d zeta  ]
      % G = [ ----     -----    ------- ]        G  = [ ----     -----    ------- ]
      %     [ d eta    d eta     d eta  ]             [ d y      d y       d y    ]
      %     [                           ]             [                           ]
      %     [ d x       d y       d z   ]             [ d xi     d eta    d zeta  ]
      %     [ ----     -----    ------- ]             [ ----     -----    ------- ]
      %     [ d zeta   d zeta    d zeta ]             [ d z      d z       d z    ]
 

       % x = a11 + a12*xi + a13*eta + a14*zeta + a15*xi*eta +
       %     a16*eta*zeta + a17*xi*zeta + a18*xi*eta*zeta
   
       % y = a21 + a22*xi + a23*eta + a24*zeta + a25*xi*eta +
       %     a26*eta*zeta + a27*xi*zeta + a28*xi*eta*zeta

       % z = a31 + a32*xi + a33*eta + a34*zeta + a35*xi*eta +
       %     a36*eta*zeta + a37*xi*zeta + a38*xi*eta*zeta

       % d[x y z]/dxi
       G(1,:) = C(2,:) + C(5,:)*eta + C(7,:)*zeta + C(8,:)*eta*zeta;
       
       % d[x y z]/deta
       G(2,:) = C(3,:) + C(5,:)*xi  + C(6,:)*zeta + C(8,:)*xi*zeta;
       
       % d[x y z]/dzeta
       G(3,:) = C(4,:) + C(6,:)*eta + C(7,:)*xi   + C(8,:)*xi*eta;
       
       % =======================================
       % Ni(X,Y,Z) = (1 + a1*X)(1+a2*Y)(1+a3*Z),
       % ---------------------------------------
       % a1=c(i,1), a2=c(i,2), a3=(i,3)
       % =======================================
       

      %                                       Np                                           
      %  |)                                  ----                                            
      %  |   -1           -1                 \      -1                -1                   
      %  |  G   grad(u) * G grad(v) d V  =    )     G (Xk) grad(u) * G (Xk) grad(v) Wk * Ve
      %  |          X           X      0     /                 Xk               Xk         
      % (|                                   ----                                          
      %   V                                  k=1                                           
      %    0
 

       rows   = 8*(k-1) + 1 : 8*k;
       grad_k = obj.D(rows, :)';

       fprintf('det(G) = %23.16f\n', det(G));
       size(grad_k)
       GD     = inv(G) * grad_k;
       L      = L + obj.Weights(k)*det(G)*GD'*GD;
  
      end
    end

    

    
    function L = makeUsingSymbolic(obj)
      
      % This is the implementation of the exact symbolic integration on 
      % the reference element [0 1]^3. It is included for comprehension
      % and debugging purposes

     x = sym('x', 'real'); 
     y = sym('y', 'real'); 
     z = sym('z', 'real'); 

     dx = sym('dx', 'real'); 
     dy = sym('dy', 'real'); 
     dz = sym('dz', 'real'); 

     c=[-1   -1    -1;  ...
         1   -1    -1;  ...
         1    1    -1;  ...
        -1    1    -1;  ...
        -1   -1     1;  ...
         1   -1     1;  ...
         1    1     1;  ...
        -1    1     1];



      %                                  Np                                           
      %  |)                             ----                                          
      %  |                              \                                             
      %  |  grad(u) * grad(v) dVe  =     )      grad(u) * grad(v) Wi * Ve
      %  |      x         x             /          Xi        Xi         
      % (|                              ----                                          
      %   Ve                            i=1                                           
      

      % x = a11 + a12*xi + a13*eta + a14*zeta + a15*xi*eta +
      %     a16*eta*zeta + a17*xi*zeta + a18*xi*eta*zeta
   
      % y = a21 + a22*xi + a23*eta + a24*zeta + a25*xi*eta +
      %     a26*eta*zeta + a27*xi*zeta + a28*xi*eta*zeta

      % z = a31 + a32*xi + a33*eta + a34*zeta + a35*xi*eta +
      %     a36*eta*zeta + a37*xi*zeta + a38*xi*eta*zeta

      % d[x y z]/dxi

      % G      = zeros(3,3);
      G(1,1) = 0.5*dx;
      
      % d[x y z]/deta
      G(2,2) = 0.5*dy;
      
      % d[x y z]/dzeta
      G(3,3) = 0.5*dz;
 
      %                                  
      %  |)                             
      %  |   -1           -1            
      %  |  G   grad(u) * G grad(v) d V 
      %  |          X           X      0
      % (|                              
      %   V                             
      %    0
 
      for i=1:8
        N(i) = 1/8*( 1 + obj.c(i,1)*x )*( 1+obj.c(i,2)*y )*( 1+obj.c(i,3)*z );
      end

      Nx  = diff(N, 'x');
      Ny  = diff(N, 'y');
      Nz  = diff(N, 'z');

      D2N = [Nx; Ny; Nz];

      GD  = inv(G) * D2N;
      R   = det(G)*GD'*GD;
      L   = int(int(int(R, 'x', -1, 1), 'y', -1, 1), 'z', -1, 1);

    end





  
    function L = makeForBrickElement(obj, dimens)

      % For brick elements of size [dx, dy, dz] the element local mass matrix is
      dx = dimens(1);
      dy = dimens(2);
      dz = dimens(3);

      
     L = [      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz);
          (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz);
         (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz);
          (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz);
          (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz);
         (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy);
       - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz);
         (dy*dz)/(36*dx) - (dx*(dy^2/18 + dz^2/18))/(dy*dz), - (dx*dz)/(36*dy) - (dy*(dx^2/36 + dz^2/36))/(dx*dz),   (dx*dz)/(36*dy) - (dy*(dx^2/18 + dz^2/18))/(dx*dz),    (dx*dz)/(18*dy) - (dy*(dx^2/9 - dz^2/18))/(dx*dz),    (dy*dz)/(18*dx) + (dx*(dy^2/18 - dz^2/9))/(dy*dz),   (dy*(dx^2/36 - dz^2/18))/(dx*dz) - (dx*dz)/(18*dy),    (dx*dz)/(18*dy) + (dy*(dx^2/18 - dz^2/9))/(dx*dz),      (dx*dz)/(9*dy) + (dy*(dx^2/9 + dz^2/9))/(dx*dz)];
       

      end



  end     % methods
end      % Laplace3DQ1

