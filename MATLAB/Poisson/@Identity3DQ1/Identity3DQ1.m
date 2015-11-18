classdef Identity3DQ1 < handle
    
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
        
               
        
        
        
        function M = makeUsingCubature(obj, number, Mesh)
            
            I = Mesh.ElementList(number, :)';
            X = Mesh.PointList(I, :);
            obj.C = obj.T \ X;
            
            %     [ d x       d x       d x   ]
            %     [ ----     -----    ------- ]
            %     [ d xi     d eta     d zeta ]
            %     [                           ]
            %     [ d y       d y       d y   ]
            % G = [ ----     -----    ------- ]
            %     [ d xi     d eta     d zeta ]
            %     [                           ]
            %     [ d z       d z       d z   ]
            %     [ ----     -----    ------- ]
            %     [ d xi     d eta     d zeta ]
            
            for   k=1:npoints
                
                xi     = obj.Points(k, 1);
                eta    = obj.Points(k, 2);
                zeta   = obj.Points(k, 3);
                
                G(1,:) = obj.C(2,:) + obj.C(5,:)*eta + obj.C(7,:)*zeta + obj.C(8,:)*eta*zeta;
                G(2,:) = obj.C(3,:) + obj.C(5,:)*xi  + obj.C(6,:)*zeta + obj.C(8,:)*xi*zeta;
                G(3,:) = obj.C(4,:) + obj.C(6,:)*eta + obj.C(7,:)*xi   + obj.C(8,:)*xi*eta;
                
                g(k,1) = det(G)*obj.Weights(k);
                
            end
            
            N = FE3DQ1(obj.Points);
            M = (g.*N)*N';
            
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
            %  |
            %  |  Ni * Nj * det(G) dV
            %  |                     0
            % (|
            %   V
            %    0
            
            
            for i=1:8
                N(i) = 1/8*( 1 + c(i,1)*x )*( 1 + c(i,2)*y )*( 1 + c(i,3)*z );
            end
            
            I   = det(G)*N'*N;
            L   = int(int(int(I, 'x', -1, 1), 'y', -1, 1), 'z', -1, 1);
            
        end
        
        
        
               
        function L = makeForBrickElement(obj, dimens)
            
            % For brick elements of size [dx, dy, dz] the element local mass matrix is
            dx = dimens(1);
            dy = dimens(2);
            dz = dimens(3);
            
            V = dx*dy*dz;
            L = [  1/27,  1/54,   1/108,  1/54,   1/54,   1/108,  1/216, 1/108;
                1/54,  1/27,   1/54,   1/108,  1/108,  1/54,   1/108, 1/216;
                1/108, 1/54,   1/27,   1/54,   1/216,  1/108,  1/54,  1/108;
                1/54,  1/108,  1/54,   1/27,   1/108,  1/216,  1/108, 1/54;
                1/54,  1/108,  1/216,  1/108,  1/27,   1/54,   1/108, 1/54;
                1/108, 1/54,   1/108,  1/216,  1/54,   1/27,   1/54,  1/108;
                1/216, 1/108,  1/54,   1/108,  1/108,  1/54,   1/27,  1/54;
                1/108, 1/216,  1/108,  1/54,   1/54,   1/108,  1/54,  1/27] * V;
            
        end
        
        
        
    end     % methods
end      % Identity3DQ1


