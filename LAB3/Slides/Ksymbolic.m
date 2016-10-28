function I = Laplace2DSymbolic()
    xi  = sym('xi', 'real'); eta = sym('eta', 'real');
    dx  = sym('dx', 'real'); dy  = sym('dy', 'real');
    u1 = sym('u1', 'real');  u2 = sym('u2', 'real');
    u3 = sym('u3', 'real');  u4 = sym('u4', 'real');
     
    c=[-1 -1; 1 -1; 1  1; -1  1];
    J(1,1) = 0.5*dx; J(2,2) = 0.5*dy;
    
    for i=1:4
        N(i) = 1/4*( 1+c(i,1)*xi )*( 1+c(i,2)*eta );
    end
    
    Nx  = diff(N, 'xi'); Ny  = diff(N, 'eta');
    dN = [Nx; Ny];
    U   = (u1*N(1) + u2*N(2) + u3*N(3) + u4*N(4))^2;
    Jd  = inv(J) * dN;
    F   = U*det(J)*Jd'*Jd;
    I   = int(int(F, 'xi', -1, 1), 'eta', -1, 1);
    I   = simplify(dx*dy*I);
end
        

