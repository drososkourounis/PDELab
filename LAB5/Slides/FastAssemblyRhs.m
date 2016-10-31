  % this is just a vector form of the source function
  f = makeSource(mesh.Points);
  I = mesh.Elements'; 
  F = f(I);
  b = Me*F;

  I = I(:);
  J = ones(mesh.N_v*mesh.N_e, 1);
  
  B = sparse(I, J, b(:), mesh.N, 1);
  
  b = full(B);

end

 
