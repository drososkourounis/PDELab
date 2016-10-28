  % this is just a vector form of the source function
  f = makeSource(Points);
  I = Elements'; 
  F = f(I);
  b = Me*F;

  I = I(:);
  J = ones(nv*nel, 1);
  
  B = sparse(I, J, b(:), n, 1);
  
  [I, J, b] = find(B);

end

 
