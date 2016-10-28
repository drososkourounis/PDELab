function [M,K,b] = assembleDiscreteOperators(mesh)
  N   = mesh.N;
  Ne  = mesh.N_e;
  M   = zeros(N,N); K = zeros(N,N);
  b   = zeros(N,1);
  for e=1:N_e
    Me = makeMe(e, mesh);
    Ke = makeKe(e, mesh);
    fp = makebe(e, mesh);
    I = mesh.Elements(e, :);
    M(I, I) = M(I, I) + Me;
    K(I, I) = K(I, I) + Ke;
    b(I)    = b(I) + Me*fp;
  end % e loop
end
