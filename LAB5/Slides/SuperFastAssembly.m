function [As,b] = assembleFast(mesh)
  Me    = makeMe(mesh);
  Ke    = makeKe(mesh);
  me    = Me(:); 
  ke    = Ke(:);
  M     = repmat(me, mesh.N_e, 1);
  K     = repmat(ke, mesh.N_e, 1);
  i=1:mesh.N_v; j=ones(1, mesh.N_v);

  % ig = [1..N_v, 1..N_v, ... 1..N_v] N_v times
  ig=repmat(i, 1, mesh.N_v);

  % jg = [1 1 ... 1, 2 2 ... 2, ... , N_v N_v ... N_v]
  jg=repmat(1:mesh.N_v, mesh.N_v, 1); jg=jg(:)';

  iA    = mesh.Elements(:,ig)';
  jA    = mesh.Elements(:,jg)';
  A     = K + M;
  As    = sparse(iA(:), jA(:), A, mesh.N, mesh.N);

