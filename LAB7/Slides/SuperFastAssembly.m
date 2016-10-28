function [As,b] = assembleFast(e, mesh)
  Me    = makeMe(e, mesh);
  Ke    = makeKe(e, mesh);
  me    = Me(:); 
  ke    = Ke(:);
  n     = mesh.N; nel = mesh.Ne; nv = mesh.Nv;
  M     = repmat(me, nel, 1);
  K     = repmat(ke, nel, 1);
  i=1:nv; j=ones(1,nv);

  % ig = [1..nv, 1..nv, ... 1..nv] nv times
  ig=repmat(i,1,nv);

  % jg = [1 1 ... 1, 2 2 ... 2, ... , nv nv ... nv]
  jg=repmat(1:nv,nv,1); jg=jg(:)';

  iA    = mesh.Elements(:,ig)';
  jA    = mesh.Elements(:,jg)';
  A     = K + M;
  As    = sparse(iA(:), jA(:), A, n, n);

