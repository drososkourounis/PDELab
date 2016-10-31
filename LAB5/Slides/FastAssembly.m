function [A] = assembleFast(mesh)
  ele_nodes= mesh.N_e*mesh.N_v^2;
  K = zeros(1, ele_nodes); iK = zeros(1, ele_nodes); 
  jK = zeros(1, ele_nodes); i = 1:mesh.N_v;
  % ig = [1..nv, 1..nv, ... 1..nv] nv times
  ig       = repmat(i,1,mesh.N_v);
  % jg = [1 1 ... 1, 2 2 ... 2, ... , nv nv ... nv]
  jg       = repmat(1:mesh.N_v, mesh.N_v,1);
  jg       = jg(:)'; b = zeros(mesh.nodes,1);
  k        = 1:mesh.N_v^2;
  for e = 1:mesh.N_e
    Me    = makeMe(mesh.dx, mesh.dy); me = Me(:);
    Ke    = makeKe(mesh.dx, mesh.dy); ke = Ke(:);
    I     = mesh.Elements(e,:); iK(k) = I(ig); jK(k) = I(jg);
    x_e   = mesh.Points(I,:); f_e   = makeSource(x_e);
    b(I)  = b(I) + Me*f_e;
    K(k)  = ke + me; k   = k + mesh.N_v^2;
  end
  A = sparse(iK, jK, K);
end


