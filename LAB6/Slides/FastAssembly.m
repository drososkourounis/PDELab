function [A] = assembleFast(Elements, Points, dx, dy, nx, ny)
  Me       = makeLocalMass(dx, dy);
  Le       = makeLocalLaplacian(dx, dy);
  me       = Me(:);
  le       = Le(:);
  n        = (nx+1)*(ny+1); nel = nx*ny; nv = 4;
  vol_size = nel*nv^2;
  K        = zeros(1, vol_size); iK   = zeros(1, vol_size);
  i        = 1:nv;               j    = ones(1, nv);
  % ig = [1..nv, 1..nv, ... 1..nv] nv times
  ig       = repmat(i,1,nv);
  % jg = [1 1 ... 1, 2 2 ... 2, ... , nv nv ... nv]
  jg       = repmat(1:nv,nv,1);
  jg       = jg(:)';
  b        = zeros(n,1);
  k        = 1:nv^2;
  
  for e = 1:nel
      I     = Elements(e,:);
      x_e   = Points(I,:); f_e   = makeSource(x_e);
      b(I)  = b(I) + Me*f_e;
      iK(k) = I(ig); jK(k) = I(jg);
      Kg(k) = le + me; k   = k + nv^2;
  end
end


