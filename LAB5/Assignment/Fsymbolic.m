function Fe = Fsymbolic()
  syms x y real;
  syms a1 a2 a3 b1 b2 b3 c1 c2 c3 real;

  T = [ a1 b1 c1; a2 b2 c2; a3 b3 c3];
  N =  T * [x; y; 1];
  f = x^2 * y^2 * (2*x-3) * (2*y-3);

  Nf = N*f;

 
  Fe = int(int(Nf, 'xi', 0, 1-eta), 'eta',0, 1)
