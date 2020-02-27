function sys = utGenExampleSystem(rho,ny,nu,nx)
% Generate a MIMO LTI with 1>= radius >=rho
foo = true;
while foo
   sys = drss(nx,ny,nu);
   foo = any(abs(pole(sys))<rho);
end
sys.Ts = 1;
end