function h = createSSAtoms(Ns,opt)
%% Create stable spline atoms.
% Generate sin(k*pi*alpha^t) atoms, where k is odd and randomly selected.
% nr: number of atoms

t = (1:Ns)'-1;
if opt.Type=="random"
   Alpha = 11;
   while Alpha>=1
      Alpha = rand(size(opt.Alpha));
   end
else
   Alpha = opt.Alpha;
end
at = (Alpha(1).^t)/2;

nr = opt.NumAtoms;
if opt.Randomize
   k = randperm(nr*10,nr);
else
   k = 1:1:nr;
end
k = 2*k-1; % odd entries only
h = sin(pi*k(ones(Ns,1),:).*at);
if ~isscalar(opt.Alpha)
   for ct = 2:numel(Alpha)
      at = (Alpha(ct).^t)/2;
      hct = sin(pi*k(ones(Ns,1),:).*at);
      h = [h,hct];
   end
end
if opt.IncludeConstant
   h = [h, ones(Ns,1)/Ns*2];
end
% delays
hh = [diag(ones(1,10),-1);zeros(Ns-11,11)];
hh = hh(:,1:end-1)./(1:10);
h = [h,hh];
%h = [h, -h];
%aa = anorm(h);
%h = h./aa';

