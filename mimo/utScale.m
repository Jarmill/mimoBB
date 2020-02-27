function  [A,Left,S] = utScale(A,Left,rc)
% rc: 'row'  or 'col'

N = [A, Left];

if strcmp(rc,'row')
   n = zeros(size(N,1),1);
   for i = 1:size(N,1)
      n(i) = norm(N(i,:));
   end
   S = diag(1./n);
   N2 = S*N;
else 
   n = zeros(size(N,2),1);
   for i = 1:size(N,2)
      n(i) = norm(N(:,i));
   end
   S = diag(1./n);
   N2 = N*S;
end
A = N2(:,1:end-1); Left = N2(:,end);
end