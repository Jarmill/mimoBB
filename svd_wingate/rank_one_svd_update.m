function [Up,Sp,Vp] = rank_one_svd_update( U, S, V, a, b, force_orth )
% function [Up,Sp,Vp] = rank_one_svd_update( U, S, V, a, b, force_orth )
%
% Given the SVD of
%
%   X = U*S*V'
%
% update it to be the SVD of
%N
%   X + ab' = Up*Sp*Vp'
%
% that is, implement a rank-one update to the SVD of X.
%
% Depending on a,b there may be considerable structure that could
% be exploited, but which this code does not.
%
% The subspace rotations involved may not preserve orthogonality due
% to numerical round-off errors.  To compensate, you can set the
% "force_orth" flag, which will force orthogonality via a QR plus
% another SVD.  In a long loop, you may want to force orthogonality
% every so often.
%
% See Matthew Brand, "Fast low-rank modifications of the thin
% singular value decomposition".
%
% D. Wingate 8/17/2007
%

  %current_rank = size( U, 2 );
  current_rank = size( V, 1); %adding columns, not rows

  % P is an orthogonal basis of the column-space
  % of (I-UU')a, which is the component of "a" that is
  % orthogonal to U.
  m = U' * a;
  p = a - U*m;
  Ra = sqrt(p'*p);
  P = (1/Ra)*p;

  % XXX this has problems if a is already in the column space of U!
  % I don't know what to do in that case.
%   if ( Ra < 1e-13 )
%     fprintf('------> Whoa! No orthogonal component of m!\n');
%   end;
%   
  % Q is an orthogonal basis of the column-space
  % of (I-VV')b.
  n = V' * b;
  q = b - V*n;
  Rb = sqrt(q'*q);
  Q = (1/Rb)*q;

  %This will fire incredibly often. Get used to it.
%   if ( Rb < 1e-13 )
%     fprintf('------> Whoa! No orthogonal component of n!\n');
%   end;
  
  %
  % Diagonalize K, maintaining rank
  %

  % XXX note that this diagonal-plus-rank-one, so we should be able
  % to take advantage of the structure!
  z = zeros( size(m) );

  %adding column (special case). can be 
  if all(n == 0)    
   %Matrix is [S m; 0 p], of a broken arrowhead form
   rho = Ra*Rb;
   K = broken_arrowhead(diag(S), m, rho);
   [tUp,tSp,tVp] = svds( K, current_rank );
  else
    K = [ S z ; z' 0 ] + [ m; Ra ]*[ n; Rb ]';
    [tUp,tSp,tVp] = svds( K, current_rank );
  end
  
 
  % Now update our matrices!
  %
  % these multiplications can get expensive
  
  Sp = tSp;

  %this is a hack. Find out why
  %happens when only one nonzero element per column of matrix.
  if Ra < 1e-9

      if issparse(U)
          [i, j, v] = find(U);
          Uaug = sparse(i, j, v, size(U, 1), size(U, 2) + 1);
          Up = Uaug*sparse(tUp);
      else
          zP = zeros(size(P));
          Uaug = [U zP];
          Up = Uaug * tUp;
      end
      %Up = Uaug * tUp;
      %Up = tUp(1:end-1, :);
  else
      Uaug = [U P];
      if issparse(U)
          %Up = Uaug * sparse(tUp);
          %bottleneck is over here. Figure it out, and/or get better code.
          Up = sparse(Uaug * tUp);
      else
          Up = Uaug * tUp;
      end
      %Up = [ U P ] * tUp;
  end
  
  if Rb < 1e-9
      Vp = [V z] * tVp;
      %Vp = tVp(:, 1:end-1);
  else
      Vp = [ V Q ] * tVp;
  end

  % The above rotations may not preserve orthogonality, so we explicitly
  % deal with that via a QR plus another SVD.  In a long loop, you may
  % want to force orthogonality every so often.

  if ( force_orth )
    [UQ,UR] = qr( Up, 0 );
    [VQ,VR] = qr( Vp, 0 );
    [tUp,tSp,tVp] = svds( UR * Sp * VR', current_rank );
    if issparse(Up)
        Up = UQ * sparse(tUp);
    else
        Up = UQ * tUp;
    end
    Vp = VQ * tVp;
    Sp = tSp;
  end;
  
return;