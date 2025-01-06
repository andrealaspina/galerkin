function P = pascal(n,k,classname)
%PASCAL Pascal matrix.
%   PASCAL(N) is the Pascal matrix of order N: a symmetric positive
%   definite matrix with integer entries, made up from Pascal's
%   triangle.  Its inverse has integer entries.
%   PASCAL(N).^r is symmetric positive semidefinite for all
%   nonnegative r.
%
%   PASCAL(N,1) is the lower triangular Cholesky factor (up to signs
%   of columns) of the Pascal matrix.  It is involutory (is its own
%   inverse).
%
%   PASCAL(N,2) is a rotated version of PASCAL(N,1), with sign flipped
%   if N is even, which is a cube root of the identity.
%
%   PASCAL(N,CLASSNAME) or PASCAL(N,K,CLASSNAME) produces a matrix of
%   class CLASSNAME.  CLASSNAME must be either 'single' or 'double'
%   (the default).

%   Reference:
%   N. J. Higham, Accuracy and Stability of Numerical Algorithms,
%   Second edition, Society for Industrial and Applied Mathematics,
%   Philadelphia, 2002, Sec. 28.4.
%
%   Nicholas J. Higham

%   Copyright 1984-2005 The MathWorks, Inc.
%   Copyright 2008-2021 Advanpix LLC.

if nargin == 2 && ischar(k)
   classname = k; k = 0;
elseif nargin < 3
   classname = 'double';
end
if nargin < 2 || isempty(k), k = 0; end

if mp.OverrideDoubleBasicArrays(), classname = 'mp'; end

P = diag((-1).^(0:cast(n,classname)-1));
P(:,1) = ones(n,1,classname);

% Generate PASCAL(N,1): the Pascal Cholesky factor (up to signs).
for j=2:n-1
    for i=j+1:n
        P(i,j) = P(i-1,j) - P(i-1,j-1);
    end
end

if k == 0
    P = P*P';
elseif k == 2
   P = rot90(P,3);
   if n/2 == round(n/2), P = -P; end
elseif k ~= 1
    error(message('MATLAB:pascal:InvalidArg2'))
end
