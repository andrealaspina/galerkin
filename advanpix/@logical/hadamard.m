function r = hadamard(n, classname)
    if nargin < 2, classname = 'double'; end    
    if mp.OverrideDoubleBasicArrays(), classname = 'mp'; end
    
    r = hadamardImpl(n, classname);
end

function H = hadamardImpl(n,classname)
%HADAMARD Hadamard matrix.
%   HADAMARD(N) is a Hadamard matrix of order N, that is,
%   a matrix H with elements 1 or -1 such that H'*H = N*EYE(N).
%   An N-by-N Hadamard matrix with N > 2 exists only if REM(N,4)=0.
%   This function handles only the cases where N, N/12 or N/20
%   is a power of 2.
%
%   HADAMARD(N,CLASSNAME) produces a matrix of class CLASSNAME.
%   CLASSNAME must be either 'single' or 'double' (the default).

%   Nicholas J. Higham
%   Copyright 1984-2005 The MathWorks, Inc.

%   Reference:
%   S. W. Golomb and L. D. Baumert, The search for Hadamard matrices,
%   Amer. Math. Monthly, 70 (1963) pp. 12-17.

if nargin < 2, classname = 'double'; end

[f,e] = log2([n n/12 n/20]);
k = find(f==1/2 & e>0);
if min(size(n)) > 1 || isempty(k)
   error(message('MATLAB:hadamard:InvalidInput'));
end
e = e(k)-1;

if k == 1        % N = 1 * 2^e;
   H = eval([classname,'(1)']);

elseif k == 2    % N = 12 * 2^e;
   H = [ones(1,12,classname); ones(11,1,classname) ...
        toeplitz([-1 -1 1 -1 -1 -1 1 1 1 -1 1],[-1 1 -1 1 1 1 -1 -1 -1 1 -1])];

elseif k == 3    % N = 20 * 2^e;
   H = [ones(1,20,classname); ones(19,1,classname)   ...
        hankel([-1 -1 1 1 -1 -1 -1 -1 1 -1 1 -1 1 1 1 1 -1 -1 1], ...
               [1 -1 -1 1 1 -1 -1 -1 -1 1 -1 1 -1 1 1 1 1 -1 -1])];
end

%  Kronecker product construction.
for i = 1:e
    H = [H  H
         H -H]; %#ok
end
end