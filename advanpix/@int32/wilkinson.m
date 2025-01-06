function W = wilkinson(n, classname)
%WILKINSON Wilkinson's eigenvalue test matrix.
%   WILKINSON(n) is J. H. Wilkinson's eigenvalue test matrix, Wn+.
%   It is a symmetric, tridiagonal matrix with pairs of nearly,
%   but not exactly, equal eigenvalues.
%   The most frequently used case is WILKINSON(21).
%   For example, WILKINSON(7) is
%
%          3  1  0  0  0  0  0
%          1  2  1  0  0  0  0
%          0  1  1  1  0  0  0
%          0  0  1  0  1  0  0
%          0  0  0  1  1  1  0
%          0  0  0  0  1  2  1
%          0  0  0  0  0  1  3
%
%   WILKINSON(CLASSNAME) produces a matrix of class CLASSNAME.
%   CLASSNAME must be either 'single' or 'double' (the default).

%   Copyright 1984-2005 The MathWorks, Inc.
%   Copyright 2008-2021 Advanpix LLC.

    if nargin < 2, classname = 'double'; end
    if mp.OverrideDoubleBasicArrays(), classname = 'mp'; end
    
    m = cast((n-1)/2,classname);
    e = ones(n-1,1,classname);
    W = diag(abs(-m:m)) + diag(e,1) + diag(e,-1);
end