function R = rosser(classname)
%ROSSER Classic symmetric eigenvalue test problem.
%   This matrix was a challenge for many matrix eigenvalue algorithms.
%   But LAPACK's DSYEV routine used in MATLAB has no trouble with it.
%   The matrix is 8-by-8 with integer elements.
%   It has:
%       * A double eigenvalue.
%       * Three nearly equal eigenvalues.
%       * Dominant eigenvalues of opposite sign.
%       * A zero eigenvalue.
%       * A small, nonzero eigenvalue.
%
%   ROSSER(CLASSNAME) produces a matrix of class CLASSNAME.
%   CLASSNAME must be either 'single' or 'double' (the default).

%   Copyright 1984-2005 The MathWorks, Inc.
%   Copyright 2008-2021 Advanpix LLC.

if nargin < 1, classname = 'double'; end
if mp.OverrideDoubleBasicArrays(), classname = 'mp'; end
    
R  = [ 611.  196. -192.  407.   -8.  -52.  -49.   29.
       196.  899.  113. -192.  -71.  -43.   -8.  -44.
      -192.  113.  899.  196.   61.   49.    8.   52.
       407. -192.  196.  611.    8.   44.   59.  -23.
        -8.  -71.   61.    8.  411. -599.  208.  208.
       -52.  -43.   49.   44. -599.  411.  208.  208.
       -49.   -8.    8.   59.  208.  208.   99. -911.
        29.  -44.   52.  -23.  208.  208. -911.   99.];
R = cast(R,classname);
