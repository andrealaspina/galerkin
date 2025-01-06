function H = hilb(n,classname)
%HILB   Hilbert matrix.
%   HILB(N) is the N by N matrix with elements 1/(i+j-1),
%   which is a famous example of a badly conditioned matrix.
%   See INVHILB for the exact inverse.
%
%   HILB(N,CLASSNAME) produces a matrix of class CLASSNAME.
%   CLASSNAME must be either 'single' or 'double' (the default).
%
%   This is also a good example of efficient MATLAB programming
%   style where conventional FOR or DO loops are replaced by
%   vectorized statements. 
%
%   See also INVHILB.
    
%   Copyright 1984-2013 The MathWorks, Inc.
%   Copyright 2008-2021 Advanpix LLC.

if nargin < 2, classname = 'double'; end    
if mp.OverrideDoubleBasicArrays(), classname = 'mp'; end

J = 1:cast(n,classname);
J = J(ones(n,1),:);
I = J';
E = ones(n,classname);
H = E./(I+J-1);

