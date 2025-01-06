function r = mprand(varargin)
%MPRAND Uniformly distributed multiple-precision random numbers.
%
%    All built-in special cases are supported:
%
%       r = mprand(n)
%       r = mprand(m,n)
%       r = mprand([m,n])
%       r = mprand(m,n,p,...)
%       r = mprand([m,n,p,...])
%       r = mprand
%       r = mprand(size(A))
%       r = mprand(..., 'double') % last argument is ignored
%       r = mprand(..., 'single') % last argument is ignored       
%    
%    Based on Mersenne Twister pseudo-random number generator algorithm.
% 
%    Precision is controlled by mp.Digits
%
%    See also mp.Digits, mprandn, RAND 

%    Copyright (c) 2006 - 2022 Advanpix LLC.

    r = mpimpl(2003,varargin{:});
end
