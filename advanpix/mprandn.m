function r = mprandn(varargin)
%MPRANDN Normally distributed pseudorandom multiple-precision numbers
%
%    All built-in special cases are supported:
%
%       r = mprandn(n)
%       r = mprandn(m,n)
%       r = mprandn([m,n])
%       r = mprandn(m,n,p,...)
%       r = mprandn([m,n,p,...])
%       r = mprandn
%       r = mprandn(size(A))
%       r = mprandn(..., 'double') % last argument is ignored
%       r = mprandn(..., 'single') % last argument is ignored       
%    
%    Precision is controlled by mp.Digits
%
%    See also mp.Digits, mprand, RANDN 

%    Copyright (c) 2006 - 2022 Advanpix LLC.

    r = mpimpl(2006,varargin{:});
end
