function r = sym2mp(s, precision)
%SYM2MP Converts sym object to multiprecision floating point object
%compatible with Multiprecision Computing Toolbox.
%
%    Parameter s should contain numeric value or mathematical expression (without symbolic variables).
% 
%    Precision is controlled by second argument ( = mp.Digits by default).
%    
%    Example:
%
%         B = vpa(rand(2,2));
%         A = sym2mp(B)
% 
%         A = 
% 
%                 0.42176128262627499143633258427144        0.79220732955955441845219411334256    
%                 0.91573552518906708996837551239878         0.9594924263929029972786111102323    
%
%
%    See also mp.Digits

%    Copyright (c) 2006 - 2022 Advanpix LLC. 

    if nargin < 2, precision = mp.Digits(); end;
    
    r = mp(zeros(size(s)), precision);
    for n=1:numel(s), r(n) = mp(char(s(n)), precision); end;
end
