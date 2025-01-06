function r = realmin( classname )
%REALMIN Smallest positive normalized floating point number.
% 
%   REALMIN(classname) returns the smallest positive normalized floating point number
%   corresponding to numeric type CLASSNAME.
%
%   Suports multiprecision numbers (CLASSNAME = 'mp') otherwise equivalent
%   to built-in REALMIN.
%
%    See also REALMAX, EPS

    if strcmpi('mp',classname)
        r = mpimpl(2001);          
    else
        r = builtin('realmin',classname);
    end
end
