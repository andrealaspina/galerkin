function r = realmax( classname )
%REALMIN Smallest positive normalized floating point number.
% 
%   REALMAX(classname) returns the largest finite floating point number
%   corresponding to numeric type CLASSNAME.
%
%   Suports multiprecision numbers (CLASSNAME = 'mp') otherwise equivalent
%   to built-in REALMAX.
%
%    See also REALMIN, EPS

    if strcmpi('mp',classname)
        r = mpimpl(2002);          
    else
        r = builtin('realmax',classname);
    end
end
