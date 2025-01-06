function r = eps( varargin )
%EPS Spacing of floating point numbers.
% 
%    EPS(classname) returns the distance from 1.0 to the next larger 
%    floating-point number of CLASSNAME type.
%  
%    Suports multiprecision floating-point type (CLASSNAME = 'mp') otherwise equivalent
%    to built-in EPS.
%
%    See also REALMAX, REALMIN

    if strcmpi('mp',varargin{end}) || strcmpi('mp',class(varargin{end}))
        r = mpimpl(2000);
    else
        r = builtin('eps',varargin{:});
    end
end
