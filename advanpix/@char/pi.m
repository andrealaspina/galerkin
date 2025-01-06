function r = pi(classname)
%PI   Mathematical constant, the ratio of a circle's circumference to its diameter, 
%    commonly approximated as 3.1415926535897932384626433832795....
%
%    PI(classname) returns Pi with the precision corresponding to 
%    requested floating-point CLASSNAME type.
%  
%    Suports multiprecision floating-point type (CLASSNAME = 'mp') otherwise equivalent
%    to built-in PI.
%
%    See also EPS, REALMAX, REALMIN

    if strcmpi('mp',classname)
        r = mp('pi');
    else
        r = builtin('pi');
    end
end
