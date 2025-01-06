function r = flintmax( classname )
%FLINTMAX Largest consecutive integer in floating point format.
%
%   FLINTMAX returns the largest consecutive integer in IEEE double
%   precision, which is 2^53. Above this value, double precision format
%   does not have integer precision, and not all integers can be represented 
%   exactly.
%
%   FLINTMAX('double') is the same as FLINTMAX.
%
%   FLINTMAX('single') returns the largest consecutive integer in
%   IEEE single precision, which is SINGLE(2^24).
%
%
%   Suports multiprecision numbers (CLASSNAME = 'mp') otherwise equivalent
%   to built-in FLINTMAX.
%
%   See also EPS, REALMAX, INTMAX.

    if strcmpi('mp',classname)
        r = 2^ceil(mp.Digits() * mp('log2(10)'));
    else
        r = builtin('flintmax',classname);
    end
end
