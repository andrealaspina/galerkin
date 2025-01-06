function varargout = nextbelow(varargin)
%NEXTBELOW Returns the next representable floating-point value towards -Infinity
%
%    Examples:
%
%    >> fprintf('%.16e\n',nextabove(1))
%    1.0000000000000002e+00
%
%    >> fprintf('%.16e\n',nextbelow(1))
%    9.9999999999999989e-01
%
%    >> fprintf('%.16e\n',nextabove(-1))
%    -9.9999999999999989e-01
%
%    >> fprintf('%.16e\n',nextbelow(-1))
%    -1.0000000000000002e+00
%
%    See also NEXTABOVE.

%    Copyright 2008-2021 Advanpix LLC.

        [varargout{1:nargout}] = mpimpl(608,varargin{:});            
end
