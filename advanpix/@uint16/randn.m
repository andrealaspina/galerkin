function varargout = randn(varargin)
    [varargout{1:nargout}] = arrayCreationOverload('randn',varargin{:});
end
