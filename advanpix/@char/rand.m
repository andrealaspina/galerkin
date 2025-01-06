function varargout = rand(varargin)
    [varargout{1:nargout}] = arrayCreationOverload('rand',varargin{:});
end
