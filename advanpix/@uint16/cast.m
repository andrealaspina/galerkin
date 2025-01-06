function r = cast(A,newclass,varargin) %#ok<INUSL>
    %CAST  Cast a variable to a different data type or class with support of multiprecision type.
    %   B = CAST(A,NEWCLASS) casts A to class NEWCLASS. A must be convertible to
    %   class NEWCLASS.
    %
    %   B = CAST(A,'like',Y) converts A to the same data type and sparsity as the 
    %   variable Y. If A and Y are both real, then B is also real. B is complex
    %   otherwise.
    %
    %   Proposed here: http://mct.userecho.com/topics/46-cast-method-may-also-be-useful-to-others/
    %
    %   Example:
    %      a = int8(5);
    %      b = cast(a,'mp');
    %
    %   See also CLASS.
    
    %   Copyright 2008-2021 Advanpix LLC.

    if nargin < 2, error('Not enough input arguments.'); end;
    if nargin > 3, error('Too many input arguments.'); end;

    sparsity = '';
    if strcmp(newclass,'like')
        if nargin < 3, error('Not enough input arguments.'); end;
        if nargin > 3, error('Too many input arguments.'); end;
        newclass = class(varargin{1}); 
        if issparse(varargin{1}), sparsity = 'sparse'; end;
    end

    % Return mp-objects for all floating-point types    
    if mp.OverrideDoubleBasicArrays() && (strcmpi(newclass,'double') || strcmpi(newclass,'single'))
       newclass = 'mp';
    end
    
    % Silent fall-back to double() is questionable, and can be used 
    % only in certain conditions (e.g. for conversion to
    % non-trivial types - 'ss', 'tf', etc. in older versions of MATLAB) 
    % 
    % Of course, it is better to implement explicit converters for
    % such types.
    %
    % We will see error message for now, uncomment try-catch only if you
    % know what you are doing.
    %try 
        r = eval([newclass,'(',sparsity,'(A)',')']);
    %catch, r = eval([newclass,'(',sparsity,'(double(A))',')']); end
end
