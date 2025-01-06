function varargout = arrayCreationOverload(funcname,varargin)
     % Generic overload for array creation functions in @double and @single
     % directories.
     %
     % Called when all numeric parameters are of standard types.
     %
     % Notes:
     % 1. The function is called only if all numeric parameters are of
     %    standard types. Otherwise (if any of parameters is of 'mp' type) -
     %    overload mp.funcname is called instead.
     %
     % 2. Case of funcname(...,'like',p), where class(p) == 'mp' 
     %    is already handled by overload in mp class itself (see above).
     %
     % 3. There is exceptional case, when zeros cannot be overriden:
     %    X = funcname            % no parameters
     %    X = funcname('like',1)  % first parameter is string
     %    MATLAB just returns scalar double without calling 'funcname'.
     %    To overload the cases, we need to put 'funcname' on path, but then
     %    we will get warnings.
     %
     % Copyright (c) 2008-2021 by Advanpix LLC.

     % Always create multiprecision array when 'mp' is required explicitely:
     % funcname(...,'mp')
     if ischar(varargin{end}) && strcmpi(varargin{end},'mp')
         [varargout{1:nargout}] = mp(builtin(funcname,varargin{1:end-1}));
     else
        if mp.OverrideDoubleBasicArrays()
            if ischar(varargin{end}) && ~(strcmpi(varargin{end},'double') || strcmpi(varargin{end},'single')) || ...
                    (nargin > 2 && ischar(varargin{end-1}) && strcmpi(varargin{end-1},'like') ...
                    && ~(isa(varargin{end},'double') || isa(varargin{end},'single')))
                
                % Call built-in only when non-'double'/'single' output is required:
                % X = funcname(___,typename), where typename ~= 'double'/'single'    
                % X = funcname(___,'like',p), where class(p) ~= 'double'/'single'
                [varargout{1:nargout}] = builtin(funcname,varargin{:});
            else
                
                % Create multiprecision array in all other cases:
                % X = funcname(n)
                % X = funcname(sz1,...,szN)
                % X = funcname(sz)
                % X = funcname(___,typename), where typename == 'double'/'single'    
                % X = funcname(___,'like',p), where class(p) == 'double'/'single'
                [varargout{1:nargout}] = mp(builtin(funcname,varargin{:}));
            end
        else
             % Call default function if no override is allowed
             [varargout{1:nargout}] = builtin(funcname,varargin{:}); 
        end
     end
end
