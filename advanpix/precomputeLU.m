classdef (InferiorClasses = {?mp}) precomputeLU
    %precomputeLU Computes and stores LU for further usage.
    %
    %    Allows LU decomposition re-use with minimum overhead.
    %    LU factors are not generated nor used explicitly - which saves
    %    memory and time.
    %
    %    Crucial for large scale sparse matrices and related iterative algorithms.
    %
    %    Usage example (pseudo-code):
    %
    %    Build matrix A
    %    A = [...]
    %    
    %    Compute factorization
    %    F = precomputeLU(A);
    %     
    %    Use the factorization of F several times
    %    for i=1:many
    %      [...]
    %      x = F\b;
    %      [...]
    %    end
    %
    %    delete(F);
    %   
    %   Suggested by Philippe Marti from University of Colorado Boulder.
    
    %   Part of Multiprecision Computing Toolbox for MATLAB.
    %   Copyright (c) 2006 - 2022 Advanpix LLC.  
    
    properties (SetAccess = public)
        id
    end
    
    methods
        
        %% Constructor
        function this = precomputeLU(A)
        %precomputeLU Compute and store LU factorization of A.
        %
        %    A must be of multiprecision type. 
        %
        %    Example:
        %
        %    F = precomputeLU(A);
 
           if nargin < 1, error('MCT:precomputeLU: Wrong number of input arguments'); end
           
           if ~ismp(A), this.id = mpimpl(11000, mp(A));
           else         this.id = mpimpl(11000,     A); 
           end;
        end
        
        %% Operator overloads
        function r = mldivide(F,B)
        %MLDIVIDE Solve systems of linear equations Ax = B for x using
        %precomputed LU decomposition of A.
        %
        %    F is precomupted LU decomposition of A.
        %
        %    See also MLDIVIDE.
           r = mpimpl(11001,F.id,B); 
        end
        
        %% Destructor
        function delete(F)
        %DELETE Frees resources allocated by precomputed LU decomposition
        %
        %    F is precomupted LU decomposition
        %
        %    Do not call this function twice on the same object
        
           if nargin < 1,error('MCT:precomputeLU:delete: Wrong number of input arguments'); end
        
           mpimpl(11002, F.id);
        end
        
    end
    
end

