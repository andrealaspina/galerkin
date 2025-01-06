classdef (InferiorClasses = {?mp}) modifiedGramSchmidt
    %ModifiedGramSchmidt Performs basis expansion using MGS
    %orthogonalization algorithm.
    %
    %    Modified Gram-Schmidt is computationally intensive algorithm requiring 
    %    element-wise access to matrix elements, etc. Weak vectorization makes it 
    %    very slow when implemented in MATLAB language, especially for custom 
    %    class types (e.g. multiprecision floating-point numbers).
    %
    %    The purpose of the class is move all the heavy computations into C++ 
    %    and provide fast MGS, heavily optimized for multiprecision case.
    %    
    %    The implementation is tuned for usage in Krylov-like iterative algorithms 
    %    allowing adding vectors one by one into basis. 
    %
    %    EXAMPLE. Krylov basis expansion routine might look like (see mpeigs.m):
    %
    %      function [Q, H] = expandKrylov(A, Q, H, sk, ek)
    %      % Expands Krylov subspace.
    %      %
    %      %  The function contruct the sk+1, sk+2, ..., ek_th column of Q.
    %      %  A * Q(:, 1:ek) = Q(:, 1:ek+1) * H
    %      % Parameters:
    %      %   sk       start index
    %      %   ek       end index
    %      % Returns:
    %      %   Q        the expanded orthornormal matrix with dimension [n x ek+1]
    %      %   H        dimension [ek+1 x ek], the upper [ek x ek] block is Hessenberg
    %      
    %        mgs = modifiedGramSchmidt(Q,H,sk,ek); % initialization
    %         
    %        for j = sk:ek
    %          v = A(mgs.lastVector());  % apply linear operator to last vector - Q(:,j)
    %          mgs.addVector(v);         % orthogonolize new vector and add it to the basis.
    %        end
    %        
    %        [Q,H] = mgs.factors();      % extract factors
    %
    %        delete(mgs);                % free resources to avoid memory leaks
    %      end    
    
    %   Part of Multiprecision Computing Toolbox for MATLAB.
    %   Copyright (c) 2006 - 2022 Advanpix LLC.  
    
    properties (SetAccess = public)
        id
    end % properties
    
    methods
        
        function this = modifiedGramSchmidt(Q,H,start,stop)
        %MODIFIEDGRAMSCHMIDT Initializes resources needed to perform
        %Modified Gram-Schmidt orthogonalization. 
        %
        %    Matrices Q and H should be pre-allocated. Number of rows in Q
        %    defines vector length.
        %
        %    start  - index of the last valid orthogonal vector stored in Q.
        %    stop+1 - maximum index to expand basis to. 
        %    Vectors start+1:stop+1 will be added to the basis in Q, total
        %    number of vectors to be added is stop-start+1.
        %
        %    Q and H will have following sizes on exit: Q(n, stop+1) and H(stop+1,stop)
        %
        %    EXAMPLE (see mpeigs.m for complete example):
        %
        %    Q = zeros(n, m+1, class_t); 
        %    H = zeros(m+1, m, class_t);
        %    Q(:,1) = v0 / norm(v0);
        %
        %    [Q, H] = expandKrylov(A, Q, H, 1, b); % build initial basis
        %
        %
        %    See also ADDVECTOR, LASTVECTOR, FACTORS.
 
           if nargin < 4, error('MCT:modifiedGramSchmidt:constructor','Wrong number of input arguments'); end
           this.id = mpimpl(20000,Q,H,start,stop);
        end
        
        function v = lastVector(this)
        %LASTVECTOR Returns last vector in the basis Q(:,end)
        %
        %    See also ADDVECTOR, FACTORS.
        
           v = mpimpl(20001,this.id); 
        end

        function r = addVector(this,v)
        %ADDVECTOR Adds candidate vector 'v' to basis Q using modified
        %Gram-Schmidt orthogonalization.
        %
        %    See also LASTVECTOR, FACTORS.
        
           r = mpimpl(20002,this.id,v); 
        end
        
        function [Q,H] = factors(this)
        %FACTORS Returns computed factors Q and H with added vectors 
        %using ADDVECTOR.
        %
        %    See also ADDVECTOR, FACTORS.
        
           [Q,H] = mpimpl(20003,this.id); 
        end
        
        %% Destructor
        function delete(this)
        %DELETE Frees allocated resources
        %
        %    'this' is an object of type modifiedGramSchmidt
        %
        %    Do not call this function twice on the same object
        
           if nargin < 1,error('MCT:modifiedGramSchmidt:delete', 'Wrong number of input arguments'); end
           mpimpl(20004, this.id);
        end
        
    end
    
end

