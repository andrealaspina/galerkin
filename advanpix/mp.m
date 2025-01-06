classdef mp
%
%    Multiprecision floating-point numeric type. 
%    Allows computing with any required precision in MATLAB.
%    
%    Detailed User's Manual can be accessed through:
%    http://www.advanpix.com/documentation/users-manual/
%
%    Complete list of implemented functions:
%    http://www.advanpix.com/documentation/function-reference/
%
%
%    Below we present basic usage rules and examples. 
%
%    ***
%    I. Precision control is done by mp.Digits() function: 
%
%    mp.Digits(34);    % set default precision to 34 decimal digits (quadruple precision)
%    mp.Digits(100);   % set default precision to 100 decimal digits
%
%    ***
%    II. Multiprecision entities are created by using special constructor function:
%
%    mp()    creates multiprecision number of default precision (with zero value).
%    mp(x)   creates multiprecision entity from numeric array, matrix, expression or other mp-object.
% 
%    Examples:
%
%    mp.Digits(34);   % set default precision to 34 digits
% 
%    a = mp('pi')     % creates 34-digits accurate 'pi' by evaluating expression
%    a = 
%       3.141592653589793238462643383279503
%
%    A = magic(5,'mp'); % creates magic square using 34-digits precision.
%    
%    B = mp(randn([3,3,3,3])) % creates 4D array of normally distributed pseudo-random numbers.
%
%    ***
%    II. Once created multiprecision objects can be used in calculations exactly the same way 
%    as built-in 'double' or 'single' precision entities:
%
%    [U,S,V] = svd(A);
%    norm(A-U*S*V',1)  
%    ans = 
%          6.08593862426366529565689029856837e-32        
%
%    [V,D] = eig(A);
%    norm(A*V - V*D,1)               
%    ans = 
%          5.238529448733281520312260003831002e-32
%
%    Please check our User's Manual for more examples and details:
%    http://www.advanpix.com/documentation/users-manual/

%    Copyright (c) 2008 - 2022 Advanpix LLC.

    properties (SetAccess = public)
        id
    end % properties
    
    methods(Static)
        
        %% Global Toolbox Settings
        function varargout = Digits(varargin)
        %MP.DIGITS Controls default precision of computations.
        % 
        %   MP.DIGITS(digits) Setups precision to be of decimal 'digits'. 
        %                     The 'digits' should be an integer indicating number of decimal digits to
        %                     be used in all subsequent multiple-precision calculations.
        %
        %   MP.DIGITS()       Shows current precision (in decimal digits).
        %
        %   MP.DIGITS(X)      Shows precision of multiprecision variable X.         
        %
        %   Internally toolbox uses mp.Digits + mp.GuardDigits decimal digits in all computations.
        %
        %   By default, toolbox uses quadruple precision - mp.Digits(34). 
        %   This can be changed in MPSTARTUP routine which runs on every start of toolbox.
        %
        %   See also MP.GUARDDIGITS, MPSTARTUP

            [varargout{1:nargout}] = mpimpl(1,varargin{:});
        end

        function varargout = GuardDigits(varargin)
        %MP.GUARDDIGITS Controls number of guard digits used in computations.
        %
        %   MP.GUARDDIGITS()       Shows number of guard digits currently used.
        %
        %   MP.GUARDDIGITS(digits) Setups number of guard digits to use. 
        %                          The 'digits' should be an integer indicating number of additional decimal digits to
        %                          be used in all multiple-precision calculations. 
        %
        %   Total precision of computations are mp.Digits + mp.GuardDigits
        %
        %   By default, toolbox uses mp.GuardDigits(0) - no guard digits. 
        %   This can be changed in MPSTARTUP routine which runs on every start of toolbox.
        %
        %   See also MP.DIGITS, MPSTARTUP
            
             [varargout{1:nargout}] = mpimpl(3,varargin{:});
        end
        
        function varargout = ExtendConstAccuracy(varargin)
        %mp.ExtendConstAccuracy Controls accuracy auto-extension of MATLAB's built-in constants.
        %
        %   mp.ExtendConstAccuracy(true) enables accuracy auto-extension for 'double' precision 
        %   constants (pi, eps, etc.) to match precision selected in the toolbox.
        %
        %   mp.ExtendConstAccuracy(false) disables this feature. 
        % 
        %   Accepts 'true' or 'false' (default) as an argument. 
        %   Returns current settings if called without arguments.
        %
        %      >> mp.Digits(50);
        % 
        %      >> mp.ExtendConstAccuracy(false);
        %      >> mp(1/3) % 1/3 is computed in double precision, stored with 50 digits of precision.
        %      ans = 
        %          0.33333333333333331482961625624739099293947219848633
        % 
        %      >> mp.ExtendConstAccuracy(true);
        %      >> mp(1/3) % 1/3 constant is detected and re-computed with 50 digits of precision
        %      ans = 
        %          0.3333333333333333333333333333333333333333333333333324        
        %
        %    Please note, recommended way of using constants is:
        %
        %      >> mp('1/3')
        %      ans = 
        %          0.3333333333333333333333333333333333333333333333333324        
        %
        %   This feature is disabled by default, and it should be used with
        %   caution, especially with small/large numbers where error
        %   of mis-detection is high.
        %   
        %   See also MP.DIGITS, MPSTARTUP
            
             [varargout{1:nargout}] = mpimpl(5,varargin{:});
        end

        function varargout = FollowMatlabNumericFormat(varargin)
        %mp.FollowMatlabNumericFormat Controls if toolbox follows numeric formatting
        %preferences in MATLAB
        %
        %   MP.FOLLOWMATLABNUMERICFORMAT(true) makes toolbox follow numeric
        %   formatting preferences in MATLAB. 
        % 
        %   MP.FOLLOWMATLABNUMERICFORMAT(false) disables this feature.
        %   Numbers are displayed with all digits of precision.
        % 
        %   Accepts 'true' or 'false' (default) as an argument. 
        %   Returns current settings if called without arguments. 
        %
        %      >> format short
        %      
        %      >> mp.FollowMatlabNumericFormat(false); % show all digits
        %      >> mp('pi')
        %      ans = 
        %          3.1415926535897932384626433832795028
        %      
        %      >> mp.FollowMatlabNumericFormat(true);  % use MATLAB formatting
        %      >> mp('pi')
        %      ans = 
        %          3.1416
        % 
        %   See also MP.DIGITS, FORMAT, MPSTARTUP        
            
             [varargout{1:nargout}] = mpimpl(6,varargin{:});
        end
        
        function varargout = OverrideDoubleBasicArrays(varargin)
        %mp.OverrideDoubleBasicArrays Controls if toolbox replaces all 
        % built-in array-creation routines so that they produce multiprecision
        % output by default. 
        %
        %   The basic array-creation routines are:
        %   ones, zeros, eye, NaN, Inf, true, false, cast, rand, randn
        %
        %   mp.OverrideDoubleBasicArrays(true) allows toolbox to replace
        %   all such routines with multiprecision versions:
        %
        %     >> mp.OverrideDoubleBasicArrays(true);
        %     >> A = zeros(3,3);
        %     >> whos A
        %     Name      Size            Bytes  Class    Attributes
        %     A         3x3               376  mp                 
        %
        %   This mode is very helpful for conversion of large projects to
        %   multiprecision with minimal changes and for doing ALL computations 
        %   using the extended precision by default. 
        %
        %   However do not use the mode if you want to run only parts of
        %   computations in multiprecision.  
        % 
        %   mp.OverrideDoubleBasicArrays(false) disables this feature.
        %   All the array-creation functions stay intact and produce
        %   'double', 'single' or else output according to intial design:
        %
        %     >> mp.OverrideDoubleBasicArrays(false);
        %     >> A = zeros(3,3);
        %     >> whos A
        %     Name      Size            Bytes  Class    Attributes
        %     A         3x3                72  double                 
        %     
        %     >> A = zeros(3,3,'mp')  % explicit conversion is required
        %                             % to create mp-output
        %
        %   Accepts 'false' (default) or 'true' as an argument. 
        %   Returns current settings if called without arguments 
        % 
        %   See also MP.DIGITS, MPSTARTUP 
            
             [varargout{1:nargout}] = mpimpl(8,varargin{:});
        end

        function varargout = NumberOfThreads(varargin)
        %mp.NumberOfThreads Sets maximum number of threads to use in computations enabled with multi-core parallelism.
        %
        %   N = 0 (default):
        %   Sets number of threads = number of real hardware cores in the system.
        %   Each thread is pinned to execute on particular hardware core for best
        %   performance. This is optimal strategy for most of the users, who runs
        %   one instance of toolbox at a time.
        %   
        %   N ~= 0:
        %   Pushes toolbox to use exactly N cores, taking into account
        %   hyper-threaded cores as well. No thread affinity is applied.
        %   
        %   This is useful if you run several toolbox instances (e.g. with parfor). 
        %   In this case compute number of threads as:
        %   
        %      N = total_number_of_cores / number_of_matlab_workers
        %   
        %   Returns current setting if called without arguments.
        % 
        %   See also MP.DIGITS, MPSTARTUP 
            
             [varargout{1:nargout}] = mpimpl(9,varargin{:});
        end
        
        function varargout = EnableImplicitExpansion(varargin)
        %mp.EnableImplicitExpansion Controls if toolbox implicitly expands dimensions of the result in arithmetic operations.
        %
        %   false (default):
        %      Toolbox strictly follows linear algebra rules and implicit expansions
        %      are disabled (standard pre-R2016b behavior). Error message is shown if operands have
        %      non-matching shapes, e.g.:
        %      
        %        >> mp.EnableImplicitExpansion(false);
        %        >> mp('[1 2]') .* mp('[1; 2]')
        %           Error using  .*  (line 1665)
        %           Matrix dimensions must agree
        %      
        %   
        %   true:
        %      Toolbox silently expands dimensions of the result if
        %      operands have different (but compatible) shapes, e.g.:
        %
        %        >> mp.EnableImplicitExpansion(true);
        %        >> mp('[1 2]') .* mp('[1; 2]')
        % 
        %          ans = 
        % 
        %             1        2    
        %             2        4    
        %
        %   Returns current setting if called without arguments.
        % 
        %   See also MP.DIGITS, MPSTARTUP 
            
             [varargout{1:nargout}] = mpimpl(20100,varargin{:});
        end
        
        %% Toolbox Information & Unit tests
        function Info()
        %MP.INFO  Dispaly information about Multiprecision Computing Toolbox.  
        %        
        %    mp.Info shows version, build date, credits, license status and other 
        %    information about Multiprecision Computing Toolbox and its components.
        %
        %    See also MP.TEST, MP.DIGITS, MP.GUARDDIGITS

            mpimpl(4);          
        end

        function Test()
        %MP.TEST Run basic tests to check consistency of Multiprecision Computing Toolbox.  
        %
        %    mp.Test() runs simple tests to check functionality and
        %    correctness of the basic methods in toolbox.
        %
        %    Before every new release of toolbox we run comprehensive set
        %    of tests in-house to check (almost) every aspect of toolbox functionality. 
	    %    Running all tests takes ~3 hours for each platform.
        %
        %    See also MP.DIGITS, MP.GUARDDIGITS
            
            mptest();  % Run tests         
        end
        
        %% Gaussian Quadrature
        function [x,w] = GaussLegendre(n, a, b)
        %MP.GAUSSLEGENDRE Compute coefficients of Gauss-Legendre quadrature with arbitrary precision. 
        %
        %    [x,w] = mp.GaussLegendre(order, a, b) computes abscissae
        %    'x' and weights 'w' of Gauss-Legendre quadrature.
        %
        %    The Gauss-Legendre quadrature rule is used as follows:
        %
        %        Integral ( a <= x <= b ) f(x) dx
        %
        %    is to be approximated by
        %
        %        sum ( 1 <= i <= n ) w(i) * f(x(i))        
        %
        %    Parameters:
        %    'n' is the number of points in the quadrature rule.
        %    'a' is the left endpoint (default -1);
        %    'b' is the right endpoint(default +1);        
        %
        %    Precision is controlled by mp.Digits and by precision of the arguments
        %
        %    See also MP.GAUSSJACOBI, MP.GAUSSLAGUERRE, MP.GAUSSHERMITE, MP.GAUSSCHEBYSHEVTYPE1, MP.GAUSSCHEBYSHEVTYPE2, MP.GAUSSGEGENBAUER
        
            if nargin < 3,          b     = mp('+1');
                if nargin < 2,      a     = mp('-1');
                    if nargin < 1
                    error(  'MCT:GaussLegendre:NotEnoughInputs',...
                            'Not enough input arguments.' );
                    end  
                end
            end
            
            if ~ismp(n), n = mp(n); end; %#ok<*NOSEL>
            if ~ismp(a), a = mp(a); end;
            if ~ismp(b), b = mp(b); end; 
        
            % Find maximum precision among all (manually supplied) arguments
            if     nargin == 1, precision = mp.Digits(n);
            elseif nargin == 2, precision = mp.Digits([n a]);
            elseif nargin == 3, precision = mp.Digits([n a b]);
            end;
            
            old_digits = mp.Digits();
            mp.Digits(precision);

            % Compute nodes & weights for default interval: [-1,1]
            %
            % Toolbox approximates roots of Legendre polynomial by Tricomi
            % formula (3.4) + accuracy refinement by few Newton steps.
            %
            % N.Hale, A. Townsend. "Fast and Accurate Computation of Gauss-Legendre and Gauss-Jacobi Quadrature Nodes and Weights"
            % http://eprints.maths.ox.ac.uk/1629/1/finalOR79.pdf

            [x, w] = mpimpl(5000, n);

            % Scale coefficients for arbitrary interval [a,b].
            if (b ~= 1 || a ~= -1)
               scale = (b-a)/2;
               shift = (a+b)/2;
               x = shift + scale.*x;
               w = w.*scale;
            end

            mp.Digits(old_digits);            
        end

        function [x,w] = GaussHermite(n, alpha, a, b)
        %MP.GAUSSHERMITE Compute coefficients of Generalized Gauss-Hermite quadrature with arbitrary precision
        %
        %    [x,w] = mp.GaussHermite(order, alpha, a, b) computes abscissas
        %    'x' and weights 'w' of Generalized Gauss-Hermite quadrature.
        %
        %    The generalized Gauss-Hermite quadrature rule is used as follows:
        %
        %        Integral ( -inf < x < +inf ) |x-a|^alpha * exp( - b * ( x - a)^2 ) f(x) dx
        %      
        %    is to be approximated using 'x' and 'w' by
        %
        %         sum ( 1 <= i <= n ) w(i) * f(x(i))      
        %
        %    Parameters:
        %    'n' is the number of points in the quadrature rule.
        %    'alpha' is the parameter for the generalized Gauss-Hermite quadrature rule. 
        %           The value of alpha may be any real value greater than -1.0. 
        %           Specifying alpha = 0.0 (default) results in the basic (non-generalized) rule.
        %    'a' is the center point (default 0);
        %    'b' is the scale factor (default 1);        
        %
        %    Precision is controlled by mp.Digits and by precision of the arguments
        %
        %    See also MP.GAUSSLEGENDRE, MP.GAUSSJACOBI, MP.GAUSSLAGUERRE, MP.GAUSSCHEBYSHEVTYPE1, MP.GAUSSCHEBYSHEVTYPE2, MP.GAUSSGEGENBAUER
        
        %    References:
        %    J. Burkardt, Generalized Gauss-Hermite Quadrature Rules:
        %       http://people.sc.fsu.edu/~jburkardt/m_src/gen_hermite_rule/gen_hermite_rule.html
           
            if nargin < 4,          b     = mp('1');
                if nargin < 3,      a     = mp('0');
                    if nargin < 2,  alpha = mp('0');
                        if nargin < 1
                            error(  'MCT:GaussHermite:NotEnoughInputs',...
                                    'Not enough input arguments.' );
                        end  
                    end
                end
            end
            
            if alpha <= -1
                error(  'MCT:GaussHermite:UnsupportedAlpha',...
                        'Argument ALPHA must be > -1');
            end
            
            if ~ismp(n),     n = mp(n); end;            
            if ~ismp(a),     a = mp(a); end;
            if ~ismp(b),     b = mp(b); end;
            if ~ismp(alpha), alpha = mp(alpha); end;
            
            % Find maximum precision among all (manually supplied) arguments
            if     nargin == 1, precision = mp.Digits(n);
            elseif nargin == 2, precision = mp.Digits([n alpha]);
            elseif nargin == 3, precision = mp.Digits([n alpha a]);
            elseif nargin == 4, precision = mp.Digits([n alpha a b]);
            end;
            
            old_digits = mp.Digits();
            mp.Digits(precision);
            
            % Equalize precision of arguments - to make sure all computations 
            % are done with the requested precision (max among all arguments).
            %
            % Function is intended to work with variable precision arguments 
            % and there are a lot of tedious checks, etc. to make this work correctly. 
            % User doesn't need to take all these into account in his/her code.
            n     = mp(    n,precision); 
            alpha = mp(alpha,precision);  
            a     = mp(    a,precision); 
            b     = mp(    b,precision); 

            if (alpha == 0 && ((n <= 512 && precision == 34)||(n <= 256 && precision <= 250)))
                
                 % Use pre-computed values.
                 [x,w] = mpimpl(5001, n);
            else
                
                % In general case, there is no accurate approximation for roots of
                % Hermite polynomials, and unfortunately we cannot use
                % fast algorithm of initial root estimation + Newton refinement.
                %
                % Thus we have to rely on Golub-Welsch algorithm instead. 
                %
                % IMPORTANT NOTE. Golub-Welsch is supposed to be O(n^2),
                % thanks to specialized eigensolvers for symmetric
                % tridiagonal matrices (D&C, MRRR). However such solvers
                % suffer from accuracy loss on small components of
                % eigenvectors. That is because the solvers use Euclidean vector norm 
                % as stopping criteria which 'ignores' outliers - smallest
                % components of a vector. It would be much better to use
                % 'Inf' norm in order to guarantee uniform accuracy of
                % all components (see email exchange with Anna Matsekh).
                % Unfortunately this is not the case and we
                % have to use full blown QR solver of O(n^3) to avoid
                % accuracy loss.
                
                % Build tridiagonal symmetric Jacobi matrix associated 
                % with the orthogonal Hermite polynomials.
                aj = zeros(n,1,'mp');
                
                k  = 1:n-1;                
                bj = sqrt((k+alpha*mod(k,2))*0.5);
                J  = diag(aj,0) + diag(bj,1) + diag(bj,-1);

                % Find eigenvectors & eigenvalues of tridiagonal Jacobi matrix
                % MRRR and D&C tridiagonal eigensolvers might loose
                % accuracy on small components of eigenvectors. 
                %
                % Thus toolbox doesn't use these ultra-fast algorithms, 
                % relying on accurate but slow QR instead. 
                % Remove 'qr' if you want higher speed, but please note
                % that accuracy might drop.
                [V,D] = eig(J,'qr');
                [x,k] = sort(diag(D));
                
                mu = gamma ((alpha + 1)/2);                
                w  = mu * subsref(V,substruct('()',{1,k}))'.^2; % V(1,k)
            end;
            
            % Adjust nodes and weights for non-default arguments.
            if (b ~= 1 || a ~= 0 || alpha ~= 0)
               scale = 1/sqrt(b);
               shift = a;
               x = shift + scale.*x;
               w = w.*scale^(alpha + 1);
            end
            
            mp.Digits(old_digits);            
        end
        
        function [x,w] = GaussGegenbauer(n, alpha, a, b)
        %MP.GAUSSGEGENBAUER Compute coefficients of Gauss-Gegenbauer quadrature with arbitrary precision
        %
        %    [x,w] = mp.GaussGegenbauer(order, alpha, a, b) computes abscissas
        %    'x' and weights 'w' of Gauss-Gegenbauer quadrature.
        %
        %    The Gauss-Gegenbauer quadrature rule is used as follows:
        %
        %        Integral ( a <= x <= b ) f(x) * ( ( x - a ) * ( b - x ) )^alpha dx
        %
        %    is to be approximated by
        %
        %        sum ( 1 <= i <= n ) w(i) * f(x(i))        
        %
        %    Parameters:
        %    'n' is the number of points in the quadrature rule.
        %    'alpha' is the exponent, which must be greater than -1 (default 1).
        %    'a' is the left endpoint (default -1);
        %    'b' is the right endpoint(default +1);        
        %
        %    Precision is controlled by mp.Digits and by precision of the arguments
        %
        %    See also MP.GAUSSLEGENDRE, MP.GAUSSJACOBI, MP.GAUSSLAGUERRE, MP.GAUSSHERMITE, MP.GAUSSCHEBYSHEVTYPE1, MP.GAUSSCHEBYSHEVTYPE2
        
        %    References:
        %    J. Burkardt, Gauss-Gegenbauer Quadrature:
        %       http://people.sc.fsu.edu/~jburkardt/m_src/gegenbauer_rule/gegenbauer_rule.html

            if nargin < 4,          b     = mp('+1');
                if nargin < 3,      a     = mp('-1');
                    if nargin < 2,  alpha = mp(' 1');                    
                        if nargin < 1
                        error(  'MCT:GaussGegenbauer:NotEnoughInputs',...
                                'Not enough input arguments.' );
                        end  
                    end
                end
            end
            
            if alpha <= -1
                error(  'MCT:GaussGegenbauer:UnsupportedAlpha',...
                        'Argument ALPHA must be > -1');
            end
            
            if ~ismp(n),     n = mp(n); end;            
            if ~ismp(a),     a = mp(a); end;
            if ~ismp(b),     b = mp(b); end;
            if ~ismp(alpha), alpha = mp(alpha); end;
        
            % Find maximum precision among all (manually supplied) arguments
            if     nargin == 1, precision = mp.Digits(n);
            elseif nargin == 2, precision = mp.Digits([n alpha]);
            elseif nargin == 3, precision = mp.Digits([n alpha a]);
            elseif nargin == 4, precision = mp.Digits([n alpha a b]);
            end;
            
            old_digits = mp.Digits();
            mp.Digits(precision);
            
            % Equalize precision of arguments - to make sure all computations 
            % are done with the requested precision (max among all arguments).
            n     = mp(    n,precision); 
            alpha = mp(alpha,precision);  
            a     = mp(    a,precision); 
            b     = mp(    b,precision); 
            
            if (alpha == 1 && ((n <= 512 && precision == 34)||(n <= 256 && precision <= 250)))
                
                 % Use pre-computed values.
                 [x,w] = mpimpl(5002, n);
            else
                
                % Golub-Welsch algorithm (see comments in mp.GaussHermite).
                aj = zeros(n,1,'mp');
                
                k  = 1:n-1;
                bj = sqrt(k.*(k+2*alpha)./(4*(k+alpha).^2-1));
                J  = diag(aj,0) + diag(bj,1) + diag(bj,-1);
                
                [V,D] = eig(J,'qr');
                [x,k] = sort(diag(D));
                
                mu = 2^(2*alpha + 1) * gamma(alpha + 1)^2 / gamma(2*alpha + 2);
                w  = mu * subsref(V,substruct('()',{1,k}))'.^2; % V(1,k)
            end;
            
            % Adjust nodes and weights for non-default arguments.
            if (b ~= 1 || a ~= -1 || alpha ~= 1)
               scale = (b-a)/2;
               shift = (a+b)/2;
               x = shift + scale.*x;
               w = w.*scale^(2*alpha+1);
            end

            mp.Digits(old_digits);                        
        end
        
        function [x,w] = GaussJacobi(n, alpha, beta, a, b)
        %MP.GAUSSJACOBI Compute coefficients of Gauss-Jacobi quadrature with arbitrary precision
        %
        %    [x,w] = mp.GaussJacobi(order, alpha, beta, a, b) computes abscissas
        %    'x' and weights 'w' of Gauss-Jacobi quadrature.
        %
        %    The Gauss-Jacobi quadrature rule is used as follows:
        %
        %        Integral ( a <= x <= b ) (b-x)^alpha (x-a)^beta f(x) dx
        %
        %    is to be approximated by
        %
        %        sum ( 1 <= i <= n ) w(i) * f(x(i))        
        %
        %    Parameters:
        %    'n' is the number of points in the quadrature rule.
        %    'alpha' is the value of the exponent of (b-x), which must be greater than -1 (default 1).
        %    'beta' is the value of the exponent of (x-a), which must be greater than -1 (default 0).
        %    'a' is the left endpoint (default -1);
        %    'b' is the right endpoint(default +1);        
        %
        %    Precision is controlled by mp.Digits and by precision of the arguments
        %
        %    See also MP.GAUSSLEGENDRE, MP.GAUSSLAGUERRE, MP.GAUSSHERMITE, MP.GAUSSCHEBYSHEVTYPE1, MP.GAUSSCHEBYSHEVTYPE2, MP.GAUSSGEGENBAUER
        
        %    References:
        %    J. Burkardt, Gauss-Jacobi Quadrature:
        %       http://people.sc.fsu.edu/~jburkardt/m_src/jacobi_rule/jacobi_rule.html

            if nargin < 5,              b     = mp('+1');
                if nargin < 4,          a     = mp('-1');
                    if nargin < 3,      beta  = mp('0');                    
                        if nargin < 2,  alpha = mp('1');                    
                            if nargin < 1
                            error(  'MCT:GaussJacobi:NotEnoughInputs',...
                                    'Not enough input arguments.' );
                            end  
                        end
                    end
                end
            end
            
            if alpha <= -1
                error(  'MCT:GaussJacobi:UnsupportedAlpha',...
                        'Argument ALPHA must be > -1');
            end
            
            if beta <= -1
                error(  'MCT:GaussJacobi:UnsupportedBeta',...
                        'Argument BETA must be > -1');
            end
            
            if ~ismp(n),     n = mp(n); end;            
            if ~ismp(a),     a = mp(a); end;
            if ~ismp(b),     b = mp(b); end;
            if ~ismp(alpha), alpha = mp(alpha); end;
            if ~ismp(beta),  beta  = mp(beta); end;    
        
            % Find maximum precision among all (manually supplied) arguments
            if     nargin == 1, precision = mp.Digits(n);
            elseif nargin == 2, precision = mp.Digits([n alpha]);
            elseif nargin == 3, precision = mp.Digits([n alpha beta]);
            elseif nargin == 4, precision = mp.Digits([n alpha beta a]);
            elseif nargin == 5, precision = mp.Digits([n alpha beta a b]);                
            end;
            
            old_digits = mp.Digits();
            mp.Digits(precision);
            
            n     = mp(    n,precision); 
            alpha = mp(alpha,precision);  
            beta  = mp( beta,precision);              
            a     = mp(    a,precision); 
            b     = mp(    b,precision); 
            
            if (alpha == 1 && beta == 0 && ((n <= 512 && precision == 34)||(n <= 256 && precision <= 250)))
                
                 % Use pre-computed values.
                 [x,w] = mpimpl(5003, n);
            else
                
                % Golub-Welsch algorithm (see comments in mp.GaussHermite).
                ab  = alpha + beta;
                abi = 2*(1:n) + ab;
                aj  = (beta^2 - alpha^2)./((abi-2).*abi);

                k   = 1:n-1;                
                abi = (2*k + ab).^2;
                bj  = sqrt(4*k.*(k + alpha).*(k + beta).*(k + ab)./((abi-1).*abi));

                J  = diag(aj,0) + diag(bj,1) + diag(bj,-1);

                [V,D] = eig(J,'qr');
                [x,k] = sort(diag(D));
                
                mu = 2^(ab + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(2 + ab);
                w  = mu * subsref(V,substruct('()',{1,k}))'.^2; % V(1,k)
            end;
            
            % Adjust nodes and weights for non-default arguments.
            if (b ~= 1 || a ~= -1 || alpha ~= 1 || beta ~= 0)
               scale = (b-a)/2;
               shift = (a+b)/2;
               x = shift + scale * x;
               w = w.*scale^(alpha+beta+1);
            end
            
            mp.Digits(old_digits);            
        end
        
        function [x,w] = GaussLaguerre(n, alpha, a, b)
        %MP.GAUSSLAGUERRE Compute coefficients of Generalized Gauss-Laguerre quadrature with arbitrary precision
        %
        %    [x,w] = mp.GaussLaguerre(order, alpha, a, b) computes abscissas
        %    'x' and weights 'w' of Generalized Gauss-Laguerre quadrature.
        %
        %    The generalized Gauss-Laguerre quadrature rule is used as follows:
        %
        %        Integral ( a < x < +inf ) (x-a)^alpha * exp(-b*(x - a)) f(x) dx
        %      
        %    is to be approximated using 'x' and 'w' by
        %
        %         sum ( 1 <= i <= n ) w(i) * f(x(i))      
        %
        %    Parameters:
        %    'n' is the number of points in the quadrature rule.
        %    'alpha' is the exponent of |x| in the weight function. 
        %           The value of alpha may be any real value greater than -1.0. 
        %           Specifying alpha = 0.0 (default) results in the basic (non-generalized) rule.
        %    'a' is the left point (default 0);
        %    'b' is the scale factor in the exponential (default 1);        
        %
        %    Precision is controlled by mp.Digits and by precision of the arguments
        %
        %    See also MP.GAUSSLEGENDRE, MP.GAUSSJACOBI, MP.GAUSSHERMITE, MP.GAUSSCHEBYSHEVTYPE1, MP.GAUSSCHEBYSHEVTYPE2, MP.GAUSSGEGENBAUER       
        
        %    References:
        %    J. Burkardt, Generalized Gauss-Laguerre Quadrature Rules:
        %       http://people.sc.fsu.edu/~jburkardt/m_src/gen_laguerre_rule/gen_laguerre_rule.html
        
            if nargin < 4,          b     = mp('1');
                if nargin < 3,      a     = mp('0');
                    if nargin < 2,  alpha = mp('0');
                        if nargin < 1
                        error(  'MCT:GaussLaguerre:NotEnoughInputs',...
                                'Not enough input arguments.' );
                        end  
                    end
                end
            end
            
            if alpha <= -1
                error(  'MCT:GaussLaguerre:UnsupportedAlpha',...
                        'Argument ALPHA must be > -1');
            end
            
            if ~ismp(n),     n = mp(n); end;            
            if ~ismp(a),     a = mp(a); end;
            if ~ismp(b),     b = mp(b); end;
            if ~ismp(alpha), alpha = mp(alpha); end;
        
            % Find maximum precision among all (manually supplied) arguments
            if     nargin == 1, precision = mp.Digits(n);
            elseif nargin == 2, precision = mp.Digits([n alpha]);
            elseif nargin == 3, precision = mp.Digits([n alpha a]);
            elseif nargin == 4, precision = mp.Digits([n alpha a b]);
            end;
            
            old_digits = mp.Digits();
            mp.Digits(precision);
       
            % Equalize precision of arguments - to make sure all computations 
            % are done with the requested precision (max among all arguments).
            n     = mp(    n,precision); 
            alpha = mp(alpha,precision);  
            a     = mp(    a,precision); 
            b     = mp(    b,precision); 
            
            if (alpha == 0 && ((n <= 512 && precision == 34)||(n <= 256 && precision <= 250)))
                
                 % Use pre-computed values in case of quadruple precision. 
                 [x,w] = mpimpl(5004, n);
            else
                
                % Golub-Welsch algorithm (see comments in mp.GaussHermite).
                aj = 2*(1:n)-1 + alpha;
                k  = 1:n-1;
                bj = sqrt(k.*(k + alpha));  
                J  = diag(aj,0) + diag(bj,1) + diag(bj,-1);
                
                [V,D] = eig(J,'qr');
                [x,k] = sort(diag(D));
                
                mu = gamma(alpha + 1.0);
                w  = mu * subsref(V,substruct('()',{1,k}))'.^2; % V(1,k)
            end;
            
            % Adjust nodes and weights for non-default arguments.
            if (b ~= 1 || a ~= 0 || alpha ~= 0)
               scale = 1/b;
               shift = a;
               x = shift + scale.*x;
               w = w.*scale^(alpha+1);
            end
            
            mp.Digits(old_digits);
        end
        
        function [x,w] = GaussChebyshevType1(n, a, b)
        %MP.GAUSSCHEBYSHEVTYPE1 Compute coefficients of Gauss-Chebyshev Type 1 quadrature with arbitrary precision 
        %
        %    [x,w] = mp.GaussChebyshevType1(order, a, b) computes abscissas
        %    'x' and weights 'w' of Gauss-Chebyshev Type 1 quadrature.
        %
        %    The Gauss-Chebyshev Type 1 quadrature rule is used as follows:
        %
        %        Integral ( a <= x <= b ) f(x) / sqrt ( ( x - a ) * ( b - x ) ) dx
        %
        %    is to be approximated by
        %
        %        sum ( 1 <= i <= order ) w(i) * f(x(i))        
        %
        %    Parameters:
        %    'n' is the number of points in the quadrature rule.
        %    'a' is the left endpoint (default -1);
        %    'b' is the right endpoint(default +1);        
        %
        %    Precision is controlled by mp.Digits and by precision of the arguments
        %
        %    See also MP.GAUSSLEGENDRE, MP.GAUSSJACOBI, MP.GAUSSLAGUERRE, MP.GAUSSHERMITE, MP.GAUSSCHEBYSHEVTYPE2, MP.GAUSSGEGENBAUER
        
        %    References:
        %    J. Burkardt, Gauss-Chebyshev Type 1 Quadrature:
        %       http://people.sc.fsu.edu/~jburkardt/m_src/chebyshev1_rule/chebyshev1_rule.html
        
            if nargin < 3,          b     = mp('+1');
                if nargin < 2,      a     = mp('-1');
                    if nargin < 1
                    error(  'MCT:GaussChebyshevType1:NotEnoughInputs',...
                            'Not enough input arguments.' );
                    end  
                end
            end
            
            if ~ismp(n), n = mp(n); end;            
            if ~ismp(a), a = mp(a); end;
            if ~ismp(b), b = mp(b); end;
        
            % Find maximum precision among all (manually supplied) arguments
            if     nargin == 1, precision = mp.Digits(n);
            elseif nargin == 2, precision = mp.Digits([n a]);
            elseif nargin == 3, precision = mp.Digits([n a b]);
            end;
            
            old_digits = mp.Digits();
            mp.Digits(precision);

            x = cos((2*(n:-1:1)-1)*mp('pi')/(2*n))';
            w = repmat(mp('pi')/n,n,1);

            % Adjust nodes and weights for non-default arguments.
            if (b ~= 1 || a ~= -1)
               scale = (b-a)/2;
               shift = (a+b)/2;
               x = shift + scale.*x;
            end
            
            mp.Digits(old_digits);            
        end

        function [x,w] = GaussChebyshevType2(n, a, b)
        %MP.GAUSSCHEBYSHEVTYPE1 Compute coefficients of Gauss-Chebyshev Type 2 quadrature with arbitrary precision
        %
        %    [x,w] = mp.GaussChebyshevType2(order, a, b) computes abscissas
        %    'x' and weights 'w' of Gauss-Chebyshev Type 2 quadrature.
        %
        %    The Gauss-Chebyshev Type 2 quadrature rule is used as follows:
        %
        %        Integral ( a <= x <= b ) f(x) * sqrt ( ( x - a ) * ( b - x ) ) dx
        %
        %    is to be approximated by
        %
        %        sum ( 1 <= i <= order ) w(i) * f(x(i))        
        %
        %    Parameters:
        %    'n' is the number of points in the quadrature rule.
        %    'a' is the left endpoint (default -1);
        %    'b' is the right endpoint(default +1);        
        %
        %    Precision is controlled by mp.Digits and by precision of the arguments
        %
        %    See also MP.GAUSSLEGENDRE, MP.GAUSSJACOBI, MP.GAUSSLAGUERRE, MP.GAUSSHERMITE, MP.GAUSSCHEBYSHEVTYPE1, MP.GAUSSGEGENBAUER       
        
        %    References:
        %    J. Burkardt, Gauss-Chebyshev Type 2 Quadrature:
        %       http://people.sc.fsu.edu/~jburkardt/m_src/chebyshev1_rule/chebyshev2_rule.html
        
            if nargin < 3,          b     = mp('+1');
                if nargin < 2,      a     = mp('-1');
                    if nargin < 1
                    error(  'MCT:GaussChebyshevType2:NotEnoughInputs',...
                            'Not enough input arguments.' );
                    end  
                end
            end
            
            if ~ismp(n), n = mp(n); end;            
            if ~ismp(a), a = mp(a); end;
            if ~ismp(b), b = mp(b); end;

            % Find maximum precision among all (manually supplied) arguments
            if     nargin == 1, precision = mp.Digits(n);
            elseif nargin == 2, precision = mp.Digits([n a]);
            elseif nargin == 3, precision = mp.Digits([n a b]);
            end;
            
            old_digits = mp.Digits();
            mp.Digits(precision);

            x = cos((n:-1:1)*mp('pi')/(n+1))';
            w = mp('pi')/(n+1)*(sin((n:-1:1)*mp('pi')/(n+1)).^2)';

            % Adjust nodes and weights for non-default arguments.
            if (b ~= 1 || a ~= -1)
               scale = (b-a)/2;
               shift = (a+b)/2;
               x = shift + scale.*x;
               w = w.*scale^2;
            end
            
            mp.Digits(old_digits);            
        end
       
        function [x,w] = GaussKronrod(n)
        %MP.GAUSSKRONROD Compute coefficients of Gauss-Kronrod quadrature with arbitrary precision. 
        %
        %    [x,w] = mp.GaussKronrod(n) computes abscissae
        %    'x' and weights 'w' of Gauss-Kronrod quadrature for [-1,1]
        %
        %    Parameters:
        %    'n' is the number of points in the quadrature rule.
        %
        %    Precision is controlled by mp.Digits and by precision of the arguments
        %
        %    See also MP.GAUSSLEGENDRE, MP.GAUSSJACOBI, MP.GAUSSLAGUERRE, MP.GAUSSHERMITE, MP.GAUSSCHEBYSHEVTYPE1, MP.GAUSSCHEBYSHEVTYPE2, MP.GAUSSGEGENBAUER
        
            if nargin < 1
                error(  'MCT:GaussKronrod:NotEnoughInputs',...
                        'Not enough input arguments.' );
            end  
            
            if ~ismp(n), n = mp(n); end;
            precision = mp.Digits(n);

            old_digits = mp.Digits();
            mp.Digits(precision);

            ab = r_jacobi01(2*n,mp(0),mp(0));
            t  = kronrod(n,ab);
            x  = 2*subsref(t,substruct('()',{':',1}))-1; %2*t(:,1)-1;
            w  = 2*subsref(t,substruct('()',{':',2}));   %2*t(:,2);

            mp.Digits(old_digits);            
        end
        
        function varargout = BernoulliNumber(n)
        %MP.BERNOULLINUMBER Compute Bernoulli numbers with arbitrary precision.   
        %
        %    B = mp.BernoulliNumber(n) computes n-th Bernoulli number(s) with
        %    any required accuracy. n can be a vector.
        %
        %    Precision is controlled by MP.DIGITS
        %
        %    See also MP.DIGITS        
        
            [varargout{1:nargout}] = mpimpl(417, n);
        end

        
        %% File I/O in textual format
        function write(A,filename)
        %WRITE Writes multiprecision matrix A to file in textual format
        %
        %    Stores values with the precision enough for exact
        %    restoration on read.
        %
        %    mp.Digits(34);
        %    
        %    A = mp(rand(25));
        %    
        %    mp.write(A,'matrix.txt');  % Write mp-matrix to the text file
        %    B = mp.read('matrix.txt'); % Read it back
        %    
        %    norm(A-B,1) % check accuracy - difference should be 0
        %
        %    See also MP.READ.

            if nargin < 2, error(message('MATLAB:narginchk:notEnoughInputs'));  end;
            if ~ismp(A), error('A must be of mp-type'); end;
            if ~ischar(filename), error('Second parameter must be a file name stored as string'); end;
            
            s = num2str(A,'%e');
            fout = fopen(filename,'w');
            fwrite(fout,s);
            fclose(fout);
        end

        function A = read(filename)
        %READ Reads multiprecision matrix from text file
        %
        %    Loads matrix from text file saved by WRITE.
        %    Matrix elements are automatically rounded to current
        %    precision mp.Digits().
        %
        %    See also MP.WRITE, MP.DIGITS.
        
            if ~ischar(filename), error('File name must be a string'); end;        
            
            s = fileread(filename);
            A = mp(['[',s,']']);
        end
        
        function varargout = cellfun(func, varargin)
        %CELLFUN Apply a function to each cell of a cell array.
        %
        %   Please be aware, this is early alpha implementation, use it with
        %   caution. Report issues/suggest improvements to support@advanpix.com
        %   
        %   See also  ARRAYFUN, STRUCTFUN, FUNCTION_HANDLE, CELL2MAT, SPFUN

          narginchk(2,inf);

          % Built-in cellfun doesn't support UniformOutput = true for custom
          % objects, so we have to disable it for a moment.
          isUniformOutput = true;
          for k=1:numel(varargin)
            if ischar(varargin{k}) && strcmp(varargin{k},'UniformOutput')
                if(k+1 <= numel(varargin))
                    isUniformOutput = logical(varargin{k+1});
                    varargin{k+1}   = false;
                end
                break;
            end
          end

          if isUniformOutput
              [varargout{1:nargout}] = builtin('cellfun', func, varargin{:}, 'UniformOutput', false);                        
              for k=1:nargout
                  A = varargout{k};
                  varargout{k} = [A{:}];  % this simple trick should work for uniform outputs
              end
          else
             [varargout{1:nargout}] = builtin('cellfun', func, varargin{:});
          end
        end        
        
        function [h,l] = getwords64(x)
        %GETWORDS64 Splits quadruple precision input into two arrays of unsigned 64-bit integers.
        %
        %    The function provides low-level access to underlying data of quadruple 
        %    precision numbers - designed to be used by developers.
        %
        %    Quadruple floating-point number occupies 128 bit in memory.
        %    Thus it can be viewed as two unsigned 64-bit integers.
        %    Function extracts these two values.  
        %    
        %    Extracted integers can be stored, transferred to
        %    another system and then quadruple floats can be
        %    reconstructed by SETWORDS64.
        %
        %    See also SETWORDS64.
        
           [h,l] = mpimpl(9050,x); 
        end

        function r = setwords64(h,l)
        %SETWORDS64 Generates quadruple precision array from two 
        %           arrays of unsigned 64-bit integers.
        %
        %    The function provides low-level access to underlying data of quadruple 
        %    precision numbers - designed to be used by developers.
        %
        %    Quadruple floating-point number occupies 128 bit in memory.
        %    Thus it can be viewed as two unsigned 64-bit integers.
        %
        %    Function builds quadruple value/array from 64-bit integer values/array.
        %
        %    See also GETWORDS64.
        
             r = mpimpl(9055,h,l); 
        end
        
    end  % Static
    
    methods(Static=true, Hidden=true)
        
        %% Untilities for internal purposes - do not call them directly
        function Init()
        %MP.INIT  Initializes Multiprecision Computing Toolbox.  
        %
        %    MP.INIT Loads and configures Multiprecision Computing Toolbox's core engine  
        %    Called automatically on first usage of any of toolbox functions.
        % 
        %    Don't need to be called explicitly by user.         
        %
        %    See also MPSTARTUP
        
            if ispc, fs = '\';
            else     fs = '/';  end %#ok<SEPEX>
            mctpath = mfilename('fullpath');
            temp = strfind( mctpath, fs );
            mctpath( temp(end) : end ) = [];
            
            mpimpl(10000, [mctpath, fs]); 
        end
        
        function r = convertListToDouble(varargin)
          for i=1:nargin, if isa(varargin{i},'mp'), varargin{i} = double(varargin{i}); end; end;
          r = varargin;
        end
         
        function r = genericArrayCreationUsingBuiltIn(funcname, varargin)
        % Basic array-creation functions (zeros, ones, eye, etc.) follow
        % the same syntax and we can handle them in generic way.
        %
        % Some information on array-creation functions for OOP:
        % http://www.mathworks.com/help/matlab/matlab_oop/class-support-for-array-creation-functions.html

         % Always create multiprecision array when 'mp' is required explicitely:
         % func(...,'mp')
         if ischar(varargin{end}) && strcmpi(varargin{end},'mp')
             args = mp.convertListToDouble(varargin{1:end-1});
             r = mp(builtin(funcname,args{:}));             
         else
            % We use aggressive strategy and override 'double'/'single'.
            if ischar(varargin{end}) && ~(strcmpi(varargin{end},'double') || strcmpi(varargin{end},'single')) || ...
                    (nargin > 2 && ischar(varargin{end-1}) && strcmpi(varargin{end-1},'like') ...
                    && ~(isa(varargin{end},'mp') || isa(varargin{end},'double') || isa(varargin{end},'single')))

                % Call built-in only when non-'mp'/'double'/'single' output is required:
                % X = func(___,typename), where typename ~= 'mp'/'double'/'single'    
                % X = func(___,'like',p), where class(p) ~= 'mp'/'double'/'single'
                args = mp.convertListToDouble(varargin{:});
                r = builtin(funcname,args{:});
            else

                % Create multiprecision array in all other cases:
                % X = func(n)
                % X = func(sz1,...,szN)
                % X = func(sz)
                % X = func(___,typename), where typename == 'mp'/'double'/'single'    
                % X = func(___,'like',p), where class(p) == 'mp'/'double'/'single'
                args = mp.convertListToDouble(varargin{:});
                r = mp(builtin(funcname,args{:}));
            end
         end
        end
        
        function [matchstring, tokenname, splitstring] = parseSpecifiers(formatSpec)
        % PARSESPECIFIERS Parses format specifiers in the format string.
        %
        %    Parses only valid format specifiers (with valid types, etc.)
        % 
        %    Do not call it explicitely.
        %    Needed for implementation of full-featured SPRINTF/NUM2STR.
        %
        %    See also SPRINTF, NUM2STR.
        %
        %    I. DevNotes for the future:
        %    This regexp covers extended printf formatting including all the
        %    MATLAB's extensions like 'identifier', e.g. %1$*5$.*2$, etc.:
        %    
        %    (?<!%)(?:%%)*%(?!%)(?:(\d+)\$)?([\+-\s#]?)(?:(?:\*+(\d+)\$)|(\d*))?(\.)?(?:(?:\*+(\d+)\$)|(\d*))?([diuoxXfeEgGcs]|(?:bx|bX|bo|bu|tx|tX|to|tu))
        %     
        %    By fields:
        %    (?<!%)(?:%%)*%(?!%)
        %    (?:(\d+)\$)?
        %    ([\+-\s#]?)
        %    (?:(?:\*+(\d+)\$)|(\d*))?
        %    (\.)?
        %    (?:(?:\*+(\d+)\$)|(\d*))?
        %    ([diuoxXfeEgGcs]|(?:bx|bX|bo|bu|tx|tX|to|tu))
        %    
        %    Explanations on capturing groups:
        %    group 1 (id)      : (\d+)
        %    group 2 (flag)    : ([\+-\s#]?)
        %    group 3 (width_id): (\d+)
        %    group 4 (width)   : (\d*)
        %    group 5           : (\.)
        %    group 6 (prec_id) : (\d+)
        %    group 7 (prec)    : (\d*)
        %    group 8 (type)    : ([diuoxXfeEgGcs]|(?:bx|bX|bo|bu|tx|tX|to|tu))
        %     
        %    However MATLAB has handicapped regexp support. Nested groups (tokens) are not supported
        %    This makes impossible to use short and elagant regexps as above :(.
        %    Maybe we will implement this on C++ level in future.
        %
        %    II. Current Implementation
        %    Due to limitations of MATLAB, now we support only standard
        %    printf specification:
        %    %[flag][width][.precision]type
        %    flag: +, -, space, # 
        %    width, precision - numbers
        %    Example: %e, %.15f, etc
        %    (?<!%)(?:%%)*%(?!%)([\+-\s#]?)(\d*)(\.)?(\d*)([diuoxXfeEgGcs]|(?:bx|bX|bo|bu|tx|tX|to|tu))
        %    
        %    Explanations on capturing groups:            
        %    group 1 (flag) : ([\+-\s#]?)
        %    group 2 (width): (\d*)
        %    group 3        : (\.)
        %    group 4 (prec) : (\d*)
        %    group 5 (type) : ([diuoxXfeEgGcs]|(?:bx|bX|bo|bu|tx|tX|to|tu))
        
            regexp_start  = '(?<!%)(?:%%)*%(?!%)';  
            regexp_flag   = '(?<flag>[\+-\s#]?)';
            regexp_width  = '(?<width>\d*)';
            regexp_dot    = '(?<dot>\.?)';
            regexp_prec   = '(?<prec>\d*)';
            regexp_type   = '(?<type>[diuoxXfeEgGcs]|(?:bx|bX|bo|bu|tx|tX|to|tu|ld|li|lo|lu|lx|lX|hd|hi|ho|hu|hx|hX))';
            
            expr = [regexp_start, regexp_flag, regexp_width, regexp_dot, regexp_prec, regexp_type];
            [matchstring, tokenname, splitstring] = regexp(formatSpec, expr, 'match', 'names', 'split');
        end
        
        function varargout = sprintf_basic(varargin)
        % SPRINTF_BASIC Basic implementation of sprintf.
        %
        %    Do not call it explicitely.
        %    Needed for implementation of full-featured SPRINTF/NUM2STR.
        %
        %    See also SPRINTF, NUM2STR.
        
            [matchstring, tokenname, splitstring] = mp.parseSpecifiers(varargin{1});
         
            % Process every format, find corresponding argument, convert to
            % string if argument of 'mp' type. Replace corresponding format
            % to '%s'
            offs = 1;            
            varargin{offs} = '';  % New format string
            for i=1:min(size(matchstring, 2),nargin-1)
                    
                if isa(varargin{offs + i}, 'mp')
                    
                   % Process 'mp' objects separately
                   
                   % 1. Analyse [flag][width][precision]type fields in future
                   % using tokenname array of structures
                   
                   % 1.1 Use default format '%e' for unsupported specifiers
                   unsupported_types = '([diuoxXcs]|(?:bx|bX|bo|bu|tx|tX|to|tu|ld|li|lo|lu|lx|lX|hd|hi|ho|hu|hx|hX))';
                   
                   index = regexp(tokenname(i).type, unsupported_types);
                   if( size(index) ~= 0)
                       
                       % In case of %i or %d we use %e for numbers with fractional part
                       % Otherwise (if fractional part exists) we use %g which gives
                       % correct result for integers, as required.
                       x = varargin{offs + i};
                       if ~isequal(x, fix(x)), modifiedFormat = regexprep(matchstring{i},unsupported_types,'e');
                       else                    modifiedFormat = regexprep(matchstring{i},unsupported_types,'g'); end %#ok<SEPEX>
                       
                       varargin{offs + i} = num2str(varargin{offs + i}, modifiedFormat);                       
                   else
                       varargin{offs + i} = num2str(varargin{offs + i}, matchstring{i});   
                   end
                  
                   % Use '%s' for converted 'mp' objects
                   varargin{offs} = [varargin{offs}, splitstring{i}, '%s'];
                else
                    
                   % Leave other arguments unchanged
                   varargin{offs} = [varargin{offs}, splitstring{i}, matchstring{i}];
                end
            end
            
            % Add last split & re-write format string to a new one
            varargin{offs} = [varargin{offs}, splitstring{end}];
            
            % Call standard fprintf
            [varargout{1:nargout}] = sprintf(varargin{:});
        end
    end
        
    % class operators and functions
    methods
        
        %% Constructor
        function this = mp(x, precision)
        %MP Create multiprecision floating-point entity (number, matrix, n-dim array) of required precision.
        %
        %    mp()    creates multiprecision number of default precision (with zero value).
        %    mp(x)   creates multiprecision entity from numeric array, matrix, expression or other mp-object.
        % 
        %    Please check our User's Manual for more detailed information:
        %    http://www.advanpix.com/documentation/users-manual/
        %
        %    Examples:
        %
        %    mp.Digits(34);   % Set default precision to 34 digits (quadruple precision)
        %    
        %    x = mp('pi/4')   % Simple expression evaluation 
        %    x = 
        %       0.7853981633974483096156608458198757
        
        %    y = mp('sqrt(2)/2')
        %    y = 
        %       0.707106781186547524400844362104849        
        %
        %    % Assemble matrix row by row:
        %    a1 = mp('[ 5/36              2/9-sqrt(15)/15   5/36-sqrt(15)/30 ]');
        %    a2 = mp('[ 5/36+sqrt(15)/24  2/9               5/36-sqrt(15)/24 ]');
        %    a3 = mp('[ 5/36+sqrt(15)/30  2/9+sqrt(15)/15   5/36             ]');
        %    a = [a1;a2;a3];                  % Concatenate rows into final matrix         
        % 
        %    A = mp(magic(5)); % converts magic square to floating-point matrix of 34-digits precision.
        %    
        %    B = randn([3,3,3,3],'mp') % creates 4D array of normally distributed pseudo-random numbers.
        %
        %    ***
        %    Once created multiprecision objects can be used in calculations exactly the same way 
        %    as built-in 'double' or 'single' precision entities:
        %
        %    [U,S,V] = svd(A);
        %    norm(A-U*S*V',1)  
        %    ans = 
        %          6.08593862426366529565689029856837e-32        
        %
        %    Please check our User's Manual for more examples and details:
        %    http://www.advanpix.com/documentation/users-manual/
        %
        %    See also MP.DIGITS        

           if nargin == 0,     this.id = mpimpl(0);             % create mp-entity of default precision
           elseif nargin == 1, this    = mpimpl(0,x);           % create mp-entity from numeric array, expression or other mp object.     
           else                this    = mpimpl(0,x,precision); %#ok<SEPEX> % create mp-entity of prescribed precision
           end 
        end % constructor

        %% Machine epsilon
        function r = eps(varargin)
        %EPS Spacing of multiprecision floating point numbers (depends on precision).
        % 
        %    EPS(x) returns positive distance from abs(x) to the next larger in magnitude multiprecision 
        %    floating point number of the same precision as x.
        %
        %    Precision is controlled by mp.Digits.
        %
        %    Usage of EPS is completely analogous to usage of built-in function EPS.
        %
        %    See also mp.Digits, REALMAX, REALMIN
      
            if nargin == 0
                r = mpimpl(2000); % eps()
            elseif nargin == 1 && strcmpi('mp',class(varargin{1}))
                r = mpimpl(2000,varargin{1});  % eps(x)
            elseif strcmpi('mp',varargin{end}) || strcmpi('mp',class(varargin{end}))
                r = mpimpl(2000); % eps('mp') or eps('like',p)
            else
                r = builtin('eps',x);
            end
        end
        
        %% Indexing operations
        function ind = subsindex(this)
        %SUBSINDEX Subscript indexing with object
        % 
        %  Convert the mp-object a to double format to be used
        %  as an index in an indexing expressions.
        
           ind = double(this)-1;
        end        
    
        function varargout = subsref(varargin)
        %SUBSREF Subscripted reference, B = A(S), or B = subsref(A,S)
        %
        %    Usage is identical to built-in function SUBSREF.
        %
        %    See also SUBSREF.
        
             [varargout{1:nargout}] = mpimpl(170, varargin{:});               
        end

        function varargout = subsasgn(varargin)
        %SUBSASGN Subscripted assignment, A(S) = B, or A = subsasgn(A,S,B)
        %
        %    Usage is identical to built-in function SUBSASGN.
        %
        %    See also SUBSASGN.

             if ~isa(varargin{1}, 'mp'), varargin{1} = mp(varargin{1}); end;                     
             if ~isa(varargin{3}, 'mp'), varargin{3} = mp(varargin{3}); end;
             [varargout{1:nargout}] = mpimpl(171, varargin{:});
        end
         
        function varargout = size(varargin)
        %SIZE Array dimensions
        %
        %    Usage is identical to built-in function SIZE.
        %
        %    See also SIZE.

            [varargout{1:nargout}] = mpimpl(172, varargin{:});  
        end
        
        function out = length(obj)
        %LENGTH Length of vector or largest array dimension
        %
        %    Usage is identical to built-in function LENGTH.
        %
        %    See also LENGTH.

            if isempty(obj), out = 0; 
            else             out = max(size(obj)); end; %#ok<SEPEX>
        end

        function varargout = end(varargin)
        %END Last index of array
        %
        %    Usage is identical to built-in function END.
        %
        %    See also END.

            [varargout{1:nargout}] = mpimpl(175, varargin{:});              
        end
        
        function varargout = numel(varargin)
        %NUMEL Number of array elements
        %
        %    Usage is identical to built-in function NUMEL.
        %
        %    See also NUMEL.
                        
            A = varargin{1};            
            if nargin > 1
                S = struct('type','()','subs',[]);
                S.subs = varargin(2:nargin);
                [varargout{1:nargout}] = mpimpl(173, A, S);                  
            else
                [varargout{1:nargout}] = mpimpl(173, A);                                  
            end    
        end
         
        function varargout = ndims(varargin)
        %NDIMS Number of array dimensions
        %
        %    Usage is identical to built-in function NDIMS.
        %
        %    See also NDIMS.
        
             [varargout{1:nargout}] = mpimpl(174, varargin{:});
        end

        function varargout = cat(varargin)
        %CAT Concatenate arrays along specified dimension
        %
        %    Usage is identical to built-in function CAT.
        %
        %    See also CAT.
        
            for i=2:nargin, if ~isa(varargin{i}, 'mp'), varargin{i} = mp(varargin{i}); end; end
            [varargout{1:nargout}] = mpimpl(177, varargin{:});
        end

        function varargout = reshape(varargin)
        %RESHAPE Reshape array
        %
        %    Usage is identical to built-in function RESHAPE.
        %
        %    See also RESHAPE.
        
            [varargout{1:nargout}] = mpimpl(178, varargin{:});
        end

        function varargout = vertcat(varargin)
        %VERTCAT Concatenate arrays vertically
        %
        %    Usage is identical to built-in function VERTCAT.
        %
        %    See also VERTCAT.
        
            [varargout{1:nargout}] = cat(1,varargin{:});
        end

        function varargout = horzcat(varargin)
        %HORZCAT Concatenate arrays horizontally
        %
        %    Usage is identical to built-in function HORZCAT.
        %
        %    See also HORZCAT.
       
            [varargout{1:nargout}] = cat(2,varargin{:});
        end
        
        function varargout = permute(varargin)
        %PERMUTE Rearrange dimensions of N-D array
        %
        %    Usage is identical to built-in function PERMUTE.
        %
        %    See also PERMUTE.
        
             [varargout{1:nargout}] = mpimpl(179, varargin{:});               
        end

        function a = ipermute(b,order)
        %IPERMUTE Inverse permute dimensions of N-D array
        %
        %    Usage is identical to built-in function IPERMUTE.
        %
        %    See also IPERMUTE.
        
            inverseorder(order) = 1:numel(order);
            a = permute(b,inverseorder);
        end
        
        function varargout = bsxfun(varargin)
        %BSXFUN  Binary Singleton Expansion Function
        %
        %    Usage is identical to built-in function BSXFUN.
        %
        %    See also BSXFUN.

            [varargout{1:nargout}] = mpbsxfun(varargin{:});            
        end
        
        function varargout = rot90(varargin)
        %ROT90 Rotate array 90 degrees.
        %
        %    Usage is identical to built-in function ROT90.
        %
        %    See also ROT90.
        
             [varargout{1:nargout}] = mprot90(varargin{:});               
        end
        
        function varargout = squeeze(varargin)
        %SQUEEZE Remove singleton dimensions.
        %
        %    Usage is identical to built-in function SQUEEZE.
        %
        %    See also SQUEEZE.
        
             [varargout{1:nargout}] = mpsqueeze(varargin{:});               
        end

        function varargout = shiftdim(varargin)
        %SHIFTDIM Shift dimensions.
        %
        %    Usage is identical to built-in function SHIFTDIM.
        %
        %    See also SQUEEZE.
        
             [varargout{1:nargout}] = mpshiftdim(varargin{:});               
        end

        function varargout = circshift(varargin)
        %CIRCSHIFT Shift array circularly.
        %
        %    Usage is identical to built-in function CIRCSHIFT.
        %
        %    See also FFTSHIFT, SHIFTDIM, PERMUTE.
        
             [varargout{1:nargout}] = mpcircshift(varargin{:});               
        end
        
        function varargout = blkdiag(varargin)
        %BLKDIAG Block diagonal concatenation of matrix input arguments.
        %
        %    Usage is identical to built-in function BLKDIAG.
        %
        %    See also BLKDIAG.
        
             [varargout{1:nargout}] = mpblkdiag(varargin{:});               
        end
        
        %% Basic Information
        function display(x) %#ok<DISPLAY>
        %DISPLAY Display multiprecision entity
        %
        %    Usage is identical to built-in function DISPLAY.
        %
        %    See also DISPLAY.
        
            disp(mpimpl(2,x,inputname(1)));
        end
        
        function disp(x)
        %DISP Display multiprecision entity        
        %
        %    Usage is identical to built-in function DISP.
        %
        %    See also DISP.
        
            disp(mpimpl(2,x,inputname(1)));
        end
        
        function r = isnan(x)
        %ISNAN Array elements that are NaN
        %
        %    Usage is identical to built-in function ISNAN.
        %
        %    See also ISNAN.
        
           r = mpimpl(800,x); 
        end
        
        function r = isinf(x)
        %ISINF Array elements that are INF        
        %
        %    Usage is identical to built-in function ISINF.
        %
        %    See also ISINF.
        
           r = mpimpl(801,x); 
        end
        
        function r = isfinite(x)
        %ISFINITE Array elements that are finite
        %
        %    Usage is identical to built-in function ISFINITE.
        %
        %    See also ISFINITE.
        
           r = mpimpl(802,x); 
        end
        
        function varargout = hasInfNaN(varargin)
        %HASINFNAN Checks if array has Inf or NaN entries.
        %
        %    HASINFNAN(X) returns true if X has any Inf or NaN entry.
        %
        %    See also ISINF, ISNAN, ISFINITE.
        
            [varargout{1:nargout}] = mpimpl(1037,varargin{:});
        end
        
        function r = isnumeric(varargin)
        %ISNUMERIC Determine if input is numeric array
        %
        %    Usage is identical to built-in function ISNUMERIC.
        %
        %    See also ISNUMERIC.
        
           r = true; 
        end

        function r = isinteger(varargin)
        %ISINTEGER True for arrays of integer data type.
        %
        %    Usage is identical to built-in function ISINTEGER.
        %
        %    See also ISA, ISNUMERIC, ISFLOAT.
        
           r = false; 
        end
        
        function r = isfloat(varargin)
        %ISFLOAT Determine if input is floating-point array
        %
        %    Usage is identical to built-in function ISFLOAT.
        %
        %    See also ISFLOAT.
        
           r = true; 
        end
        
        function r = isempty(x)
        %ISEMPTY Determine whether array is empty
        %
        %    Usage is identical to built-in function ISEMPTY.
        %
        %    See also ISEMPTY.
        
           r = numel(x) == 0; 
        end
        
        function r = isa(x,category)
        %ISA Returns true if category is 'mp', 'float' or 'numeric' and false otherwise
        %
        %   See also ISNUMERIC, ISLOGICAL, ISCHAR, ISCELL, ISSTRUCT, ISFLOAT,
        %            ISINTEGER, ISOBJECT, ISJAVA, ISSPARSE, ISREAL, CLASS.
        
            if ischar(category)
                r = strcmpi(class(x),'mp') && any(strcmpi(category,{'mp','float','numeric'}));
            else
                error('Must be a text scalar.');
            end
        end

        %% Arithemtic operators        
        function r = minus(x,y)
        %MINUS Minus
        %
        %    Usage is identical to built-in function MINUS.
        %
        %    See also MINUS.
        
           r = mpimpl(10,x,y); 
        end
        
        function r = plus(x,y)
        %PLUS Plus
        %
        %    Usage is identical to built-in function PLUS.
        %
        %    See also PLUS.
        
           r = mpimpl(11,x,y); 
        end
        
        function r = uplus(x)
        %UPLUS Unary plus
        %
        %    Usage is identical to built-in function UPLUS.
        %
        %    See also UPLUS.
        
           r = x;
        end
        
        function r = uminus(x)
        %UMINUS Unary minus
        %
        %    Usage is identical to built-in function UMINUS.
        %
        %    See also UMINUS.
        
           r = mpimpl(160,x); 
        end

        function r = times(x,y)
        %TIMES Array multiply
        %
        %    Usage is identical to built-in function TIMES.
        %
        %    See also TIMES.
        
           r = mpimpl(15,x,y);
        end
       
        function r = mtimes(x,y)
        %MTIMES Matrix multiplication
        %
        %    Usage is identical to built-in function MTIMES.
        %
        %    See also MTIMES.
        
           r = mpimpl(20,x,y); 
        end
        
        function r = rdivide(x,y)
        %RDIVIDE Right array division
        %
        %    Usage is identical to built-in function RDIVIDE.
        %
        %    See also RDIVIDE.
        
           r = mpimpl(25,x,y); 
        end
        
        function r = ldivide(x,y)
        %LDIVIDE Left array division
        %
        %    Usage is identical to built-in function LDIVIDE.
        %
        %    See also LDIVIDE.
        
           r = mpimpl(30,x,y); 
        end
        
        function r = mrdivide(B,A)
        %MRDIVIDE Solve systems of linear equations xA = B for x
        %
        %    Usage is identical to built-in function MRDIVIDE.
        %
        %    See also MRDIVIDE.
           if isscalar(A)
               r = mpimpl(25,B,A); 
           else
               r = (A'\B')';
           end;
        end
        
        function x = mldivide(A,B)
        %MLDIVIDE Solve systems of linear equations Ax = B for x
        %
        %    Details of implementation can be found in the article:
        %    http://www.advanpix.com/2016/10/07/architecture-of-linear-systems-solver/
        %
        %    Usage is identical to built-in function MLDIVIDE.
        %
        %    See also MLDIVIDE.
        
           w = warning('query','MATLAB:nearlySingularMatrix');
           if strcmp(w.state,'off')
              % Skip rcond estimation if warning is disabled.
              % This gives us better performance, since rcond estimation 
              % takes 20%-50% of overall time in solver.
              x = mpimpl(40,A,B,'norcond');
           else
              x = mpimpl(40,A,B); 
           end;
        end
        
        function r = mpower(x,y)
        %^ Matrix power
        %
        %    Usage is identical to built-in function MPOWER.
        %
        %    See also MPOWER.
        
            if ~isscalar(x) && ~isscalar(y)
                error('MCT:mpower:IncorrectInputs',...
                           'Inputs must be a scalar and a square matrix.' );
                
            end
        
            if isscalar(x) && ~isscalar(y)
                 r = expm(log(x)*y);
            end
            
            if ~isscalar(x) && isscalar(y)
                if isreal(y) && y >= 0 && fix(y) == y
                    r = mpimpl(45,x,y); % fast binary squaring in C++
                else
                    r = expm(logm(x)*y);
                end
            end
            
            if isscalar(x) && isscalar(y)            
                r = mpimpl(310,x,y); % call element-wise power (.^)
            end
        end
        
        function r = ctranspose(x)
        %CTRANSPOSE Complex conjugate transpose
        %
        %    Usage is identical to built-in function CTRANSPOSE.
        %
        %    See also CTRANSPOSE.
        
           r = mpimpl(50,x); 
        end
        
        function r = transpose(x)
        %TRANSPOSE Transpose
        %
        %    Usage is identical to built-in function TRANSPOSE.
        %
        %    See also TRANSPOSE.
        
           r = mpimpl(55,x); 
        end

        function varargout = colon(varargin)
        %COLON Create vectors, array subscripting, and for-loop iterators
        %
        %    Usage is identical to built-in function COLON.
        %
        %    See also COLON.
        
             [varargout{1:nargout}] = mpimpl(400, varargin{:});               
        end
        
        %% Type Conversion
        function r = double(x)
        %DOUBLE Convert to double precision.
        %
        %    Usage is identical to built-in function DOUBLE.
        %
        %    See also DOUBLE.
        
           r = mpimpl(402,x); 
        end
        
        function r = single(x)
        %SINGLE Convert to single precision.
        %
        %    Usage is identical to built-in function SINGLE.
        %
        %    See also SINGLE.
        
           r =  mpimpl(9000,x); 
        end
        
        function r = int8(x)
        %INT8 Convert to signed 8-bit integer.
        %
        %    Usage is identical to built-in function INT8.
        %
        %    See also INT8.
        
           r =  mpimpl(9010,x);  
        end
        
        function r = uint8(x)
        %UINT8 Convert to unsigned 8-bit integer.
        %
        %    Usage is identical to built-in function UINT8.
        %
        %    See also UINT8.
        
           r = mpimpl(9015,x); 
        end
        
        function r = int16(x)
        %INT16 Convert to signed 16-bit integer.
        %
        %    Usage is identical to built-in function INT16.
        %
        %    See also INT16.
        
           r = mpimpl(9020,x); 
        end
        
        function r = uint16(x)
        %UINT16 Convert to unsigned 16-bit integer.
        %
        %    Usage is identical to built-in function UINT16.
        %
        %    See also UINT16.
        
           r = mpimpl(9025,x);  
        end
        
        function r = int32(x)
        %INT32 Convert to signed 32-bit integer.
        %
        %    Usage is identical to built-in function INT32.
        %
        %    See also INT32.
        
           r = mpimpl(9030,x);  
        end
        
        function r = uint32(x)
        %UINT32 Convert to unsigned 32-bit integer.
        %
        %    Usage is identical to built-in function UINT32.
        %
        %    See also UINT32.
        
           r = mpimpl(9035,x); 
        end
        
        function r = int64(x)
        %INT64 Convert to signed 64-bit integer.
        %
        %    Usage is identical to built-in function INT64.
        %
        %    See also INT64.
        
           r = mpimpl(9040,x);  
        end
        
        function r = uint64(x)
        %UINT64 Convert to unsigned 64-bit integer.
        %
        %    Usage is identical to built-in function UINT64.
        %
        %    See also UINT64.
        
           r = mpimpl(9045,x); 
        end
        
        function bits = getbits(x)
        %GETBITS Retrieves raw data as array of unsigned 64-bit integers.
        %
        %    The function provides low-level access to underlying data of 
        %    multiple-precision floating-point numbers - designed to be used by developers.
        %
        %    Multiple-precision floating-point number are stored in memory 
        %    as array of unsigned 64-bit integers. Function returns uint64 matrix 
        %    where i-th row corresponds to x(i) as
        % 
        %    data(i,:) = [precision, sign, exponent, dummy1, dummy2, mantissa(1)...mantissa(n)]
        %
        %    Please use GETWORDS64 to get access to bits of quadruple precision numbers.
        % 
        %    Extracted integers can be stored, transferred to
        %    another system and then multiple-precision floats can be
        %    reconstructed by SETBITS (not implemented yet).
        %
        %    See also GETWORDS64, SETWORDS64, SETBITS.
            
           bits = mpimpl(9060,x); 
        end
        
        %% Trigonometric Functions
        function r = sin(x) 
        %SIN Sine of argument in radians.
        %
        %    Usage is identical to built-in function SIN.
        %
        %    See also SIN.
        
            r = mpimpl(200,x); 
        end
        
        function r = cos(x)
        %COS Cosine of argument in radians.
        %
        %    Usage is identical to built-in function COS.
        %
        %    See also COS.
        
           r = mpimpl(201,x); 
        end
        
        function r = tan(x)
        %TAN Tangent of argument in radians.
        %
        %    Usage is identical to built-in function TAN.
        %
        %    See also TAN.
        
           r = mpimpl(202,x); 
        end
        
        function r = sec(x)
        %SEC Secant of argument in radians.
        %
        %    Usage is identical to built-in function SEC.
        %
        %    See also SEC.
        
           r = mpimpl(203,x); 
        end
        
        function r = csc(x)
        %CSC Cosecant of argument in radians.
        %
        %    Usage is identical to built-in function CSC.
        %
        %    See also CSC.
        
           r = mpimpl(204,x); 
        end
        
        function r = cot(x)
        %COT Cotangent of argument in radians.
        %
        %    Usage is identical to built-in function COT.
        %
        %    See also COT.
        
           r = mpimpl(205,x); 
        end
        
        function r = acos(x)
        %ACOS Inverse cosine, result in radians.
        %
        %    Usage is identical to built-in function ACOS.
        %
        %    See also ACOS.
        
           r = mpimpl(206,x); 
        end
        
        function r = asin(x)
        %ASIN Inverse sine, result in radians.
        %
        %    Usage is identical to built-in function ASIN.
        %
        %    See also ASIN.
        
           r = mpimpl(207,x); 
        end
        
        function r = atan(x)
        %ATAN Inverse tangent, result in radians.
        %
        %    Usage is identical to built-in function ATAN.
        %
        %    See also ATAN.
        
           r = mpimpl(208,x); 
        end
        
        function r = acot(x)
        %ACOT Inverse cotangent, result in radian.
        %
        %    Usage is identical to built-in function ACOT.
        %
        %    See also ACOT.
        
           r = mpimpl(220,x); 
        end
        
        function r = asec(x)
        %ASEC  Inverse secant, result in radians.
        %
        %    Usage is identical to built-in function ASEC.
        %
        %    See also ASEC.
        
           r = mpimpl(221,x); 
        end
        
        function r = acsc(x)
        %ACSC Inverse cosecant, result in radian.
        %
        %    Usage is identical to built-in function ACSC.
        %
        %    See also ACSC.
        
           r = mpimpl(222,x); 
        end
        
        function r = cosh(x)
        %COSH Hyperbolic cosine.
        %
        %    Usage is identical to built-in function COSH.
        %
        %    See also COSH.
        
           r = mpimpl(209,x); 
        end

        function r = sinh(x)
        %SINH Hyperbolic sine.
        %
        %    Usage is identical to built-in function SINH.
        %
        %    See also SINH.
        
           r = mpimpl(210,x); 
        end
       
        function r = tanh(x)
        %TANH Hyperbolic tangent.
        %
        %    Usage is identical to built-in function TANH.
        %
        %    See also TANH.
        
           r = mpimpl(211,x); 
        end
        
        function r = sech(x)
        %SECH Hyperbolic secant.
        %
        %    Usage is identical to built-in function SECH.
        %
        %    See also SECH.
        
           r = mpimpl(212,x); 
        end
        
        function r = csch(x)
        %CSCH Hyperbolic cosecant.
        %
        %    Usage is identical to built-in function CSCH.
        %
        %    See also CSCH.
        
           r = mpimpl(213,x); 
        end
        
        function r = coth(x)
        %COTH Hyperbolic cotangent.
        %
        %    Usage is identical to built-in function COTH.
        %
        %    See also COTH.
        
           r = mpimpl(214,x); 
        end
        
        function r = acosh(x)
        %ACOSH Inverse hyperbolic cosine.
        %
        %    Usage is identical to built-in function ACOSH.
        %
        %    See also ACOSH.
        
           r = mpimpl(215,x); 
        end
        
        function r = asinh(x)
        %ASINH Inverse hyperbolic sine.
        %
        %    Usage is identical to built-in function ASINH.
        %
        %    See also ASINH.
        
           r = mpimpl(216,x); 
        end
        
        function r = atanh(x)
        %ATANH Inverse hyperbolic tangent.
        %
        %    Usage is identical to built-in function ATANH.
        %
        %    See also ATANH.
        
           r = mpimpl(217,x); 
        end
        
        function r = acoth(x)
        %ACOTH Inverse hyperbolic cotangent.
        %
        %    Usage is identical to built-in function ACOTH.
        %
        %    See also ACOTH.
        
           r = mpimpl(223,x); 
        end
        
        function r = asech(x)
        %ASECH Inverse hyperbolic secant.
        %
        %    Usage is identical to built-in function ASECH.
        %
        %    See also ASECH.
        
           r = mpimpl(224,x); 
        end
        
        function r = acsch(x)
        %ACSCH Inverse hyperbolic cosecant.
        %
        %    Usage is identical to built-in function ACSCH.
        %
        %    See also ACSCH.

           r = mpimpl(225,x); 
        end
        
        function r = atan2(y,x)
        %ATAN2 Four quadrant inverse tangent.
        %
        %    Usage is identical to built-in function ATAN2.
        %
        %    See also ATAN2
        
           r = mpimpl(218,y,x); 
        end
        
        function r = hypot(x,y)
        %HYPOT Robust computation of the square root of the sum of squares.
        %
        %    Usage is identical to built-in function HYPOT.
        %
        %    See also HYPOT
            
           r = mpimpl(219,x,y); 
        end
        
        %% Matrix Functions
        function varargout = funm(varargin)
        %FUNM Evaluate general matrix function.
        %
        %    Usage is identical to built-in function FUNM.
        %
        %    See also FUNM
        
            [varargout{1:nargout}] = mpfunm(varargin{:});
        end
        
        function varargout = sqrtm(varargin)
        %SQRTM Matrix square root.
        %
        %    Usage is identical to built-in function SQRTM.
        %
        %    See also SQRTM
        
            [varargout{1:nargout}] = mpsqrtm(varargin{:});
        end
        
        function varargout = sqrtm_tri(varargin)
        %SQRTM_TRI Square root of upper triangular matrix.
        %
        %    Heavily optimized routine to compute square root of upper triangular matrix.
        %    Only upper triangular part of input matrix is accessed.
        %
        %    See also SQRTM
        
            [varargout{1:nargout}] = mpimpl(4007,varargin{:});
        end
        
        function varargout = expm(varargin)
        %EXPM Matrix exponential.
        %
        %    Usage is identical to built-in function EXPM.
        %
        %    See also EXPM

            [varargout{1:nargout}] = mpexpm(varargin{:});                    
        end
       
        function r = logm(x)
        %LOGM Matrix logarithm.
        %
        %    Usage is identical to built-in function LOGM.
        %
        %    See also LOGM
        
           r = mpfunm(x,@log);
        end
        
        function r = sinm(x)
        %SINM Matrix sine.
        %
        %    Usage is identical to built-in function SINM.
        %
        %    See also SINM
        
           r = mpfunm(x,@sin);  %r = mpimpl(4003,x);
        end
        
        function r = cosm(x)
        %COSM Matrix cosine.
        %
        %    Usage is identical to built-in function COSM.
        %
        %    See also COSM
        
           r = mpfunm(x,@cos);  %r = mpimpl(4004,x);
        end
        
        function r = sinhm(x)
        %SINHM Matrix hyperbolic sine.
        %
        %    Usage is identical to built-in function SINHM.
        %
        %    See also SINHM
        
           r = mpfunm(x,@sinh); %r = mpimpl(4005,x);        
        end
        
        function r = coshm(x)
        %COSHM Matrix hyperbolic cosine.
        %
        %    Usage is identical to built-in function COSHM.
        %
        %    See also COSHM
        
           r = mpfunm(x,@cosh); %r = mpimpl(4006,x);
        end
        
        function varargout = arrayfun(func, varargin)
        %ARRAYFUN Apply a function to each element of an array.
        %
        %   Please be aware, this is early alpha implementation, use it with
        %   caution. Report issues/suggest improvements to support@advanpix.com
        %   
        %   See also  CELLFUN, STRUCTFUN, FUNCTION_HANDLE, CELL2MAT, SPFUN
           
            [varargout{1:nargout}] = mparrayfun(func, varargin{:});         
        end        
        
        function varargout = accumarray(varargin)
        %ACCUMARRAY Construct an array by accumulation.
        %
        %   Usage is identical to built-in function ACCUMARRAY.
        %
        %   Please be aware, this is early alpha implementation and can be
        %   slow for large arrays. Report issues/suggest improvements to support@advanpix.com
        %
        %   See also FULL, SPARSE, SUM, FUNCTION_HANDLE.
        
           [varargout{1:nargout}] = mpaccumarray(varargin{:}); 
        end
        
        %% Exponential Functions
        function r = exp(x)
        %EXP Exponential.
        %
        %    Usage is identical to built-in function EXP.
        %
        %    See also EXP
        
           r = mpimpl(300,x); 
        end
        
        function r = expm1(x)
        %EXPM1 Compute EXP(X)-1 accurately.
        %
        %    Usage is identical to built-in function EXPM1.
        %
        %    See also EXPM1
        
           r = mpimpl(301,x); 
        end
        
        function r = log(x)
        %LOG Natural logarithm.
        %
        %    Usage is identical to built-in function LOG.
        %
        %    See also LOG
        
           r = mpimpl(302,x); 
        end
        
        function r = log10(x)
        %LOG10 Common (base 10) logarithm.
        %
        %    Usage is identical to built-in function LOG10.
        %
        %    See also LOG10
        
           r = mpimpl(303,x); 
        end
        
        function r = log1p(x)
        %LOG1P Compute LOG(1+X) accurately.
        %
        %    Usage is identical to built-in function LOG1P.
        %
        %    See also LOG1P
        
           r = mpimpl(304,x); 
        end
        
        function r = log2(x)
        %LOG2 Base 2 logarithm and dissect floating point number.
        %
        %    Usage is identical to built-in function LOG2.
        %
        %    See also LOG2
        
           r = mpimpl(305,x); 
        end
        
        function r = nextpow2(x)
        %NEXTPOW2 Next higher power of 2.
        %
        %    Usage is identical to built-in function NEXTPOW2.
        %
        %    See also NEXTPOW2
        
           r = mpimpl(306,x); 
        end
        
        function r = nthroot(x,N)
        %NTHROOT Real n-th root of real numbers.
        %
        %    Usage is identical to built-in function NTHROOT.
        %
        %    See also NTHROOT
        
           r = mpimpl(307,x,N); 
        end
        
        function r = pow2(x,e)
        %POW2 Base 2 power and scale floating point number.
        %
        %    Usage is identical to built-in function POW2.
        %
        %    See also POW2
        
           if nargin == 1,     r = mpimpl(308,x); 
           elseif nargin == 2, r = mpimpl(309,x,e); end 
        end
        
        function r = power(x,y)
        %.^  Array power.
        %
        %    Usage is identical to built-in function POWER.
        %
        %    See also POWER

           r = mpimpl(310,x,y);
        end
        
        function r = reallog(x)
        %REALLOG Real logarithm.
        %
        %    Usage is identical to built-in function REALLOG.
        %
        %    See also REALLOG
        
           r = mpimpl(311,x); 
        end
        
        function r = realpow(x,y)
        %REALPOW Real power.
        %
        %    Usage is identical to built-in function REALPOW.
        %
        %    See also REALPOW
        
           r = mpimpl(312,x,y);
        end
        
        function r = realsqrt(x)
        %REALSQRT Real square root.
        %
        %    Usage is identical to built-in function REALSQRT.
        %
        %    See also REALSQRT
        
           r = mpimpl(313,x); 
        end
        
        function r = sqrt(x)
        %SQRT Square root.
        %
        %    Usage is identical to built-in function SQRT.
        %
        %    See also SQRT
        
           r = mpimpl(314,x); 
        end
        
        %% Signal Processing Functions
        function varargout = sinc(varargin)
        %SINC Sin(pi*x)/(pi*x) function.
        %
        %    Usage is identical to built-in function SINC.
        %
        %    See also SIN, COS.
        
            [varargout{1:nargout}] = mpsinc(varargin{:});        
        end
        
        %% Error and related functions
        function r = erf(x)
        %ERF Error function.
        %
        %    Usage is identical to built-in function ERF.
        %
        %   See also ERFC, ERFCX, ERFINV, ERFCINV.
            
           r = mpimpl(403,x); 
        end
        
        function r = erfc(x)
        %ERFC Complementary error function.
        %
        %    Usage is identical to built-in function ERFC.
        %
        %   See also ERF, ERFCX, ERFINV, ERFCINV.
            
           r = mpimpl(404,x); 
        end
        
        function r = erfinv(x)
        %ERFINV Inverse error function.
        %
        %    Usage is identical to built-in function ERFINV.
        %
        %   See also ERF, ERFC, ERFCX, ERFCINV.            
        
           r = mpimpl(440,x); 
        end
        
        function r = erfcinv(x)
        %ERFCINV Inverse complementary error function.
        %
        %    Usage is identical to built-in function ERFCINV.
        %
        %   See also ERF, ERFC, ERFCX, ERFINV.            
        
           r = mpimpl(441,x); 
        end
        
        function varargout = norminv(varargin)
        %NORMINV Inverse of the normal cumulative distribution function (cdf).
        %
        %    Usage is identical to built-in function NORMINV.
        %
        %   See also ERFINV, ERFCINV
        
           [varargout{1:nargout}] = mpnorminv(varargin{:}); 
        end

        function r = erfi(x)
        %ERFI  Imaginary error function.
        %
        %    Usage is identical to built-in function ERFI.
        %
        %    See also ERF, ERFC, ERFINV, ERFCX, ERFCINV
            
           r = mpimpl(419,x); 
        end

        function r = fresnels(x)
        %FRESNELS Computes the Fresnel integral S(z).  
        %
        %    Usage is identical to built-in function FRESNELS.
        %
        %    See also FRESNELC.
        
           r = mpimpl(420,x); 
        end

        function r = fresnelc(x)
        %FRESNELC Computes the Fresnel integral C(z).
        %
        %    Usage is identical to built-in function FRESNELC.
        %
        %    See also FRESNELS.
        
           r = mpimpl(421,x); 
        end
        
        %% Gamma and related functions
        function r = gamma(z)
        %GAMMA Gamma function.
        %
        %    GAMMA(Z) evaluates gamma function for each element Z.
        %    Argument Z can be real or complex. 
        %
        %    Built-in function GAMMA are limited only to real arguments Z.
        %    Otherwise usage is identical.
        %
        %   See also GAMMALN, GAMMAINC, GAMMAINCINV, PSI.
            
           r = mpimpl(405,z); 
        end
        
        function varargout = gammainc(varargin)
        %GAMMAINC Incomplete gamma function.
        %
        %   The function accepts arbitrary complex or real arguments.
        %   Options 'scaledlower' / 'scaledupper' are not supported.
        %   Otherwise usage is identical to built-in function GAMMAINC.
        %
        %   See also GAMMAINCINV, GAMMA, GAMMALN, PSI.
            
            [varargout{1:nargout}] = mpimpl(418,varargin{:});
        end
        
        function r = gammaln(z)
        %GAMMALN Logarithm of gamma function.
        %
        %    GAMMALN(Z) evaluates logarithm of gamma function for each element Z.
        %    Argument Z can be real or complex. 
        %
        %    Built-in function GAMMALN are limited only to real arguments Z.
        %    Otherwise usage is identical.
        %
        %    See also GAMMA, GAMMAINC, GAMMAINCINV, PSI.
        
           r = mpimpl(406,z); 
        end
        
        function varargout = legendre(varargin)
        %LEGENDRE Associated Legendre function.
        %
        %    Usage is identical to built-in function LEGENDRE.
        
           [varargout{1:nargout}] = mplegendre(varargin{:}); 
        end
        
        function r = psi(varargin)
        %PSI  Psi (polygamma) function.
        %
        %   Y = PSI(Z) evaluates the digamma function for each element of Z.
        %   Argument Z can be real or complex. The digamma function, 
        %   is the logarithmic derivative of the gamma function: 
        %
        %      psi(z) = digamma(z) = d(log(gamma(z)))/dz = (d(gamma(z))/dz)/gamma(z).
        %
        %   Y = PSI(S,Z) evaluates generalized polygamma function. If S is a nonnegative 
        %   integer, this is simply the S-order derivative of the digamma function.
        %   Generalization to any order values S is due to: 
        %        O. Espinosa and V. Moll, "A generalized polygamma function", 
        %        Integral Transforms and Special Functions (2004), 101-115.
        %
        %   See also GAMMA, GAMMALN, GAMMAINC, GAMMAINCINV.
            
           r = mpimpl(407,varargin{:}); 
        end

        %% Airy Functions
        function varargout = airy(varargin)
        %AIRY  Airy functions.
        %
        %    Real and complex arguments are supported.
        %    Usage is identical to built-in function AIRY.
        %
        %    See also BESSELH, BESSELI, BESSELJ, BESSELK, BESSELY.

            [varargout{1:nargout}] = mpimpl(448,varargin{:});                    
        end
        
        %% Bessel Functions
        function varargout = besselj(varargin)
        %BESSELJ Bessel function of the first kind.
        %
        %    Usage is identical to built-in function BESSELJ.
        %
        %    See also BESSELJ.
        
            [varargout{1:nargout}] = mpimpl(424, varargin{:});
        end

        function varargout = bessely(varargin)
        %BESSELY Bessel function of the second kind.
        %
        %    Usage is identical to built-in function BESSELY.
        %
        %    See also BESSELY
        
            [varargout{1:nargout}] = mpimpl(425, varargin{:});
        end

        function varargout = besseli(varargin)
        %BESSELI Modified Bessel function of the first kind.
        %
        %    Usage is identical to built-in function BESSELI.
        %
        %    See also BESSELI.
        
            [varargout{1:nargout}] = mpimpl(426, varargin{:});
        end

        function varargout = besselk(varargin)
        %BESSELK Modified Bessel function of the second kind.
        %
        %    Usage is identical to built-in function BESSELK.
        %
        %    See also BESSELK.
        
            [varargout{1:nargout}] = mpimpl(427, varargin{:});
        end

        function varargout = besselh(varargin)
        %BESSELH Inverse cosine, result in radians.
        %
        %    Usage is identical to built-in function BESSELH.
        %
        %    See also BESSELH.
        
            [varargout{1:nargout}] = mpimpl(428, varargin{:});
        end

        function y = beta(z,w)
        %BETA Beta function.
        %
        %   Usage is identical to built-in function BETA.
        %
        %   See also BETAINC, BETALN.
        
            y = exp(betaln(z,w));
        end
        
        function y = betaln(z,w)
        %BETALN Logarithm of beta function.
        %
        %   Usage is identical to built-in function BETALN.
        %
        %   See also BETAINC, BETA.
        
            y = gammaln(z)+gammaln(w)-gammaln(z+w);
        end

        function varargout = betainc(varargin)
        %BETAINC Incomplete beta function.
        %
        %   The function accepts arbitrary complex or real arguments.
        %   Otherwise usage is identical to built-in function BETAINC.
        %
        %   See also BETALN, BETA.
        
            [varargout{1:nargout}] = mpimpl(445, varargin{:});
        end
        
        function varargout = ellipj(varargin)
        %ELLIPJ Jacobi elliptic functions.
        %
        %    Usage is identical to built-in function ELLIPJ.
        %
        %    See also ELLIPKE.
        
           [varargout{1:nargout}] = mpellipj(varargin{:}); 
        end

        function varargout = ellipke(varargin)
        %ELLIPKE Complete elliptic integral.
        %
        %    Real and complex arguments are supported.
        %    Otherwise usage is identical to built-in function ELLIPKE.
        %
        %    See also ELLIPJ.
        
            [varargout{1:nargout}] = mpimpl(449,varargin{:});
        end
        
        %% Zeta Functions
        function r = zeta(varargin)
        %ZETA  Riemann zeta function.
        %
        %   Z = zeta(Z) returns the Riemann zeta function.
        %
        %   See also HURWITZZETA.
        
            if nargin==1
                r = mpimpl(408,varargin{1}); 
            else
               error('%s\n%s','Two-parameter call to ZETA is not supported in multiprecision mode.',...
                              'Please request it using support@advanpix.com');
            end
        end

        function r = hurwitzZeta(varargin)
        %HURWITZZETA Hurwitz zeta function.
        %
        %   Y = hurwitzZeta(S,Z) computes the Hurwitz zeta function.
        %
        %   See also ZETA.
            
            r = mpimpl(446,varargin{:});             
        end
        
        function varargout = LerchPhi(varargin)
        % LERCHPHI Lerch transcendent function
        
           error('%s\n%s','Lerch transcendent is not provided in multiprecision toolbox.',...
                          'Please request it using support@advanpix.com');
        end

        %% Integral Functions En(z), n = 1 by default.
        function r = expint(varargin)
        % EXPINT Generalized exponential integral.
        %
        %   Y = EXPINT(S,Z) computes generalized exponential integral for each
        %   element of S and Z. Generalized exponential integral is defined as:
        %
        %             E(s,z) = int(exp(-z*t)/t^s,t=1..Inf)
        %
        %   Y = EXPINT(X) is equivalent to built-in function EXPINT.
        %
        %   See also EINT, COSINT, SININT, LOGINT, COSHINT, SINHINT.

            if nargin==1,  r = mpimpl(409, 1, varargin{1});
            else           r = mpimpl(409, varargin{1}, varargin{2});  end %#ok<SEPEX>
        end

        function r = eint(z)
        % EINT Exponential integral Ei(z).
        %
        %   Y = EINT(Z) computes exponential integral for each
        %   element of Z. Exponential integral is defined as:
        %
        %             Ei(z) = -int(exp(-t)/t,t=-x..Inf)
        %
        %   See also EXPINT, COSINT, SININT, LOGINT, COSHINT, SINHINT.
            
            r = mpimpl(410, z);
        end
       
        function r = cosint(z)
        %COSINT Cosine integral function.
        %
        %  COSINT(x) = eulergamma + log(x) + int((cos(t)-1)/t,t,0,x),
        %  where eulergamma denotes the Euler-Mascheroni constant. 
        %
        %  See also EINT, EXPINT, SININT, LOGINT, COSHINT, SINHINT.
            
            r = mpimpl(412, z);
        end

        function r = sinint(z)
        %SININT Sine integral function.
        %
        %  SININT(x) = int(sin(t)/t,t,0,x).
        %
        %  See also EINT, EXPINT, COSINT, LOGINT, COSHINT, SINHINT.
            
            r = mpimpl(413, z);
        end
        
        function r = logint(z)
        % LOGINT Logarithmic integral
        %
        %  Y = LOGINT(Z) is logarithmic integral of Z:
        %
        %        LOGINT(x) = li(z) = Ei(log(z))
        %
        %  See also EINT, EXPINT, COSINT, SININT, COSHINT, SINHINT.
           
            r = mpimpl(411, z);
        end
        
        function r = coshint(z)
        % COSHINT  Hyperbolic cosine integral function.
        %
        %  Y = COSHINT(Z) computes the hyperbolic cosine integral Chi(z). 
        %
        %  See also EINT, EXPINT, COSINT, SININT, LOGINT, SINHINT.
            
            r = mpimpl(414, z);
        end

        function r = sinhint(z)
        % SINHINT  Hyperbolic sine integral function.
        %
        %  Y = SINHINT(Z) computes the hyperbolic sine integral Shi(z). 
        %
        %  See also EINT, EXPINT, COSINT, SININT, LOGINT, COSHINT.
            
            r = mpimpl(415, z);
        end

        %% Hypergeometric functions
        function r = hypergeom(a, b, z)
        %HYPERGEOM Generalized hypergeometric function.
        %
        %    Computes generalized hypergeometric function:
        %
        %         pFq([a_1,...a_p],[b_1,...b_q],z)
        %
        %    Any combination of real and complex arguments are supported.
        %    Vectors a and b can have arbitrary length (p and q).
        %
        %    See also KUMMERU, KUMMERM.
        
            r = mpimpl(416, a, b, z);
        end

        function r = kummerM(a, b, z, scale)
        %KUMMERM Kummer's (confluent hypergeometric) function M(a,b,z).
        %
        %    Computes Kummer's (confluent hypergeometric) function
        %
        %         M(a,b,z) = 1F1(a,b,z),          scale = 0 (default)
        %         M(a,b,z) = 1F1(a,b,z)/gamma(b), scale = 1 
        %
        %    KummerM and KummerU functions solve the differential equation
        %
        %                 zw''+(b-z)w'-aw=0 
        %
        %    See also KUMMERU, HYPERGEOM.
        
            if nargin==3 
                r = mpimpl(422, a, b, z);
            elseif nargin==4  
                r = mpimpl(422, a, b, z, scale);
            end
        end

        function r = kummerU(a, b, z)
        %KUMMERU Tricomi's (confluent hypergeometric) function U(a,b,z).
        %
        %    Computes Tricomi's (confluent hypergeometric) function
        %
        %    U(a,b,z) = gamma(1-b)/gamma(a+1-b)*M(a,b,z)+gamma(b-1)/gamma(a)*z^(1-b)*M(a+1-b,2-b,z)
        %
        %    KummerM and KummerU functions solve the differential equation
        %
        %                   zw''+(b-z)w'-aw=0 
        %
        %    Usage is identical to built-in function KUMMERU.
        %
        %    See also KUMMERM, HYPERGEOM.
        
            r = mpimpl(423, a, b, z);
        end
        
        %% Lambert W function
        function varargout = lambertw(varargin)
        %LAMBERTW Lambert's W function.
        %
        %   Computes the Lambert W function, which solves the equation W*exp(W) = Z.
        %
        %   W = LAMBERTW(Z) solves W*exp(W) = Z.
        %   W = LAMBERTW(K,Z) is the K-th branch of this multi-valued function.
        %
        %   Real and complex arguments are supported.
        %
        
           [varargout{1:nargout}] = mpimpl(447,varargin{:});         
        end
        
        %% Orthogonal Polynomials
        function varargout = legendreP(varargin)
        %LEGENDREP Legendre polynomials.
        %
        %    LEGENDREP(N,X) computes values of polynomial PN(x), where:
        %      N - polynomial degree, must be nonnegative integer matrix or n-dim array.
        %      X - matrix or n-dim array of real values.
        %
        %    LEGENDREP(N,X,'all') computes values of all polynomials up to Nth degree in X points:
        %      N - maximum degree must be nonnegative scalar.
        %      X - real vector of length M.
        %
        %     P0(X(1)), P1(X(1)), P2(X(1)),...,PN(X(1))
        %     P0(X(2)), P1(X(2)), P2(X(2)),...,PN(X(2))
        %     ...
        %     P0(X(M)), P1(X(M)), P2(X(M)),...,PN(X(M))
        %
        %    See also CHEBYSHEVT, CHEBYSHEVU, CHEBYSHEVV, CHEBYSHEVW, HERMITH.
         
           [varargout{1:nargout}] = mpimpl(430,varargin{:}); 
        end

        function varargout = chebyshevT(varargin)
        %CHEBYSHEVT Chebyshev polynomials of the first kind.
        %
        %    Usage syntax is identical to legendreP function.
        %    Please see its help for details.
        %
        %    See also LEGENDREP, CHEBYSHEVU, CHEBYSHEVV, CHEBYSHEVW, HERMITEH.
         
           [varargout{1:nargout}] = mpimpl(431,varargin{:}); 
        end

        function varargout = chebyshevU(varargin)
        %CHEBYSHEVU Chebyshev polynomials of the second kind.
        %
        %    Usage syntax is identical to legendreP function.
        %    Please see its help for details.
        %
        %    See also LEGENDREP, CHEBYSHEVT, CHEBYSHEVV, CHEBYSHEVW, HERMITEH.
         
           [varargout{1:nargout}] = mpimpl(432,varargin{:}); 
        end

        function varargout = chebyshevV(varargin)
        %CHEBYSHEVV Chebyshev polynomials of the third kind.
        %
        %    Usage syntax is identical to legendreP function.
        %    Please see its help for details.
        %
        %    See also LEGENDREP, CHEBYSHEVT, CHEBYSHEVU, CHEBYSHEVW, HERMITEH.
         
           [varargout{1:nargout}] = mpimpl(435,varargin{:}); 
        end
        
        function varargout = chebyshevW(varargin)
        %CHEBYSHEVW Chebyshev polynomials of the fourth kind.
        %
        %    Usage syntax is identical to legendreP function.
        %    Please see its help for details.
        %
        %    See also LEGENDREP, CHEBYSHEVT, CHEBYSHEVU, CHEBYSHEVV, HERMITEH.
         
           [varargout{1:nargout}] = mpimpl(436,varargin{:}); 
        end
        
        function varargout = hermiteH(varargin)
        %HERMITEH Hermite polynomials.
        %
        %    Usage syntax is identical to built-in function hermiteH.
        %    See legendreP function for description of the last parameter.
        %
        %    See also LEGENDREP, CHEBYSHEVT, CHEBYSHEVU, CHEBYSHEVV,
        %    CHEBYSHEVW.
         
           [varargout{1:nargout}] = mpimpl(433,varargin{:}); 
        end

        function varargout = gegenbauerC(varargin)
        %GEGENBAUERC Gegenbauer polynomials.
        %    Y = gegenbauerC(N,A,X[,'all']) is the N-th ultraspherical polynomial.
        %
        %    Usage syntax is identical to built-in function gegenbauerC.
        %    See legendreP function for description of the last parameter.
        %
        %    See also LEGENDREP, CHEBYSHEVT, CHEBYSHEVU, CHEBYSHEVV,
        %    CHEBYSHEVW, HERMITEH, JACOBIP.
         
           [varargout{1:nargout}] = mpimpl(434,varargin{:}); 
        end

        function varargout = jacobiP(varargin)
        %JACOBIP Jacobi polynomials.
        %    Y = JACOBIP(N,A,B,X[,'all']) is the N-th Jacobi polynomial.
        %
        %    Usage syntax is identical to built-in function jacobiP.
        %    See legendreP function for description of the last parameter.
        %
        %    See also LEGENDREP, CHEBYSHEVT, CHEBYSHEVU, CHEBYSHEVV,
        %    CHEBYSHEVW, HERMITEH, LAGUERREL, GEGENBAUERC.
         
           [varargout{1:nargout}] = mpimpl(437,varargin{:}); 
        end
        
        function Y = laguerreL(n, a, x, scope)
        %LAGUERREL  Laguerre's L function and Laguerre polynomials.
        %
        %    Usage syntax is identical to built-in function laguerreL.
        %    See legendreP function for description of the last parameter.
        %
        %    See also LEGENDREP, CHEBYSHEVT, CHEBYSHEVU, CHEBYSHEVV,
        %    CHEBYSHEVW, HERMITEH, JACOBIP, GEGENBAUERC.
        
           if nargin < 2
              error(  'MCT:laguerreL:NotEnoughInputs',...
                      'Not enough input arguments.' );
           elseif n < 0
              error('Generalized Laguerre function is not implemented yet. Request it by email: support@advanpix.com');
           else
              if nargin == 2                   % (n,x)
                 Y = mpimpl(438,n,0,a);
              elseif nargin == 3 && ischar(x)  % (n,x,'all')
                 Y = mpimpl(438,n,0,a,x);
              elseif nargin == 3 && ~ischar(x) % (n,a,x)
                 Y = mpimpl(438,n,a,x);
              elseif nargin == 4               % (n,a,x,scope)
                Y = mpimpl(438,n,a,x,scope);                 
              end
           end
        end
        
        function R = zernikeR(varargin)
        %ZERNIKER Radial Zernike polynomial.
        %
        %    Computes radial Zernike polynomials: 
        %
        %        R(n,m,r) = r^m * JacobiP((n-m)/2,0,m,2*r^2-1)
        %
        %    where m and n are nonnegative integers:
        %
        %         m = 1,3,5,...,n, n - odd
        %         m = 0,2,4,...,n, n - even
        %
        %    Optional last argument can be 'all'. Then function computes
        %    polynomial values for all degrees up to n.
        %    
        %    See also LEGENDREP, CHEBYSHEVT, CHEBYSHEVU, CHEBYSHEVV,
        %    CHEBYSHEVW, HERMITEH, JACOBIP, KOORNWINDERK.
         
           R = mpimpl(439,varargin{:}); 
        end

        function varargout = koornwinderK(varargin)   
        %KOORNWINDER Koornwinder orthogonal polynomials.
        %
        %     Computes Koornwinder polynomials:
        %
        %         K(n,j,a,b) = LegendreP(j,(a-b)/(a+b)) * (a+b)^j * JacobiP(n-j,2j+1,0,1-2a-2b)
        %
        %     Optional last argument can be 'all'. Then function computes
        %     polynomial values for all degrees up to n.
        %    
        %    See also LEGENDREP, CHEBYSHEVT, CHEBYSHEVU, CHEBYSHEVV,
        %    CHEBYSHEVW, HERMITEH, JACOBIP, ZERNIKER.
        
           [varargout{1:nargout}] = mpimpl(442,varargin{:});
        end
               
        %% Complex Numbers
        function varargout = abs(varargin)
        %ABS Absolute value.
        %
        %    Usage is identical to built-in function ABS
        %
        %    See also SIGN, ANGLE, UNWRAP, HYPOT.
         
           [varargout{1:nargout}] = mpimpl(500,varargin{:}); 
        end
        
        function varargout = sign(varargin)
        %SIGN Signum function.
        %
        %    Usage is identical to built-in function SIGN
        %
        %    See also ABS.
         
           [varargout{1:nargout}] = mpimpl(501,varargin{:}); 
        end
        
        function varargout = isreal(varargin)
        %ISREAL True for real array.
        %
        %    Usage is identical to built-in function ISREAL
        %
        %    See also REAL, IMAG, COMPLEX, I, J.
         
           [varargout{1:nargout}] = mpimpl(502,varargin{:}); 
        end
        
        function varargout = conj(varargin)
        %CONJ Complex conjugate.
        %
        %    Usage is identical to built-in function CONJ
        %
        %    See also REAL, IMAG, COMPLEX, I, J.
         
           [varargout{1:nargout}] = mpimpl(503,varargin{:}); 
        end
        
        function varargout = cplxpair(varargin)
        %CPLXPAIR Sort numbers into complex conjugate pairs.            
        %
        %    Usage is identical to built-in function ISREAL
        %
        %    See also REAL, IMAG, COMPLEX, I, J.
         
            [varargout{1:nargout}] = mpcplxpair(varargin{:});            
        end
        
        function varargout = unwrap(varargin)
        %UNWRAP Unwrap phase angle.
        %
        %    Usage is identical to built-in function UNWRAP
        %
        %    See also ANGLE, ABS.
         
            [varargout{1:nargout}] = mpunwrap(varargin{:});            
        end
            
        function varargout = angle(varargin)
        %ANGLE Phase angle.
        %
        %    Usage is identical to built-in function ANGLE.
        %
        %    See also I, J, IMAG, CONJ, ABS, REAL, ISREAL.
        
           [varargout{1:nargout}] = mpimpl(504,varargin{:}); 
        end
        
        function r = complex(x,y)
        %COMPLEX Create complex array.
        %
        %    Usage is identical to built-in function COMPLEX.
        %
        %    See also I, J, IMAG, CONJ, ANGLE, ABS, REAL, ISREAL.
        
           if nargin == 1,  r = x;
           else             r = mpimpl(505,x,y);  %#ok<SEPEX>
           end;
        end
        
        function varargout = imag(varargin)
        %IMAG Complex imaginary part.
        %
        %    Usage is identical to built-in function IMAG.
        %
        %    See also REAL, ISREAL, CONJ, ANGLE, ABS.
        
           [varargout{1:nargout}] = mpimpl(506,varargin{:}); 
        end
        
        function varargout = real(varargin)
        %REAL Complex real part.
        %
        %    Usage is identical to built-in function REAL.
        %
        %    See also ISREAL, IMAG, CONJ, ANGLE, ABS.
        
           [varargout{1:nargout}] = mpimpl(507,varargin{:}); 
        end
       
        %% Rounding and Remainder        
        function r = ceil(x)
           r = mpimpl(600,x); 
        end
        
        function r = fix(x)
           r = mpimpl(601,x); 
        end
        
        function r = floor(x)
           r = mpimpl(602,x); 
        end
        
        function varargout = idivide(varargin)
            [varargout{1:nargout}] = mpimpl(603,varargin{:});            
        end
        
        function r = mod(x,y)
           r = mpimpl(604,x,y);
        end
        
        function r = rem(x,y)
           r = mpimpl(605,x,y);
        end
        
        function r = round(x,N,type) %#ok<INUSD>
        %ROUND  rounds towards nearest decimal or integer         
        %
        %    Usage is identical to built-in function ROUND
        %    except third parameter, which is always assumed to be
        %    'decimal'.
        %
        %    See also FLOOR, CEIL, FPRINTF.

           if nargin < 3,          type  = 'decimal'; %#ok<NASGU>
               if nargin < 2,      N     = 0;
                   if nargin < 1
                   error(  'MATLAB:minrhs',...
                           'Not enough input arguments.' );
                   end  
                end
           end
        
           if isscalar(N) && isreal(N) && fix(N) == N
                if N == 0
                  r = mpimpl(606,x);
               else
                   if N > 0
                       r = mpimpl(606,x*10^N)/10^N;
                   else
                       r = mpimpl(606,x/10^(-N))*10^(-N);                       
                   end;
               end
           else
               error(  'MATLAB:round:SecondArg',...
                       'The second input must be a real integer scalar.' );
           end
        end
        
        function varargout = nextabove(varargin)
        %NEXTABOVE Returns the next representable floating-point value towards +Infinity
        %
        %     Example A:
        %
        %     >> nextabove(mp(1))
        %     ans = 
        %         1.00000000000000000000000000000000019
        %     >> nextbelow(mp(1))
        %     ans = 
        %         0.999999999999999999999999999999999904
        %
        %     >> nextabove(mp(-1))
        %     ans = 
        %         -0.999999999999999999999999999999999904
        %     >> nextbelow(mp(-1))
        %     ans = 
        %         -1.00000000000000000000000000000000019        
        %
        %     Example B:
        %
        %     >> mp.Digits(34);
        %     >> mp.FollowMatlabNumericFormat(true);
        %     >> format longE;
        %
        %     >> nextabove(mp(0))
        %     ans = 
        %         6.475175119438025110924438958227647e-4966
        %
        %     >> nextbelow(mp(0))
        %     ans = 
        %         -6.475175119438025110924438958227647e-4966
        %
        %    See also NEXTBELOW

            [varargout{1:nargout}] = mpimpl(607,varargin{:});            
        end

        function varargout = nextbelow(varargin)
        %NEXTBELOW Returns the next representable floating-point value towards -Infinity
        %
        %    See NEXTABOVE for usage examples.

            [varargout{1:nargout}] = mpimpl(608,varargin{:});            
        end
        
        %% Discrete Math        
        function varargout = factorial(varargin)
        %FACTORIAL Factorial of input.
        %
        %    Usage is identical to built-in function FACTORIAL.
        %
        %    See also FACTORIAL

            [varargout{1:nargout}] = mpimpl(700,varargin{:});            
        end

        function varargout = isprime(varargin)
        %ISPRIME Determine which array elements are prime.
        %
        %    Usage is identical to built-in function ISPRIME.
        %
        %    See also ISPRIME

            [varargout{1:nargout}] = mpimpl(701,varargin{:});            
        end
        
        function varargout = primes(varargin)
        %PRIMES Prime numbers less than or equal to input value.
        %
        %    Usage is identical to built-in function PRIMES.
        %
        %    See also PRIMES

            [varargout{1:nargout}] = mpimpl(702,varargin{:});            
        end

        function varargout = factor(varargin)
        %FACTOR Prime factors.
        %
        %    Usage is identical to built-in function FACTOR.
        %
        %    See also FACTOR

            [varargout{1:nargout}] = mpfactor(varargin{:});            
        end

        function varargout = gcd(varargin)
        %GCD Greatest common divisor.
        %
        %    Usage is identical to built-in function GCD.
        %
        %    See also GCD

            [varargout{1:nargout}] = mpimpl(703,varargin{:});            
        end
        
        function varargout = nextprime(varargin)
        %NEXTPRIME Next prime greater than input value.
        %
        %    See also PREVPRIME, PRIMES.

            [varargout{1:nargout}] = mpimpl(704,varargin{:});            
        end
        
        function varargout = prevprime(varargin)
        %PREVPRIME Previous prime less than input value.
        %
        %    See also NEXTPRIME, PRIMES.

            [varargout{1:nargout}] = mpimpl(705,varargin{:});            
        end

        function varargout = lcm(varargin)
        %LCM Least common multiple.
        %
        %    Usage is identical to built-in function LCM.
        %
        %    See also LCM

            [varargout{1:nargout}] = mplcm(varargin{:});            
        end
        
        %% Logical Operations
        function r = islogical(varargin)
        %ISLOGICAL Determine if input is logical array.
        %
        %    Usage is identical to built-in function ISLOGICAL.
        %
        %    See also ISLOGICAL
            
            r = false;
        end
        
        function r = logical(x)
        %LOGICAL Convert numeric values to logical.
        %
        %    Usage is identical to built-in function LOGICAL.
        %
        %    See also LOGICAL
            
            r = logical(double(x));
        end
        
        function r = not(x)
        %NOT Find logical NOT.
        %
        %    Usage is identical to built-in function NOT.
        %
        %    See also NOT
            
            r = mpimpl(803,x);
        end
        
        function r = lt(x, y)
        %LT Less than.
        %
        %    Usage is identical to built-in function LT.
        %
        %    See also LT
            
            r = mpimpl(804,x,y);
        end
        
        function r = gt(x, y)
        %GT Greater than.
        %
        %    Usage is identical to built-in function GT.
        %
        %    See also GT
            
            r = mpimpl(805,x,y);
        end
        
        function r = le(x, y)
        %LE Less than or equal.
        %
        %    Usage is identical to built-in function LE.
        %
        %    See also LE
            
            r = mpimpl(806,x,y);
        end
        
        function r = ge(x, y)
        %GE Greater than or equal.
        %
        %    Usage is identical to built-in function GE.
        %
        %    See also GE
            
            r = mpimpl(807,x,y);
        end
        
        function r = ne(x, y)
        %NE Not equal.
        %
        %    Usage is identical to built-in function NE.
        %
        %    See also NE
            
            r = mpimpl(808,x,y);
        end
        
        function r = eq(x, y)
        %EQ Equal.
        %
        %    Usage is identical to built-in function EQ.
        %
        %    See also EQ
            
            r = mpimpl(809,x,y);
        end
        
        function r = and(x, y)
        %AND Logical and.
        %
        %    Usage is identical to built-in function AND.
        %
        %    See also and
        
            r = mpimpl(810,x,y);
        end
        
        function r = or(x, y)
        %OR Logical OR.
        %
        %    Usage is identical to built-in function OR.
        %
        %    See also OR.
            
            r = mpimpl(811,x,y);
        end

        function r = xor(x, y)
        %XOR Logical XOR.
        %
        %    Usage is identical to built-in function XOR.
        %
        %    See also XOR.
            
            r = mpimpl(812,x,y);
        end
        
        function varargout = isequal(varargin)
        %ISEQUAL Determine array equality.
        %
        %    Usage is identical to built-in function ISEQUAL.
        %
        %    See also ISEQUAL.
            
            [varargout{1:nargout}] = mpimpl(813, varargin{:});
        end
        
        function varargout = isequaln(varargin)
        %ISEQUALN Determine array equality, treating NaN values as equal.
        %
        %    Usage is identical to built-in function ISEQUALN.
        %
        %    See also ISEQUALN.
            
            [varargout{1:nargout}] = mpimpl(814, varargin{:});
        end

        function r = all(A, dim)
        %ALL True if all elements of a vector are nonzero.
        %
        %    Usage is identical to built-in function ALL.
        %
        %    See also all
        
            if nargin == 1,     r = all(A~=0);          
            elseif nargin == 2, r = all(A~=0, dim); end 
        end
        
        function r = any(A, dim)
        %ANY True if any element of a vector is a nonzero number or is logical 1 (TRUE).  
        %    ANY ignores entries that are NaN (Not a Number).
        %
        %    Usage is identical to built-in function ANY.
        %
        %    See also any
        
            if nargin == 1,     r = any(A~=0);          
            elseif nargin == 2, r = any(A~=0, dim); end 
        end
        
         %% Specialized Matrices - with 'mp' arguments
        function varargout = compan(varargin)
            [varargout{1:nargout}] = mpcompan(varargin{:});            
        end
        
        function varargout = hankel(varargin)
            [varargout{1:nargout}] = mphankel(varargin{:});            
        end
        
        function varargout = vander(varargin)
            [varargout{1:nargout}] = mpvander(varargin{:});            
        end
         
        function varargout = toeplitz(varargin)
            [varargout{1:nargout}] = mptoeplitz(varargin{:});            
        end
        
        %% Linear Algebra
        function varargout = lu(varargin)
        %LU LU factorization.
        %
        %    Usage syntax is identical to built-in function LU.
        %
        %    See also CHOL, ILU, QR.
            
            [varargout{1:nargout}] = mpimpl(1000, varargin{:});
        end
       
        function varargout = svd(varargin)
        %SVD Singular value decomposition.
        %
        %    Usage syntax is identical to built-in function SVD.
        %
        %    See also SVDS, GSVD.
            
            [varargout{1:nargout}] = mpimpl(1001, varargin{:});
        end

        function varargout = gsvd(varargin)
        %GSVD Generalized Singular Value Decomposition.
        %
        %    Usage syntax is identical to built-in function GSVD.
        %
        %    Please note, the function is partially implemented using MATLAB language 
        %    and hence might be slow for large matrices.
        %    Let us know if this is the case - we will prioritize 
        %    its implementation in C++.
        %
        %    See also SVD.
            
            [varargout{1:nargout}] = mpgsvd(varargin{:});
        end
        
        function varargout = csd(varargin)
        %CSD Cosine Sine Decomposition.
        %
        %  Full (2 by 2):
        %    
        %    [C,S,U1,V1,U2,V2] = CSD(X11,X12,X21,X22)
        %    
        %                                           [  I  0  0 |  0  0  0 ]
        %               Q    M-Q                    [  0  C  0 |  0 -S  0 ]
        %            [ X11 | X12 ] P    [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]^T
        %        X = [-----------]    = [---------] [---------------------] [---------]   
        %            [ X21 | X22 ] M-P  [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
        %                                           [  0  S  0 |  0  C  0 ]
        %                                           [  0  0  I |  0  0  0 ]
        %  Thin (2 by 1):
        %    
        %    [C,S,U1,V1,U2] = CSD(X11,X21)
        %    
        %                                        [  I  0  0 ]
        %                 Q                      [  0  C  0 ]
        %              [ X11 ] P     [ U1 |    ] [  0  0  0 ]
        %          X = [-----]     = [---------] [----------] V1^T
        %              [ X21 ] M-P   [    | U2 ] [  0  0  0 ]
        %                                        [  0  S  0 ]
        %                                        [  0  0  I ]
        %
        %
        %    See also SVD, GSVD.
            
            [varargout{1:nargout}] = mpimpl(1041, varargin{:});
        end

        function varargout = qr(varargin)
        %QR  Orthogonal-triangular decomposition.
        %
        %    Usage syntax is identical to built-in function QR.
        %
        %    See also LU, NULL, ORTH, QRDELETE, QRINSERT, QRUPDATE.
            
            [varargout{1:nargout}] = mpimpl(1015, varargin{:});
        end
        
        function varargout = qrdelete(varargin)
        %QRDELETE Delete a column or row from QR factorization.
        %   [Q1,R1] = QRDELETE(Q,R,J) returns the QR factorization of the matrix A1,
        %   where A1 is A with the column A(:,J) removed and [Q,R] = QR(A) is the QR
        %   factorization of A. Matrices Q and R can also be generated by 
        %   the "economy size" QR factorization [Q,R] = QR(A,0). 
        %
        %   Usage syntax is identical to built-in function QRDELETE.
        %
        %   See also QR, QRINSERT, QRUPDATE, PLANEROT.
            
            [varargout{1:nargout}] = mpqrdelete(varargin{:});
        end
        
        function varargout = qrinsert(varargin)
        %QRINSERT Insert a column or row into QR factorization.
        %   [Q1,R1] = QRINSERT(Q,R,J,X) returns the QR factorization of the matrix A1,
        %   where A1 is A=Q*R with an extra column, X, inserted before A(:,J). If A has
        %   N columns and J = N+1, then X is inserted after the last column of A.
        %
        %   Usage syntax is identical to built-in function QRINSERT.
        %
        %   See also QR, QRDELETE, QRUPDATE, PLANEROT.
            
            [varargout{1:nargout}] = mpqrinsert(varargin{:});
        end

        function varargout = qrupdate(varargin)
        %QRUPDATE Rank 1 update to QR factorization.
        %   If [Q,R] = QR(A) is the original QR factorization of A, then
        %   [Q1,R1] = QRUPDATE(Q,R,U,V) returns the QR factorization of A + U*V',
        %   where U and V are column vectors of appropriate lengths.
        %
        %   Usage syntax is identical to built-in function QRUPDATE.
        %
        %   See also QR, QRDELETE, QRINSERT, PLANEROT.
            
            [varargout{1:nargout}] = mpimpl(1039,varargin{:});
        end
        
        function varargout = planerot(varargin)
        %PLANEROT Givens plane rotation.
        %   [G,Y] = PLANEROT(X), where X is a 2-component column vector,
        %   returns a 2-by-2 orthogonal matrix G so that Y = G*X has Y(2) = 0.
        %
        %    Usage syntax is identical to built-in function LDL.
        %
        %   See also QRINSERT, QRDELETE.
            
            [varargout{1:nargout}] = mpplanerot(varargin{:});
        end
        
        function varargout = chol(varargin)
        %CHOL Cholesky factorization.
        %
        %    Usage syntax is identical to built-in function CHOL.
        %
        %    See also CHOLUPDATE, ICHOL, LDL, LU.
            
            [varargout{1:nargout}] = mpimpl(1017, varargin{:});
        end

        function varargout = cholupdate(varargin)
        %CHOLUPDATE Rank 1 update to Cholesky factorization.
        %   If R = CHOL(A) is the original Cholesky factorization of A, then
        %   R1 = CHOLUPDATE(R,X) returns the upper triangular Cholesky factor of A + X*X',
        %   where X is a column vector of appropriate length.  CHOLUPDATE uses only the
        %   diagonal and upper triangle of R.  The lower triangle of R is ignored.
        %
        %   Usage syntax is identical to built-in function CHOLUPDATE.
        %
        %   See also CHOL..
            
            [varargout{1:nargout}] = mpimpl(1040,varargin{:});
        end
        
        function varargout = ldl(varargin)
        %LDL Block LDL' factorization for Hermitian indefinite matrices.
        %
        %    Usage syntax is identical to built-in function LDL.
        %
        %   See also CHOL, LU, QR.
            
            [varargout{1:nargout}] = mpimpl(1034, varargin{:});
        end
        
        function varargout = eig(varargin)
        %EIG Eigenvalues and eigenvectors.
        %
        %    Usage syntax is identical to built-in function EIG.
        %
        %    In case of symmetric/Hermitian and real tridiagonal standard eigenproblems
        %    toolbox accepts second optional parameter - eig(A,method).
        % 
        %    Values for the parameter:
        %    'mr' - Multiple Relatively Robust Representations (MRRR).
        %    'dc' - Divide and conquer algorithm.
        %    'qr' - Implicit QL or QR method (Pal-Walker-Kahan for eigenvalues only).
        %     
        %    By default, method='mr' if no eigenvectors required, 'qr' otherwise.
        %
        %    In case of very ill-conditioned eigenvectors, QR might
        %    give more accurate results working with the same precision.
        %    MRRR and D&C might require increase in precision to provide
        %    the same accuracy.
        % 
        %    More details of implementation can be found in the article:
        %    http://www.advanpix.com/2016/10/20/architecture-of-eigenproblem-solver/
        %
        %    See also CONDEIG, EIGS, ORDEIG.
            
            [varargout{1:nargout}] = mpimpl(1016, varargin{:});
        end

        function varargout = eigs(varargin)
        %EIGS  Find a few eigenvalues and eigenvectors of a matrix using 
        %      Krylov-Schur decomposition implemented using arbitrary
        %      precision precision.
        %   
        %   Working modes and options are equivalent to classic EIGS. The only
        %   difference is that we are using Krylov-Schur decomposition instead of Arnoldi (ARPACK).
        %   Krylov-Schur is reportedly more numerically stable, and combined with
        %   extended precision allows solving sensitive or ill-conditioned
        %   eigenproblems. See MPEIGS.M source code for more details and references.
        %
        %   Algorithm is implemented using MATLAB language and thus it is (very)slow.
        %   After short period tests and feedback from users we will move it to C++.
        %   Estimated speed-up will be ~25-100 times. Stay tuned.
        %
        %   See also EIG, SVDS, FUNCTION_HANDLE.
            
            [varargout{1:nargout}] = mpeigs(varargin{:});
        end
        
        function varargout = condeig(varargin)
        %CONDEIG Condition number with respect to eigenvalues.
        %
        %    Usage syntax is identical to built-in function CONDEIG.
        %
        %    See also COND, EIG, ORDEIG.
            
            [varargout{1:nargout}] = mpcondeig(varargin{:});
        end
        
        function varargout = qz(varargin)
        %QZ Factorization for generalized eigenvalues.
        %
        %    Usage syntax is identical to built-in function QZ.
        %
        %    See also ORDQZ, ORDEIG, EIG.
                    
            [varargout{1:nargout}] = mpimpl(1026, varargin{:});
        end

        function varargout = ordqz(varargin)
        %ORDQZ Reorder eigenvalues in QZ factorization.
        %
        %    Usage syntax is identical to built-in function ORDQZ.
        %
        %    See also QZ, ORDEIG, EIG.
                    
            [varargout{1:nargout}] = mpimpl(1027, varargin{:});
        end

        function varargout = ordeig(varargin)
        %ORDEIG Eigenvalues of quasi-triangular matrices.
        %
        %    Usage syntax is identical to built-in function ORDEIG.
        %
        %    See also QZ, ORDQZ, SCHUR, ORDSCHUR.
                    
            [varargout{1:nargout}] = mpimpl(1028, varargin{:});
        end
        
        function varargout = schur(varargin)
        %SCHUR Schur decomposition.
        %
        %    Usage syntax is identical to built-in function SCHUR.
        %
        %    See also ORDSCHUR, QZ, ORDEIG, ORDQZ, EIG.
            
            [varargout{1:nargout}] = mpimpl(1018, varargin{:});
        end
        
        function varargout = ordschur(varargin)
        %ORSCHUR Reorder eigenvalues in Schur factorization.
        %
        %    Usage syntax is identical to built-in function ORDSCHUR.
        %
        %    In addition to standard KEYWORD options we also support:
        %
        %       'ref' Real eigenvalues first (top-left conner)
        %       'cef' Complex eigenvalues first (top-left conner)
        %
        %    See also SCHUR, QZ, ORDEIG, ORDQZ, EIG.
            
            [varargout{1:nargout}] = mpimpl(1022, varargin{:});
        end
        
        function varargout = hess(varargin)
        %HESS Hessenberg form of matrix.
        %
        %    Usage syntax is identical to built-in function HESS.
        %
        %    See also EIG, QZ, SCHUR.
            
            [varargout{1:nargout}] = mpimpl(1021, varargin{:});
        end
        
        function varargout = trisylv(varargin)
        %TRISYLV Solve triangular Sylvester equation.
        %
        %   DEPRECATED. Please use SYLVESTER_TRI instead.
        %
        %   X = TRISYLV(A,B,C) solves the Sylvester equation
        %   A*X + X*B = C, where A and B are square upper triangular matrices.
        %   
        %   See also SYLVESTER_TRI, SYLVESTER.        
        
            [varargout{1:nargout}] = mpimpl(1023, varargin{:});
        end
        
        function varargout = sort(varargin)
        %SORT Sort in ascending or descending order.
        %
        %    Usage syntax is identical to built-in function SORT.
        %
        %    See also ISSORTED, SORTROWS, MIN, MAX, MEAN, MEDIAN, UNIQUE.
            
            [varargout{1:nargout}] = mpimpl(1100, varargin{:});
        end
        
        function varargout = sortrows(varargin)
        %SORTROWS Sort rows in ascending order.
        %
        %    Usage syntax is identical to built-in function SORTROWS.
        %
        %    See also SORT, ISSORTED.
            
            [varargout{1:nargout}] = mpsortrows(varargin{:});
        end

        function varargout = det(varargin)
        %DET Determinant.
        %
        %    Usage syntax is identical to built-in function DET.
        %
        %    See also COND.
            
            [varargout{1:nargout}] = mpimpl(1002,varargin{:});
        end
        
        function varargout = inv(varargin)
        %INV Matrix inverse.
        %
        %    Usage syntax is identical to built-in function INV.
        %
        %    See also SLASH, PINV, COND, CONDEST.
            
            [varargout{1:nargout}] = mpimpl(1003,varargin{:});
        end
        
        function [X,R] = linsolve(A,b,opts)
        %LINSOLVE Solve linear system A*X=B.
        %
        %    Toolbox has mature MLDIVIDE enabled with automatic
        %    detection of matrix properties and set of optimal algorithms
        %    for each matrix type. Hence LINSOLVE just calls MLDIVIDE.
        %
        %    Cost of transpose is negligible compared to complexity of
        %    direct solvers especilly in case of extended precision.  
        %
        %    Usage syntax is identical to built-in function LINSOLVE.
        %
        %    See also MLDIVIDE, SLASH.
            
            if nargin < 2, error(message('MATLAB:narginchk:notEnoughInputs'));  end;
            if ~(isnumeric(A) && isnumeric(b)), error('First and second arguments must be numeric.');  end;
            
            transa = false;
            if nargin == 3 && isfield(opts,'TRANSA'), transa = opts.TRANSA; end;

            % Let MLDIVIDE do its job - detect matrix type and apply the best
            % matching algorithm to solve it.
            
            if nargout > 1
                w1 = warning('off','MATLAB:nearlySingularMatrix');
                w2 = warning('off','MATLAB:singularMatrix');
                w3 = warning('off','MATLAB:rankDeficientMatrix');

                % Force mldivide to compute rcond, estimate rank (RRQR), etc.
                if transa, [X,R] = mpimpl(40,A',b);  
                else       [X,R] = mpimpl(40,A, b); end; %#ok<SEPEX>

                warning(w1);
                warning(w2);
                warning(w3);                
            else
                % Only solution is required
                if transa, X = A'\b;  
                else       X = A \b; end; %#ok<SEPEX>
            end
        end

        function varargout = issymmetric(varargin)
        %ISSYMMETRIC Determine whether a matrix is real or complex symmetric.
        %
        %    Usage syntax is identical to built-in function ISSYMMETRIC.
        %
        %    See also ISHERMITIAN.
            
            [varargout{1:nargout}] = mpimpl(1029,varargin{:});
        end

        function varargout = ishermitian(varargin)
        %ISHERMITIAN Determine whether a matrix is real symmetric or complex Hermitian.
        %
        %    Usage syntax is identical to built-in function ISHERMITIAN.
        %
        %    See also ISSYMMETRIC.
            
            [varargout{1:nargout}] = mpimpl(1030,varargin{:});
        end
        
        function varargout = istril(varargin)
        %ISTRIL Determine whether a matrix is lower triangular.
        %
        %    Usage syntax is identical to built-in function ISTRIL.
        %
        %    See also ISDIAG, ISTRIU, DIAG, TRIL, TRIU.
            
            [varargout{1:nargout}] = mpimpl(1031,varargin{:});
        end

        function varargout = istriu(varargin)
        %ISTRIU Determine whether a matrix is upper triangular.
        %
        %    Usage syntax is identical to built-in function ISTRIU.
        %
        %    See also ISDIAG, ISTRIL, DIAG, TRIL, TRIU.
            
            [varargout{1:nargout}] = mpimpl(1032,varargin{:});
        end

        function varargout = isdiag(varargin)
        %ISDIAG Determine whether a matrix is diagonal.
        %
        %    Usage syntax is identical to built-in function ISDIAG.
        %
        %    See also ISTRIL, ISTRIU, DIAG, TRIL, TRIU.
            
            [varargout{1:nargout}] = mpimpl(1033,varargin{:});
        end

        function varargout = isbanded(varargin)
        %ISBANDED Determine whether a matrix is within a specific bandwidth.
        %
        %    Usage syntax is identical to built-in function ISBANDED.
        %
        %    See also BANDWIDTH, ISDIAG, ISTRIL, ISTRIU.
            
            [varargout{1:nargout}] = mpimpl(1035,varargin{:});
        end
        
        function varargout = isschur(varargin)
        %ISSCHUR Determine whether T has Schur form conformed with LAPACK
        %   ISSCHUR(T) returns true if matrix T is upper quasi-triangular for real
        %   matrix, or upper triangular for complex matrix.
        %
        %   ISSCHUR(T, 'complex') returns true if T is upper triangular.
        %
        %   ISSCHUR(T, 'real') is same as ISSCHUR(T).
        %   
        %   Usage syntax is identical to built-in function ISSCHUR 
        %   located in matlab.internal.math.isschur
            
            [varargout{1:nargout}] = mpimpl(1036,varargin{:});
        end

        function varargout = bandwidth(varargin)
        %[LOWER,UPPER] = BANDWIDTH(X) returns the lower bandwidth LOWER, and the upper 
        %   bandwidth UPPER of matrix X.
        %   
        %   LOWER = BANDWIDTH(X,'lower') returns the lower bandwidth LOWER of matrix X.
        % 
        %   UPPER = BANDWIDTH(X,'upper') returns the upper bandwidth UPPER of matrix X.
        %
        %   Usage syntax is identical to built-in function BANDWIDTH 
        % 
        %   See also ISBANDED, ISDIAG, ISTRIL, ISTRIU.
        %   
        
            [varargout{1:nargout}] = mpimpl(1038,varargin{:});
        end
        
        function varargout = pinv(varargin)
        %PINV Pseudoinverse.
        %
        %    Usage syntax is identical to built-in function PINV.
        %
        %    See also RANK.
            
            [varargout{1:nargout}] = mpimpl(1004,varargin{:});
        end
        
        function varargout = null(varargin)
        %NULL  Null space.            
        %
        %    Usage is identical to built-in function NULL.
        %
        %    See also SVD, ORTH, RANK, RREF.
       
               [varargout{1:nargout}] = mpnull(varargin{:});
        end
        
        function varargout = orth(varargin)
        %ORTH  Orthonormal basis for range of matrix
        %
        %    Usage is identical to built-in function ORTH.
        %
        %    See also SVD, RANK, NULL.
        
               [varargout{1:nargout}] = mporth(varargin{:});
        end
        
        function r = rank(A,tol)
        %RANK Matrix rank.
        %
        %    Usage syntax is identical to built-in function RANK.
        %
        %    See also SVD, COND.
            
            s = svd(A);
            if nargin==1, tol = max(size(A)) * eps(max(s)); end
            r = sum(s > tol);
        end

        function varargout = subspace(varargin)
        %SUBSPACE Angle between subspaces.
        %   SUBSPACE(A,B) finds the angle between two subspaces specified by the
        %   columns of A and B.
        %
        %    Usage syntax is identical to built-in function SUBSPACE.
        %
        %    See also ORTH, RREF.
        
            [varargout{1:nargout}] = mpsubspace(varargin{:});
        end
        
        function varargout = trace(varargin)
        %TRACE Sum of diagonal elements.
        %
        %    Usage syntax is identical to built-in function TRACE.
        %
        %    See also RANK, COND.
        
            [varargout{1:nargout}] = mpimpl(1006,varargin{:});
        end
      
        function varargout = norm(varargin)
        %NORM Matrix or vector norm.
        %
        %    Usage syntax is identical to built-in function NORM.
        %
        %    See also COND, RCOND, CONDEST, NORMEST, HYPOT.
        
            [varargout{1:nargout}] = mpimpl(1007,varargin{:});
        end
        
        function c = cond(A, p)
        %COND Condition number with respect to inversion.  
        %
        %    Usage syntax is identical to built-in function COND.
        %
        %    See also RCOND, CONDEST, CONDEIG, NORM, NORMEST.
        
            if nargin < 2, p = 2; end;
            
            if issparse(A)
                % Use 'double' precision for sparse matrices, since
                % we do not have full-featured sparse LU yet.
                warning(message('MATLAB:cond:SparseNotSupported'))
                c = condest(double(A));
            else
                
                % Matrix must be square if p!=2
                [m, n] = size(A);
                if m~=n && ~isequal(p,2)
                    error(message('MATLAB:cond:normMismatchSizeA')); 
                end
                
                % Standard code for condition number estimation
                if p == 2
                    s = svd(A);
                    if any(s == 0)
                        c = mp('Inf');
                    else
                        c = max(s)./min(s);
                        if isempty(c)
                            c = mp('0');
                        end
                    end
                else
                    c = norm(A,p) * norm(inv(A),p);
                end
             end
        end
        
        function varargout = rcond(varargin)
        %RCOND Reciprocal condition number.
        %
        %    Usage is identical to built-in function RCOND.
        %
        %    See also COND, CONDEST, NORM, NORMEST, RANK, SVD.
            
            [varargout{1:nargout}] = mpimpl(1024,varargin{:});
        end
        
        function varargout = rref(varargin)
        %RREF   Reduced row echelon form.
        %
        %   Usage is identical to built-in function RREF.
        %
        %   See also RANK, ORTH, NULL, QR, SVD.
        
               [varargout{1:nargout}] = mprref(varargin{:});
        end
        
        function varargout = balance(varargin)
        %BALANCE Diagonal scaling to improve eigenvalue accuracy.
        %
        %    Usage is identical to built-in function BALANCE.
        %
        %    See also BALANCE
        
            [varargout{1:nargout}] = mpimpl(1025,varargin{:});
        end
        
        function varargout = diag(varargin)
               [varargout{1:nargout}] = mpimpl(1008,varargin{:});
        end
        
        function varargout = triu(varargin)
               [varargout{1:nargout}] = mpimpl(1009,varargin{:});
        end

        function varargout = tril(varargin)
               [varargout{1:nargout}] = mpimpl(1010,varargin{:});
        end
       
        function varargout = sum(varargin)
               [varargout{1:nargout}] = mpimpl(1011,varargin{:});
        end
       
        function varargout = prod(varargin)
               [varargout{1:nargout}] = mpimpl(1012,varargin{:});
        end

        function varargout = cumprod(varargin)
        %CUMPROD Cumulative product of elements.
        %   Y = CUMPROD(X) computes the cumulative product along the first non-singleton
        %   dimension of X. Y is the same size as X.
        %
        %    Usage is identical to built-in function CUMPROD.
        %
        %    See also PROD, CUMSUM, CUMMIN, CUMMAX.
            
               [varargout{1:nargout}] = mpimpl(3000,varargin{:});
        end
        
        function varargout = cumsum(varargin)
        %CUMSUM Cumulative sum of elements.
        %   Y = CUMSUM(X) computes the cumulative sum along the first non-singleton
        %   dimension of X. Y is the same size as X.
        %
        %   Usage is identical to built-in function CUMSUM.
        %
        %   See also SUM, MOVSUM, CUMPROD, CUMMIN, CUMMAX.
            
               [varargout{1:nargout}] = mpimpl(3001,varargin{:});
        end
        
        function varargout = cummin(varargin)
        %CUMMIN Cumulative smallest component.
        %   Y = CUMMIN(X) computes the cumulative smallest component of X along
        %   the first non-singleton dimension of X. Y is the same size as X.
        %
        %   Usage is identical to built-in function CUMMIN.
        %
        %   See also MIN, MOVMIN, CUMMAX, CUMSUM, CUMPROD.
            
               [varargout{1:nargout}] = mpimpl(3005,varargin{:});
        end
        
        function varargout = cummax(varargin)
        %CUMMAX Cumulative largest component.
        %   Y = CUMMAX(X) computes the cumulative largest component of X along
        %   the first non-singleton dimension of X. Y is the same size as X.
        %
        %   Usage is identical to built-in function CUMMAX.
        %
        %   See also MIN, MOVMIN, CUMMIN, CUMSUM, CUMPROD.
            
               [varargout{1:nargout}] = mpimpl(3006,varargin{:});
        end
        
        function varargout = dot(varargin)
        %DOT  Vector dot product.
        %   C = DOT(A,B) returns the scalar product of the vectors A and B.
        %   A and B must be vectors of the same length.  When A and B are both
        %   column vectors, DOT(A,B) is the same as A'*B.
        %
        %   Usage is identical to built-in function DOT.
        %
        %   See also CROSS, SUM, KRON.
            
               [varargout{1:nargout}] = mpdot(varargin{:});
        end
   
        function varargout = cross(varargin)
        %CROSS  Vector cross product.
        %   C = CROSS(A,B) returns the cross product of the vectors
        %   A and B.  That is, C = A x B.  A and B must be 3 element
        %   vectors.
        %
        %   Usage is identical to built-in function CROSS.
        %
        %   See also DOT, KRON.
            
               [varargout{1:nargout}] = mpcross(varargin{:});
        end
        
        function varargout = kron(varargin)
        %KRON Kronecker tensor product
        %
        %    Usage is identical to built-in function KRON.
        %
        %    See also CROSS, DOT, HANKEL, TOEPLITZ.
            
               [varargout{1:nargout}] = mpimpl(3004,varargin{:});
        end
        
        %% Basic Statistics 
        function varargout = mean(varargin)
        %MEAN Average or mean value.
        %
        %    Usage is identical to built-in function MEAN.
        %
        %    See also MEDIAN, STD, MIN, MAX, VAR, COV, MODE.
        
               [varargout{1:nargout}] = mpimpl(2004,varargin{:});
        end
        
        function varargout = std(varargin)
        %STD Standard deviation.
        %
        %    Usage is identical to built-in function STD.
        %
        %    See also COV, MEAN, VAR, MEDIAN, CORRCOEF.
            
               [varargout{1:nargout}] = mpimpl(2005,varargin{:});
        end
        
        function varargout = max(varargin)
        %MAX Largest component.
        %
        %    Usage is identical to built-in function MAX.
        %
        %    See also MIN, CUMMIN, MEDIAN, MEAN, SORT.
            
               [varargout{1:nargout}] = mpimpl(1013,varargin{:});
        end
    
        function varargout = min(varargin)
        %MIN Smallest component.
        %
        %    Usage is identical to built-in function MIN.
        %
        %    See also MAX, CUMMIN, MEDIAN, MEAN, SORT.
            
               [varargout{1:nargout}] = mpimpl(1014,varargin{:});
        end

        function varargout = del2(varargin)
        %DEL2 Discrete Laplacian.
        %
        %    Usage is identical to built-in function DEL2.
        %
        %    See also GRADIENT, DIFF.
        
               [varargout{1:nargout}] = mpdel2(varargin{:});
        end
        
        function varargout = diff(varargin)
        %DIFF Difference and approximate derivative.
        %
        %    Usage is identical to built-in function DIFF.
        %
        %    See also GRADIENT, SUM, PROD.
        
               [varargout{1:nargout}] = mpdiff(varargin{:});
        end
        
        function varargout = gradient(varargin)
        %GRADIENT Approximate gradient.
        %
        %    Usage is identical to built-in function GRADIENT.
        %
        %    See also DIFF, DEL2.
        
               [varargout{1:nargout}] = mpgradient(varargin{:});
        end
        
        %% Sparse matrices
        function r = issparse(x)
        %ISSPARSE Determine whether input is sparse
        %
        %    Usage is identical to built-in function ISSPARSE.
        %
        %    See also ISSPARSE.
        
           r = mpimpl(8501, x); 
        end
        
        function r = nnz(x)
        %NNZ Number of nonzero matrix elements
        %
        %    Usage is identical to built-in function NNZ.
        %
        %    See also NNZ.
        
           r = mpimpl(8502, x); 
        end

        function r = nonzeros(x)
        %NONZEROS Nonzero matrix elements
        %
        %    Usage is identical to built-in function NONZEROS.
        %
        %    See also NONZEROS.
        
           r = mpimpl(8503, x); 
        end
        
        function r = nzmax(x)
        %NZMAX Amount of storage allocated for nonzero matrix elements
        %
        %    Usage is identical to built-in function NZMAX.
        %
        %    See also NZMAX.
           
           r = mpimpl(8504, x); 
        end
        
        function varargout = full(varargin)
        %FULL Convert sparse matrix to full matrix
        %
        %    Usage is identical to built-in function FULL.
        %
        %    See also FULL
            
            [varargout{1:nargout}] = mpimpl(8505, varargin{:});
        end
        
        function varargout = sparse(varargin)
        %SPARSE Create sparse matrix
        %
        %    Usage is identical to built-in function SPARSE.
        %
        %    See also SPARSE
            
            [varargout{1:nargout}] = mpimpl(8508, varargin{:});
        end

        function varargout = colamd(S,varargin)
        %COLAMD Column approximate minimum degree permutation.
        %
        %    Usage is identical to built-in function COLAMD.
        %
        %    See also COLAMD

            [row,col] = find(S);
            [varargout{1:nargout}] = colamd(sparse(row,col,1), varargin{:});
        end
        
        function varargout = amd(S,varargin)
        %AMD  Approximate minimum degree permutation.
        %
        %    Usage is identical to built-in function AMD.
        %
        %    See also AMD

            [row,col] = find(S);
            [varargout{1:nargout}] = amd(sparse(row,col,1), varargin{:});
        end

        function varargout = symamd(S,varargin)
        %SYMAMD Symmetric approximate minimum degree permutation.
        %
        %    Usage is identical to built-in function SYMAMD.
        %
        %    See also SYMAMD

            [row,col] = find(S);
            [varargout{1:nargout}] = symamd(sparse(row,col,1), varargin{:});
        end

        function varargout = symrcm(S,varargin)
        %SYMRCM Symmetric reverse Cuthill-McKee permutation.
        %
        %    Usage is identical to built-in function SYMRCM.
        %
        %    See also SYMRCM

            [row,col] = find(S);
            [varargout{1:nargout}] = symrcm(sparse(row,col,1), varargin{:});
        end
        
        function r = spones(x)
        %SPONES Replace nonzero sparse matrix elements with ones
        %
        %    Usage is identical to built-in function SPONES.
        %
        %    See also SPONES.
           
           r = mpimpl(8506, x); 
        end
        
        function r = spfun(fun, s)
        %SPFUN Apply function to nonzero sparse matrix elements
        %
        %    Usage is identical to built-in function SPFUN.
        %
        %    See also SPFUN.
           
            if ~ismatrix(s), error(message('MCT:spfun:ndInput')); end
            [i,j,x] = find(s);
            [m,n] = size(s);
            r = sparse(i,j,feval(fun,x),m,n);
        end
        
        function varargout = find(varargin)
        %FIND Find indices and values of nonzero elements
        %
        %    Usage is identical to built-in function FIND.
        %
        %    See also FIND
            
            [varargout{1:nargout}] = mpimpl(8000, varargin{:});
        end
        
        function varargout = spdiags(varargin)
        %SPDIAGS Sparse matrix formed from diagonals.
        %
        %    Usage is identical to built-in function SPDIAGS.
        %
        %    See also SPDIAGS
            
            [varargout{1:nargout}] = mpspdiags(varargin{:});
        end
        
        %% System functions
        function  classname = superiorfloat(varargin)  
            classname = 'mp'; 
        end

        %% Numerical Integration
        function varargout = integral(varargin)
        %INTEGRAL  Numerically evaluate integral.
        %
        %    Usage is identical to built-in function INTEGRAL.
        %    except that we also provide error bound as second output:
        %
        %    [q,errbnd] = integral(...)
        %
        %    See also QUADGK, INTEGRAL2, INTEGRAL3, TRAPZ, FUNCTION_HANDLE
            
            [varargout{1:nargout}] = mpintegral(varargin{:});
        end

        function varargout = integral2(varargin)
        %INTEGRAL2  Numerically evaluate double integral.
        %
        %    Usage is identical to built-in function INTEGRAL2,
        %    except that we also provide error bound as second output:
        %
        %    [q,errbnd] = integral2(...)
        %
        %    See also QUADGK, INTEGRAL, INTEGRAL3, TRAPZ, FUNCTION_HANDLE
            
            [varargout{1:nargout}] = mpintegral2(varargin{:});
        end
        
        function varargout = integral3(varargin)
        %INTEGRAL3  Numerically evaluate triple integral.
        %
        %    Usage is identical to built-in function INTEGRAL3,
        %    except that we also provide error bound as second output:
        %
        %    [q,errbnd] = integral2(...)
        %
        %    See also QUADGK, INTEGRAL, INTEGRAL2, TRAPZ, FUNCTION_HANDLE
            
            [varargout{1:nargout}] = mpintegral3(varargin{:});
        end
        
        function varargout = quad2d(varargin)
        %QUAD2D    Numerically evaluate double integral over a planar region.
        %
        %    Usage is identical to built-in function QUAD2D.
        %
        %    See also QUADGK, INTEGRAL, INTEGRAL2, TRAPZ, FUNCTION_HANDLE
            
            [varargout{1:nargout}] = mpquad2d(varargin{:});
        end
        
        function varargout = trapz(varargin)
        %TRAPZ  Trapezoidal numerical integration.
        %
        %    Usage is identical to built-in function TRAPZ.
        %
        %    See also SUM, CUMSUM, CUMTRAPZ, INTEGRAL.
            
            [varargout{1:nargout}] = mptrapz(varargin{:});
        end

        function varargout = cumtrapz(varargin)
        %CUMTRAPZ Cumulative trapezoidal numerical integration.
        %
        %    Usage is identical to built-in function CUMTRAPZ.
        %
        %    See also CUMSUM, TRAPZ, INTEGRAL.
            
            [varargout{1:nargout}] = mpcumtrapz(varargin{:});
        end
        
        function varargout = quadgk(varargin)
        %QUADGK Numerically evaluate integral, adaptive Gauss-Kronrod quadrature
        %
        %    Usage is identical to built-in function QUADGK.
        %
        %    See also INTEGRAL
            
            [varargout{1:nargout}] = mpquadgk(varargin{:});
        end
        
        function varargout = integralCalc(varargin)
        %INTEGRALCALC Helper function for INTEGRAL and INTEGRAL2.
        %
        %    This routine is used by INTEGRAL, INTEGRAL2 and INTEGRAL3. 
        %    Usually you do not need to call it explicitly.
        %
        %    See also INTEGRAL, INTEGRAL2, INTEGRAL3.
            
            [varargout{1:nargout}] = mpintegralCalc(varargin{:});
        end
        
        function varargout = integral2Calc(varargin)
        %INTEGRAL2CALC Helper function for INTEGRAL2.
        %
        %    This routine is used by INTEGRAL2. 
        %    Usually you do not need to call it explicitly.
        %
        %    See also INTEGRAL, INTEGRAL2.
            
            [varargout{1:nargout}] = mpintegral2Calc(varargin{:});
        end
        
        function varargout = quad(varargin)
        %QUAD Numerically evaluate integral, adaptive Simpson quadrature
        %
        %    Usage is identical to built-in function QUAD.
        %
        %    See also QUAD
            
            [varargout{1:nargout}] = mpquad(varargin{:});
        end

        function varargout = dblquad(varargin)
        %DBLQUAD Numerically evaluate double integral over rectangle
        %
        %    Usage is identical to built-in function DBLQUAD.
        %
        %    See also DBLQUAD
            
            [varargout{1:nargout}] = mpdblquad(varargin{:});
        end

        function varargout = triplequad(varargin)
        %TRIPLEQUAD Numerically evaluate triple integral
        %
        %    Usage is identical to built-in function TRIPLEQUAD.
        %
        %    See also TRIPLEQUAD
            
            [varargout{1:nargout}] = mptriplequad(varargin{:});
        end
        
        function Q = quadgl(funfcn,a,b,n,varargin)
        %QUADGL  Numerically evaluate integral by fixed order Gauss-Legendre quadrature.   
        %   Q = QUADGL(FUN,A,B,N) tries to approximate the integral of scalar-valued
        %   function FUN from A to B by Gauss-Legendre quaradure of order N.
        %   FUN is a function handle. The function Y=FUN(X) should accept a vector 
        %   argument X and return a vector result Y, the integrand evaluated at each element of X.
        % 
        %   Intended to be used with multiple-precision numbers A, B. 
        %   FUN should be able to work with multiple-precision arguments. 
       
            narginchk(3,inf);        
            f = fcnchk(funfcn);
           
            if ~ismp(a), a = mp(a); end;
            if ~ismp(b), b = mp(b); end;
            
            if nargin < 4 || isempty(n), n = 16; end;
            if ~isscalar(a) || ~isscalar(b)
                error('MCT:quadgl:scalarLimits',...
                    'The limits of integration must be scalars.');
            end
            
            [x,w] = mp.GaussLegendre(n);
            d = (b-a)/2;
            c = (b+a)/2;
            Q = d*sum(w.*f(d.*x+c));
        end
        
        %% ODE
        function varargout = ode45(varargin)
            [varargout{1:nargout}] = mpode45(varargin{:});
        end
        
        function varargout = ode113(varargin)
            [varargout{1:nargout}] = mpode113(varargin{:});
        end
        
        function varargout = ode15s(varargin)
            [varargout{1:nargout}] = mpode15s(varargin{:});
        end
        
        function varargout = daeic12(varargin)
            [varargout{1:nargout}] = mpdaeic12(varargin{:});
        end
        
        function varargout = daeic3(varargin)
            [varargout{1:nargout}] = mpdaeic3(varargin{:});
        end
        
        function varargout = odearguments(varargin)
            [varargout{1:nargout}] = mpodearguments(varargin{:});
        end
        
        function varargout = odeevents(varargin)
            [varargout{1:nargout}] = mpodeevents(varargin{:});
        end
        
        function varargout = odemass(varargin)
            [varargout{1:nargout}] = mpodemass(varargin{:});
        end
        
        function varargout = odejacobian(varargin)
            [varargout{1:nargout}] = mpodejacobian(varargin{:});
        end

        function varargout = odenumjac(varargin)
            [varargout{1:nargout}] = mpodenumjac(varargin{:});
        end
        
        function varargout = odeplot(varargin)
            [varargout{1:nargout}] = mpodeplot(varargin{:});
        end
        
        function varargout = ntrp45(varargin)
            [varargout{1:nargout}] = mpntrp45(varargin{:});
        end
        
        function varargout = ntrp113(varargin)
            [varargout{1:nargout}] = mpntrp113(varargin{:});
        end
        
        function varargout = ntrp15s(varargin)
            [varargout{1:nargout}] = mpntrp15s(varargin{:});
        end
        
        function varargout = odefinalize(varargin)
            [varargout{1:nargout}] = mpodefinalize(varargin{:});
        end
        
        function varargout = odezero(varargin)
            [varargout{1:nargout}] = mpodezero(varargin{:});
        end
       
        %% Fast Fourier Transform
        function varargout = fft(varargin)
        %FFT Discrete Fourier transform.     
        %
        %    Usage is identical to built-in function FFT
        %
        %    See also FFT2, FFTN, FFTSHIFT, FFTW, IFFT, IFFT2, IFFTN.
        
        % In case if dim > ndims(A) and there is n, then: 
        % B = A;
        % for i=1:(n-1), B = cat(dim,B,A);end;

            [varargout{1:nargout}] = mpimpl(7000,varargin{:});
        end
        
        function varargout = fft2(varargin)
        %FFT2 Two-dimensional discrete Fourier Transform.            
        %
        %    Usage is identical to built-in function FFT2
        %
        %   See also FFT, FFTN, FFTSHIFT, FFTW, IFFT, IFFT2, IFFTN.
        
            [varargout{1:nargout}] = mpfft2(varargin{:});
        end
        
        function Y = fftn(X,varargin)
        %FFTN N-dimensional discrete Fourier Transform.           
        %
        %    Usage is identical to built-in function FFTN
        %    except that two-parameter version Y = fftn(X,siz)
        %    is not supported.
        %
        %   See also FFT, FFT2, FFTSHIFT, FFTW, IFFT, IFFT2, IFFTN.

            if nargin > 1
               error('%s\n%s','Two-parameter call to FFTN is not supported in multiprecision mode.',...
                              'Please request it using support@advanpix.com');
            end
                
            Y = X;
            for p = 1:length(size(X))
                Y = fft(Y,[],p);
            end        
        end
        
        function varargout = ifft(varargin)
        %IFFT Inverse discrete Fourier transform.            
        %
        %    Usage is identical to built-in function FFT2
        %
        %   See also FFT, FFT2, FFTN, FFTSHIFT, FFTW, IFFT2, IFFTN.
        
            [varargout{1:nargout}] = mpimpl(7001,varargin{:});
        end
        
        function varargout = ifft2(varargin)
        %IFFT2 Two-dimensional inverse discrete Fourier transform.            
        %
        %    Usage is identical to built-in function FFT2
        %
        %   See also FFT, FFT2, FFTN, FFTSHIFT, FFTW, IFFT, IFFTN.
        
            [varargout{1:nargout}] = mpifft2(varargin{:});
        end
        
        function Y = ifftn(X,varargin)
        %IFFTN N-dimensional inverse discrete Fourier transform.           
        %
        %    Usage is identical to built-in function IFFTN
        %    except that we support only one-parameter call Y = ifftn(X)
        %
        %   See also FFT, FFT2, FFTSHIFT, FFTW, IFFT, IFFT2, FFTN.

            if nargin > 1
               error('%s\n%s','Two-parameter call to IFFTN is not supported in multiprecision mode.',...
                              'Please request it using support@advanpix.com');
            end
                
            Y = X;
            for p = 1:length(size(X))
                Y = ifft(Y,[],p);
            end        
        end
        
        
        %% Polynomial roots
        function varargout = roots(varargin)
        %ROOTS Polynomial roots.
        %
        %    Usage is identical to built-in function ROOTS
        %    
        %
        %    See also ROOTS.
            
            [varargout{1:nargout}] = mproots(varargin{:});
        end

        function varargout = poly(varargin)
        %PLOY Polynomial with specified roots.
        %
        %    Usage is identical to built-in function POLY
        %    
        %
        %    See also POLY, ROOTS, CONV.
            
            [varargout{1:nargout}] = mppoly(varargin{:});
        end
        
        function varargout = polyder(varargin)
        %POLYDER Differentiate polynomial.
        %
        %    Usage is identical to built-in function POLYDER
        %    
        %
        %    See also POLYINT, CONV, DECONV.
            
            [varargout{1:nargout}] = mppolyder(varargin{:});
        end
        
        function varargout = polyeig(varargin)
        %POLYEIG Polynomial eigenvalue problem.
        %
        %    Implemented in quadruple precision only.
        %
        %    Usage is identical to built-in function POLYEIG
        %    
        %
        %    See also EIG, COND, CONDEIG.
            
            [varargout{1:nargout}] = mppolyeig(varargin{:});
        end
        
        function varargout = polyval(varargin)
        %POLYVAL Evaluate polynomial.
        %
        %
        %    Usage is identical to built-in function POLYVAL
        %    
        %
        %    See also POLYFIT, POLYVALM.
            
            [varargout{1:nargout}] = mppolyval(varargin{:});
        end
        
        function varargout = polyvalm(varargin)
        %POLYVALM Evaluate polynomial with matrix argument.
        %
        %
        %    Usage is identical to built-in function POLYVALM
        %    
        %
        %    See also POLYVAL, POLYFIT.
            
            [varargout{1:nargout}] = mppolyvalm(varargin{:});
        end

        function varargout = polyfit(varargin)
        %POLYFIT Evaluate polynomial.
        %
        %
        %    Usage is identical to built-in function POLYFIT
        %    
        %
        %    See also POLY, POLYVAL, ROOTS.
            
            [varargout{1:nargout}] = mppolyfit(varargin{:});
        end

        function varargout = polyint(varargin)
        %POLYINT Integrate polynomial analytically.
        %
        %
        %    Usage is identical to built-in function POLYINT
        %    
        %
        %   See also RESIDUEPOLYDER, POLYVAL, POLYVALM, POLYFIT.
            
            [varargout{1:nargout}] = mppolyint(varargin{:});
        end
        
        function varargout = residue(varargin)
        %RESIDUE Partial-fraction expansion (residues).
        %
        %
        %    Usage is identical to built-in function RESIDUEZ
        %    
        %
        %    See also POLY, ROOTS, DECONV.
            
            [varargout{1:nargout}] = mpresidue(varargin{:});
        end

        function varargout = residuez(varargin)
        %RESIDUEZ Z-transform partial-fraction expansion.            
        %
        %
        %    Usage is identical to built-in function RESIDUE
        %    
        %
        %    See also RESIDUE, POLY, ROOTS, DECONV.
            
            [varargout{1:nargout}] = mpresiduez(varargin{:});
        end
        
        function varargout = mpoles(varargin)
        %MPOLES Identify repeated poles & their multiplicities.
        %
        %
        %    Usage is identical to built-in function MPOLES
        %    
        %
        %   See also RESIDUE.
            
            [varargout{1:nargout}] = mpmpoles(varargin{:});
        end
        
        function varargout = resi2(varargin)
        %RESI2  Residue of a repeated pole.
        %
        %
        %    Usage is identical to built-in function RESI2
        %    
        %
        %   See also RESIDUE.
            
            [varargout{1:nargout}] = mpresi2(varargin{:});
        end

        function varargout = conv(varargin)
        %CONV Convolution and polynomial multiplication.
        %
        %    Usage is identical to built-in function CONV
        %
        %    See also CONV.
            
            [varargout{1:nargout}] = mpimpl(7002, varargin{:});
        end

        function varargout = filter(varargin)
        %FILTER One-dimensional digital filter.
        %
        %    Usage is identical to built-in function FILTER
        %
        %    See also CONV.
            
            [varargout{1:nargout}] = mpimpl(7003, varargin{:});
        end

        function varargout = deconv(varargin)
        %DECONV Deconvolution and polynomial division.
        %
        %    Usage is identical to built-in function DECONV
        %
        %    See also DECONV.
            
            [varargout{1:nargout}] = mpdeconv(varargin{:});
        end
        
        %% Array Creation Functions
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
        
        function r = ones(varargin)
        %ONES Create array of all ones.
        %
        %    This overload is called for all cases when one of the parameters are of 'mp' type: 
        %    (M N, x or else):
        %
        %    ones([M N ... ],'mp')
        %    ones([M N ... ],'like',x), when isa(x,'mp') = true
        %    ones([M N ... ])
        %
        %    Usage is identical to built-in function ONES.
        %    
        %    See also ZEROS, EYE.

            r = mp.genericArrayCreationUsingBuiltIn('ones',varargin{:});
        end

        function r = zeros(varargin)
        %ZEROS Create array of all zeros.
        %
        %    This overload is called for all cases when one of the parameters are of 'mp' type: 
        %    (M N, x or else):
        %
        %    zeros([M N ... ],'mp')
        %    zeros([M N ... ],'like',x), when isa(x,'mp') = true
        %    zeros([M N ... ])
        %
        %    Usage is identical to built-in function ZEROS.
        %    
        %    See also ONES, EYE.
         
            r = mp.genericArrayCreationUsingBuiltIn('zeros',varargin{:});
        end

        function r = eye(varargin)
        %EYE Identity matrix.
        %
        %    This overload is called for all cases when one of the parameters are of 'mp' type: 
        %    (M N, x or else):
        %
        %    eye([M N ... ],'mp')
        %    eye([M N ... ],'like',x), when isa(x,'mp') = true
        %    eye([M N ... ])
        %
        %    Usage is identical to built-in function EYE.
        %    
        %    See also ZEROS, ONES.
        
            r = mp.genericArrayCreationUsingBuiltIn('eye',varargin{:});
        end

        function r = Inf(varargin)
        %INF Infinity.
        %
        %    This overload is called for all cases when one of the parameters are of 'mp' type: 
        %    (M N, x or else):
        %
        %    Inf([M N ... ],'mp')
        %    Inf([M N ... ],'like',x), when isa(x,'mp') = true
        %    Inf([M N ... ])
        %
        %    Usage is identical to built-in function INF.
        %    
        %    See also ZEROS, ONES, EYE, NAN.
        
            r = mp.genericArrayCreationUsingBuiltIn('Inf',varargin{:});
        end
        
        function r = NaN(varargin)
        %NaN Not-a-Number.
        %
        %    This overload is called for all cases when one of the parameters are of 'mp' type: 
        %    (M N, x or else):
        %
        %    NaN([M N ... ],'mp')
        %    NaN([M N ... ],'like',x), when isa(x,'mp') = true
        %    NaN([M N ... ])
        %
        %    Usage is identical to built-in function NAN.
        %    
        %    See also ZEROS, ONES, EYE, INF.
        
            r = mp.genericArrayCreationUsingBuiltIn('NaN',varargin{:});
        end
        
        function r = rand(varargin)
        %RAND Uniformly distributed pseudorandom numbers.
        %
        %    Use mprand(...) for full precision random numbers.
        %
        %    Usage is identical to built-in function RAND.
        %    
        %    See also RANDN.
            
            r = mp.genericArrayCreationUsingBuiltIn('rand',varargin{:});            
        end

        function r = randn(varargin)
        %RANDN Normally distributed pseudorandom numbers.
        %
        %    Use mprandn(...) for full precision random numbers.
        %
        %    Usage is identical to built-in function RANDN.
        %    
        %    See also RAND.
            
            r = mp.genericArrayCreationUsingBuiltIn('randn',varargin{:});                        
        end
        
        function r = randi(varargin)
        %RANDI Pseudorandom integers from a uniform discrete distribution.
        %
        %    Usage is identical to built-in function RANDI.
        %    
        %    See also RAND, RANDN.
            
            r = mp.genericArrayCreationUsingBuiltIn('randi',varargin{:});                        
        end

        function varargout = linspace(varargin)
        %LINSPACE Linearly spaced vector.
        %
        %    Usage is identical to built-in function LINSPACE.
        %
        %    See also LOGSPACE, COLON.
        
               [varargout{1:nargout}] = mplinspace(varargin{:});
        end

        function varargout = logspace(varargin)
        %LOGSPACE Logarithmically spaced vector.
        %
        %    Usage is identical to built-in function LOGSPACE.
        %
        %    See also LINSPACE, COLON.
        
               [varargout{1:nargout}] = mplogspace(varargin{:});
        end
        
        function varargout = meshgrid(varargin)
        %MESHGRID Cartesian grid in 2-D/3-D space
        %
        %    Usage is identical to built-in function MESHGRID.
        %
        %    See also SURF, SLICE, NDGRID.
        
               [varargout{1:nargout}] = mpmeshgrid(varargin{:});
        end

        function varargout = ndgrid(varargin)
        %NDGRID Rectangular grid in N-D space
        %
        %    Usage is identical to built-in function NDGRID.
        %
        %    See also MESHGRID, SLICE, INTERPN.
        
               [varargout{1:nargout}] = mpndgrid(varargin{:});
        end
        
        %% Interpolation
        function varargout = interp1(varargin)
        %INTERP1 1-D interpolation (table lookup)
        %
        %   Usage is identical to built-in function INTERP1.
        %    
        %   See also INTERPFT, SPLINE, PCHIP, INTERP2, INTERP3, INTERPN, PPVAL.
            for i=1:nargin, if (~ischar(varargin{i}) && ~isa(varargin{i}, 'mp')), varargin{i} = mp(varargin{i}); end; end
            [varargout{1:nargout}] = mpinterp1(varargin{:});
        end
        
        function varargout = pchip(varargin)
        %PCHIP Piecewise Cubic Hermite Interpolating Polynomial.
        %
        %   Usage is identical to built-in function PCHIP.
        %    
        %   See also INTERP1, SPLINE, PPVAL, UNMKPP.
            for i=1:nargin, if ~isa(varargin{i}, 'mp'), varargin{i} = mp(varargin{i}); end; end
            [varargout{1:nargout}] = mppchip(varargin{:});
        end
        
        function varargout = pwch(varargin)
        %PWCH Piecewise cubic Hermite interpolation
        %
        %   Usage is identical to built-in function PWCH.
        %    
        %   See also SPLINE, INTERP1, PCHIP, PPVAL, UNMKPP.
            for i=1:nargin, if ~isa(varargin{i}, 'mp'), varargin{i} = mp(varargin{i}); end; end
            [varargout{1:nargout}] = mppwch(varargin{:});
        end
        
        function varargout = spline(varargin)
        %SPLINE Cubic spline data interpolation.
        %
        %   Usage is identical to built-in function SPLINE.
        %    
        %   See also INTERP1, PCHIP, PPVAL, MKPP, UNMKPP.
            for i=1:nargin, if ~isa(varargin{i}, 'mp'), varargin{i} = mp(varargin{i}); end; end
            [varargout{1:nargout}] = mpspline(varargin{:});
        end
        
        function varargout = chckxy(varargin)
        %CHCKXY Check and adjust input for SPLINE and PCHIP
        %
        %   Usage is identical to built-in function CHCKXY.
        %    
        %   See also PCHIP, SPLINE.
            for i=1:nargin, if ~isa(varargin{i}, 'mp'), varargin{i} = mp(varargin{i}); end; end
            [varargout{1:nargout}] = mpchckxy(varargin{:});
        end
        
        function varargout = mkpp(varargin)
        %MKPP Make piecewise polynomial.
        %
        %   Usage is identical to built-in function MKPP.
        %    
        %   See also UNMKPP, PPVAL, SPLINE.
            [varargout{1:nargout}] = mpmkpp(varargin{:});
        end
        
        function varargout = unmkpp(varargin)
        %UNMKPP Supply details about piecewise polynomial.
        %
        %   Usage is identical to built-in function UNMKPP.
        %    
        %   See also MKPP, PPVAL, SPLINE.
            [varargout{1:nargout}] = mpunmkpp(varargin{:});
        end
        
        function varargout = ppval(varargin)
        %PPVAL Evaluate piecewise polynomial.
        %
        %   Usage is identical to built-in function PPVAL.
        %    
        %   See also SPLINE, PCHIP, INTERP1, MKPP, UNMKPP.
            [varargout{1:nargout}] = mpppval(varargin{:});
        end
        
        %% Iterative solvers
        function varargout = bicg(varargin)
        %BICG BiConjugate Gradients Method.
        %
        %   Usage is identical to built-in function BICG.
        %    
        %   See also BICGSTAB, BICGSTABL, CGS, LSQR, MINRES, PCG, QMR, SYMMLQ,
        %   TFQMR, ILU, FUNCTION_HANDLE.
        
            [varargout{1:nargout}] = mpbicg(varargin{:});
        end
        
        function varargout = bicgstab(varargin)
        %BICGSTAB BiConjugate Gradients Stabilized Method.
        %
        %   Usage is identical to built-in function BICGSTAB.
        %    
        %   See also BICG, BICGSTABL, CGS, LSQR, MINRES, PCG, QMR, SYMMLQ,
        %   TFQMR, ILU, FUNCTION_HANDLE.
        
            [varargout{1:nargout}] = mpbicgstab(varargin{:});
        end
        
        function varargout = bicgstabl(varargin)
        %BICGSTABL BiConjugate Gradients Stabilized(l) Method.
        %
        %   Usage is identical to built-in function BICGSTABL.
        %    
        %   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ,
        %   TFQMR, ILU, FUNCTION_HANDLE.
        
            [varargout{1:nargout}] = mpbicgstabl(varargin{:});
        end

        function varargout = cgs(varargin)
        %CGS Conjugate Gradients Squared Method.
        %
        %   Usage is identical to built-in function CGS.
        %    
        %   See also BICG, BICGSTAB, BICGSTABL, LSQR, MINRES, PCG, QMR, SYMMLQ,
        %   TFQMR, ILU, FUNCTION_HANDLE.
        
            [varargout{1:nargout}] = mpcgs(varargin{:});
        end
        
        function varargout = gmres(varargin)
        %GMRES Generalized Minimum Residual Method.
        %
        %   Usage is identical to built-in function GMRES.
        %    
        %   See also BICG, BICGSTAB, BICGSTABL, CGS, LSQR, MINRES, PCG, QMR, SYMMLQ,
        %   TFQMR, ILU, FUNCTION_HANDLE.
        
            [varargout{1:nargout}] = mpgmres(varargin{:});
        end

        function varargout = pcg(varargin)
        %PCG Preconditioned Conjugate Gradients Method.
        %
        %   Usage is identical to built-in function PCG.
        %    
        %   See also BICG, BICGSTAB, BICGSTABL, LSQR, CGS, MINRES, QMR, SYMMLQ,
        %   TFQMR, ILU, FUNCTION_HANDLE.
        
            [varargout{1:nargout}] = mpcgs(varargin{:});
        end
        
        function varargout = minres(varargin)
        %MINRES Minimum Residual Method.
        %
        %   Usage is identical to built-in function MINRES.
        %    
        %   See also BICG, BICGSTAB, BICGSTABL, CGS, GMRES, LSQR, PCG, QMR, SYMMLQ,
        %   TFQMR, ICHOL, ILU, FUNCTION_HANDLE.
        
            [varargout{1:nargout}] = mpminres(varargin{:});
        end
        
        function varargout = iterchk(varargin)
        %ITERCHK Checks arguments to iterative methods.
        %
        %   Usage is identical to built-in function ITERCHK.
        %    
        %   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ.
        
            [varargout{1:nargout}] = mpiterchk(varargin{:});
        end
        
        function varargout = iterapp(varargin)
        %ITERAPP Apply matrix operator to vector and error gracefully.
        %
        %   Usage is identical to built-in function ITERAPP.
        %    
        %   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ.
        
            [varargout{1:nargout}] = mpiterapp(varargin{:});
        end

        function varargout = itermsg(varargin)
        %ITERMSG Apply matrix operator to vector and error gracefully.
        %
        %   Usage is identical to built-in function ITERMSG.
        %    
        %   See also BICG, BICGSTAB, BICGSTABL, CGS, GMRES, LSQR, MINRES, PCG, QMR,
        %   SYMMLQ, TFQMR.
        
            [varargout{1:nargout}] = mpitermsg(varargin{:});
        end

        %% Sylvester, Lyapunov and Riccati solvers && selected control routines
        function varargout = sylvester_tri(varargin)
        %SYLVESTER_TRI Solve (Quasi-)Triangular Sylvester Equations.
        %   X = SYLVESTER_TRI(A,B,C) solves the (quasi-) triangular Sylvester equation A*X + X*B = C,
        %   where A is a m-by-m matrix, B is a n-by-n matrix, and X and C are
        %   m-by-n matrices. 
        %   
        %   If A, B and C are real, A and B must be in Schur canonical form (upper quasi-triangular). 
        %   Otherwise, A and B must be upper triangular. SYLVESTER_TRI assumes
        %   these properties and does not perform check.
        %
        %    See also SYLVESTER
        
            [varargout{1:nargout}] = mpimpl(4008,varargin{:});
        end
        
        function varargout = sylvester(varargin)
        %SYLVESTER Solve Sylvester Equations.
        %   X = SYLVESTER(A,B,C) solves the Sylvester equation A*X + X*B = C,
        %   where A is a m-by-m matrix, B is a n-by-n matrix, and X and C are
        %   m-by-n matrices.
        %
        %   The algorithm is based on Schur decomposition due to Bartels and Stewart [1] 
        %   and costs O(m^3+n^3) operations [2]. See [3] for its numerical
        %   stability.
        % 
        %   References:
        %
        %   [1] R. H. Bartels & G. W. Stewart, Solution of the matrix equation 
        %       AX+XB = C, Comm. ACM, 15, 820–826, 1972.
        %   [2] J. D. Gardiner, A. J. Laub, J. J. Amato, & C. B. Moler, Solution of the
        %       Sylvester matrix equation AXB^T + CXD^T = E, ACM Transactions on
        %       Mathematical Software (TOMS), 18(2), 223-231, 1992.
        %   [3] N. J. Higham, Accuracy and Stability of Numerical
        %       Algorithms, Chapter 16, 2002, 2nd edition, SIAM.
        %
        %    See also SYLVESTER_TRI
        
            [varargout{1:nargout}] = mpsylvester(varargin{:});
        end

        function varargout = lyap(varargin)
        %LYAP  Solve continuous-time Lyapunov equations.
        %
        %   Usage is identical to function LYAP from control toolbox.
        %
        %   X = LYAP(A,Q) solves the Lyapunov matrix equation:
        %
        %       A*X + X*A' + Q = 0
        %
        %   X = LYAP(A,B,C) solves the Sylvester equation:
        %
        %       A*X + X*B + C = 0
        %
        %   X = LYAP(A,Q,[],E) solves the generalized Lyapunov equation:
        %
        %       A*X*E' + E*X*A' + Q = 0    where Q is symmetric
        %
        %   References:
        %
        %     [1] Barraud, A.Y.                   
        %         A numerical algorithm to solve A'XA - X = Q.
        %         IEEE Trans. Auto. Contr., AC-22, pp. 883-885, 1977.
        %
        %     [2] Bartels, R.H. and Stewart, G.W.  
        %         Solution of the matrix equation A'X + XB = C.
        %         Comm. A.C.M., 15, pp. 820-826, 1972.
        %
        %     [3] Hammarling, S.J.
        %         Numerical solution of the stable, non-negative definite
        %         Lyapunov equation.
        %         IMA J. Num. Anal., 2, pp. 303-325, 1982.
        %
        %     [4] Penzl, T. Numerical solution of generalized Lyapunov equations.
        %         Advances in Comp. Math., vol. 8, pp. 33-48, 1998.
        %
        %   See also LYAPCHOL, DLYAP.
        
            [varargout{1:nargout}] = mplyap(varargin{:});
        end

        function varargout = dlyap(varargin)
        %DLYAP  Solve discrete Lyapunov equations.
        %
        %   Usage is identical to function LYAP from control toolbox.
        %
        %   X = DLYAP(A,Q) solves the discrete Lyapunov matrix equation:
        %
        %       A*X*A' - X + Q = 0
        %
        %   X = DLYAP(A,B,C) solves the Sylvester equation:
        %
        %       A*X*B - X + C = 0
        %
        %   X = DLYAP(A,Q,[],E) solves the generalized discrete Lyapunov equation:
        %
        %       A*X*A' - E*X*E' + Q = 0
        %
        %   References:
        %
        %     [1] Barraud, A.Y.                   
        %         A numerical algorithm to solve A'XA - X = Q.
        %         IEEE Trans. Auto. Contr., AC-22, pp. 883-885, 1977.
        %
        %     [2] Bartels, R.H. and Stewart, G.W.  
        %         Solution of the matrix equation A'X + XB = C.
        %         Comm. A.C.M., 15, pp. 820-826, 1972.
        %
        %     [3] Hammarling, S.J.
        %         Numerical solution of the stable, non-negative definite
        %         Lyapunov equation.
        %         IMA J. Num. Anal., 2, pp. 303-325, 1982.
        %
        %     [4] Penzl, T. Numerical solution of generalized Lyapunov equations.
        %         Advances in Comp. Math., vol. 8, pp. 33-48, 1998.
        %
        %   See also DLYAPCHOL, LYAP.
        
            [varargout{1:nargout}] = mpdlyap(varargin{:});
        end
        
        function varargout = lyapcheckin(varargin)
        %LYAPCHECKIN Validates input arguments to LYAP and DLYAP.
        %
        %   See also LYAP, DLYAP.
        
            [varargout{1:nargout}] = mplyapcheckin(varargin{:});
        end
        
        function varargout = aebalance(varargin)
        %AEBALANCE  Balances pair of matrices
        %
        %   See also LYAP, DLYAP.
        
            [varargout{1:nargout}] = mpaebalance(varargin{:});
        end
        
        function varargout = care(varargin)
        %CARE Solve continuous-time algebraic Riccati equations.
            [varargout{1:nargout}] = mpcare(varargin{:});
        end
        
        function varargout = dare(varargin)
        %DARE Solve discrete-time algebraic Riccati equations.
            [varargout{1:nargout}] = mpdare(varargin{:});
        end
        
        function varargout = gcare(varargin)
        %GCARE Generalized solver for continuous algebraic Riccati equations.
            [varargout{1:nargout}] = mpgcare(varargin{:});
        end
        
        function varargout = gdare(varargin)
        %GDARE Generalized solver for discrete algebraic Riccati equations.
            [varargout{1:nargout}] = mpgdare(varargin{:});
        end
        
        function varargout = arecheckin(varargin)
        %ARECHECKIN Checks input arguments to CARE and DARE.
            [varargout{1:nargout}] = mparecheckin(varargin{:});
        end
        
        function varargout = arecheckout(varargin)
        %ARECHECKOUT Checks for proper extraction of stable invariant subspace.
            [varargout{1:nargout}] = mparecheckout(varargin{:});
        end
        
        function varargout = arefact2x(varargin)
        %AREFACT2X Computes Riccati solution X = D*(X2/X1)*D from X1,X2,D.
            [varargout{1:nargout}] = mparefact2x(varargin{:});
        end
        
        function varargout = lrscale(varargin)
        %LRSCALE  Applies left and right scaling matrices.
            [varargout{1:nargout}] = mplrscale(varargin{:});
        end
        
        function varargout = arescale(varargin)
        %ARESCALE  Computes diagonal scaling for GARE.
            [varargout{1:nargout}] = mparescale(varargin{:});
        end
        
        function varargout = mscale(varargin)
        %MSCALE  Row/column scaling of 2D matrix.
            [varargout{1:nargout}] = mpmscale(varargin{:});
        end
        
        %% Overloads for Graphics
        function varargout = axis(varargin)
        %AXIS  Control axis scaling and appearance.
            for i=1:nargin, if isa(varargin{i},'mp'), varargin{i} = double(varargin{i}); end; end;        
            [varargout{1:nargout}] = axis(varargin{:});
        end

        function varargout = histogram(varargin)
        %HISTOGRAM  Plots a histogram.
            for i=1:nargin, if isa(varargin{i},'mp'), varargin{i} = double(varargin{i}); end; end;        
            [varargout{1:nargout}] = histogram(varargin{:});
        end
        
        function varargout = text(varargin)
        %TEXT   Text annotation.
            for i=1:nargin, if isa(varargin{i},'mp'), varargin{i} = double(varargin{i}); end; end;        
            [varargout{1:nargout}] = text(varargin{:});
        end
        
        %% Optimization
        function varargout = fzero(varargin)
        %FZERO  Single-variable nonlinear zero finding. 
        %
        %   Usage is identical to built-in function FZERO.
        %    
        %   See also ROOTS, FMINBND, FUNCTION_HANDLE.
        
            [varargout{1:nargout}] = mpfzero(varargin{:});
        end

        function varargout = fminsearch(varargin)
        %FMINSEARCH Multidimensional unconstrained nonlinear minimization (Nelder-Mead).
        %
        %   Usage is identical to built-in function FMINSEARCH.
        %    
        %   See also OPTIMSET, FMINBND, FUNCTION_HANDLE.
        
            [varargout{1:nargout}] = mpfminsearch(varargin{:});
        end
        
        function varargout = isoptimargdbl(varargin)
        %ISOPTIMARGDBL returns true when input arguments are double or single.
        %
        %   Usage is identical to built-in function ISOPTIMARGDBL.
        %    
        %   See also OPTIMSET, FMINBND, FUNCTION_HANDLE.
        
            [varargout{1:nargout}] = mpisoptimargdbl(varargin{:});
        end
        
        function varargout = fsolve(varargin)
        %FSOLVE solves systems of nonlinear equations of several variables.
        %
        %   Usage is identical to built-in function FSOLVE.
        %    
        %   See also OPTIMOPTIONS, LSQNONLIN.
        
            [varargout{1:nargout}] = mpfsolve(varargin{:});
        end
        
        function varargout = lsqnonneg(varargin)
        %LSQNONNEG Linear least squares with nonnegativity constraints.
        %
        %   X = LSQNONNEG(C,d) returns the vector X that minimizes NORM(d-C*X)
        %   subject to X >= 0. C and d must be real.
        %
        %   Usage is identical to built-in function lsqnonneg.
        %    
        %   See also OPTIMOPTIONS, LSQNONLIN.
        
            [varargout{1:nargout}] = mplsqnonneg(varargin{:});
        end

        %% Internal/shared routines needed for high-level optimization solvers 
        function varargout = snls(varargin)
        %SNLS  Sparse nonlinear least squares solver.
        %
        %   Usage is identical to built-in function SLNS.
        %    
        %   See also FSOLVE.
        
            [varargout{1:nargout}] = mpsnls(varargin{:});
        end
        
        function varargout = sfdnls(varargin)
        %SFDNLS  SFDNLS Estimate Jacobian via finite differences
        %
        %   Usage is identical to built-in function SFDNLS.
        %    
        %   See also FSOLVE.
        
            [varargout{1:nargout}] = mpsfdnls(varargin{:});
        end

        function varargout = definev(varargin)
        %DEFINEV Scaling vector and derivative
        %
        %   Usage is identical to built-in function DEFINEV.
        %    
        %   See also FSOLVE.
        
            [varargout{1:nargout}] = mpdefinev(varargin{:});
        end
        
        function varargout = trdog(varargin)
        %TRDOG Reflected (2-D) trust region trial step (box constraints)
        %
        %   Usage is identical to built-in function TRDOG.
        %    
        %   See also FSOLVE.
        
            [varargout{1:nargout}] = mptrdog(varargin{:});
        end

        function varargout = pcgr(varargin)
        %PCGR	Preconditioned conjugate gradients
        %
        %   Usage is identical to built-in function PCGR.
        %    
        %   See also FSOLVE.
        
            [varargout{1:nargout}] = mppcgr(varargin{:});
        end

        function varargout = trust(varargin)
        %TRUST	Exact soln of trust region problem
        %
        %   Usage is identical to built-in function TRUST.
        %    
        %   See also FSOLVE.
        
            [varargout{1:nargout}] = mptrust(varargin{:});
        end
        
        function varargout = aprecon(varargin)
        %APRECON Banded preconditioner function for least-squares problems.
        %
        %   Usage is identical to built-in function TRUST.
        %    
        %   See also FSOLVE.
        
            [varargout{1:nargout}] = mpaprecon(varargin{:});
        end

        function varargout = perturbTrustRegionReflective(varargin)
        %perturbTrustRegionReflective Perturb point from bounds
        %
        %   Usage is identical to built-in function TRUST.
        %    
        %   See also FSOLVE.
        
            [varargout{1:nargout}] = mpperturbTrustRegionReflective(varargin{:});
        end
        
        function varargout = trustnleqn(varargin)
        %TRUSTNLEQN Trust-region dogleg nonlinear systems of equation solver.
        %
        %   Usage is identical to built-in function TRUSTNLEQN.
        %    
        %   See also FSOLVE, SNLS.
        
            [varargout{1:nargout}] = mptrustnleqn(varargin{:});
        end
        
        function varargout = levenbergMarquardt(varargin)
        %levenbergMarquardt Levenberg-Marquardt solver of non-linear least squares problems.
        %
        %   Usage is identical to built-in function levenbergMarquardt.
        %    
        %   See also FSOLVE, SNLS.
        
            [varargout{1:nargout}] = mplevenbergMarquardt(varargin{:});
        end
        
        function y = nonzerosign(x)
        %NONZEROSIGN Signum function excluding zero
        %   For each element of X, NONZEROSIGN(X) returns 1 if the element
        %   is greater than or equal to zero and -1 if it is
        %   less than zero. NONZEROSIGN differs from SIGN in that NONZEROSIGN(0)
        %   returns 1, while SIGN(0) returns 0.
        
            y = ones(size(x),'mp');
            y(x < 0) = mp(-1);        
        end
        %% Formatted output
        function s = num2str(x, f)
        %NUM2STR Convert number to string.
        %
        %    Usage is identical to built-in function NUM2STR
        %
        %    See also SPRINTF.
        
        %
        %    TODO:
        %    1. Better handling of %i and %d specifiers when x is a matrix.
        %       See the sprintf_basic for more details. 
        %       Currently matrix conversion is implemented on C++ level with spacing, etc.
        %       Not sure what way is better: (a) to re-implement it in Matlab
        %       or (b) extend C++ version with %i, %d support.
        %
        %       Once this will be implemented - then we can simplify
        %       sprintf_basic (remove unsupported_types handling), since
        %       all the low-level stuff will happen here.
        %
        %    2. Add currently unsupported format specifiers - hex, etc.        
        %
        %    3. Complex values should be converted as in Matlab (imag part
        %       under real).
        %    4. Alignment and spaces. See source code for built-in num2str
        %       for more details.
        
            if isempty(x) ,  s = '';      return;  end
            if issparse(x),  x = full(x);          end
            
            if     nargin < 2  , s = mpimpl(6000,x);   return;
            elseif isnumeric(f), s = mpimpl(6000,x,f); return;
            elseif ischar(f)
                
                if isempty(strfind(f,'%')), error(message('MATLAB:num2str:fmtInvalid', f)); end; %#ok<STREMP>
                
                [match, tokenname, split] = mp.parseSpecifiers(f);
                argc = size(match, 2);
                
                if argc==0, error(message('MATLAB:num2str:fmtInvalid', f)); end;
                
                if argc==1
                       % Only one format specifier 
                       % x might be a matrix or scalar.
                       % all recursive calls from matrix-enabled routines,
                       % sprintf, num2str end here.
                       
                       unsupported_types = '([diuoxXcs]|(?:bx|bX|bo|bu|tx|tX|to|tu|ld|li|lo|lu|lx|lX|hd|hi|ho|hu|hx|hX))';
                       if regexp(tokenname(1).type, unsupported_types)
                           % Just use %g for all unsupported types for now
                           % Later we will need to add proper handling for
                           % integer types.
                           modifiedFormat = regexprep(match{1},unsupported_types,'g');                       
                           s = mpimpl(6000,x,[split{1} modifiedFormat split{2}]);
                       else
                           s = mpimpl(6000,x,f);
                       end
                       return;
                end;
                
                [rows,cols] = size(x);
                if cols > argc
                    fspec = '';
                    for k=1:floor(cols/argc), fspec = [fspec f]; end; %#ok<AGROW>
                    for k=1:mod  (cols,argc), fspec = [fspec split{k} match{k}]; end; %#ok<AGROW>
                else
                    fspec = f;
                end
                
                s = '';
                argv = cell(1,cols);                
                for k=1:rows 
                    for p=1:cols, argv{p} = subsref(x,substruct('()',{k,p})); end;
                    
                    if k<rows,  s = [s mp.sprintf_basic([fspec '\n'],argv{:})];  %#ok<AGROW>
                    else        s = [s mp.sprintf_basic(fspec,       argv{:})]; end; %#ok<SEPEX,AGROW>
                end
            else
                error(message('MATLAB:num2str:invalidSecondArgument'))        
            end
        end
        
        function varargout = sprintf(varargin)
        % SPRINTF Write formatted data to string.      
        %
        %    Usage is identical to built-in function SPRINTF
        %
        %    See also NUM2STR.
               
            % Iterate over every argc input elements 
            % and call sprintf_basic for each assembly.
            [matchstring, tokenname, splitstring] = mp.parseSpecifiers(varargin{1}); %#ok<ASGLU>
            argc = size(matchstring, 2);
            
            if argc==0 
                % let the built-in to handle the singular cases                
                [varargout{1:nargout}] = builtin('sprintf',varargin{1});                
                return;                
            end;
                
            argv = cell(1,argc);
            
            s = '';
            n =  1;            
            for k=2:nargin
                object = varargin{k};
                if isnumeric(object)
                    object = real(object); 
                    for p=1:numel(object)
                        if ismp(object),  argv{n} = subsref(object,substruct('()',{p}));
                        else              argv{n} = object(p); end; %#ok<SEPEX>

                        if n==argc
                            s = [s mp.sprintf_basic(varargin{1},argv{:})]; %#ok<AGROW>
                            n = 1;
                        else
                            n = n+1;                        
                        end
                    end
                else
                    argv{n} = object;
                    
                    if n==argc
                        s = [s mp.sprintf_basic(varargin{1},argv{:})]; %#ok<AGROW>
                        n = 1;
                    else
                        n = n+1;                        
                    end
                end
            end
            
            if n > 1
                partial = '';
                for k=1:n-1, partial = [partial splitstring{k} matchstring{k}]; end; %#ok<AGROW>
                s = [s mp.sprintf_basic(partial,argv{1:n-1})];
            end;
            
            varargout{1} = s;
            for k=2:nargout, varargout{k} = ''; end; %#ok<AGROW>
        end
        
        function fprintf(varargin)
        %FPRINTF Write formatted data to text file.        
        %
        %    Usage is identical to built-in function FPRINTF
        %
        %    See also SPRINTF, NUM2STR.
            if isnumeric(varargin{1}), fprintf(varargin{1}, sprintf(varargin{2:end}));                
            else                       fprintf(sprintf(varargin{1:end})); end %#ok<SEPEX>
        end
        
    end % methods
   
end % classdef
