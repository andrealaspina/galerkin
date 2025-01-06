function varargout = mpeigs(varargin)
%EIGS  Find a few eigenvalues and eigenvectors of a matrix using 
%      Krylov-Schur decomposition implemented using arbitrary
%      precision precision.
%   
%   Routine is capable of working with arbitrary precision floating-point 
%   numerics provided by Advanpix Multiprecision Computing Toolbox.
%
%   Working modes and options are equivalent to classic EIGS. The only
%   difference is that we are using Krylov-Schur decomposition instead of Arnoldi (ARPACK).
%   Krylov-Schur is reportedly more numerically stable, and combined with
%   extended precision allows solving sensitive or ill-conditioned
%   eigenproblems. See references below.
%
%   Algorithm is implemented using MATLAB language and thus it is slow.
%   After short period tests and feedback from users we will move it to C++.
%   Estimated speed-up will be ~25-100 times. Stay tuned.
%
%   D = EIGS(A) returns a vector of A's 6 largest magnitude eigenvalues.
%   A must be square and should be large and sparse.
%
%   [V,D] = EIGS(A) returns a diagonal matrix D of A's 6 largest magnitude
%   eigenvalues and a matrix V whose columns are the corresponding
%   eigenvectors.
%
%   [V,D,FLAG,ITER] = EIGS(A) also returns a convergence flag and number of iterations. 
%   If FLAG is 0 then all the eigenvalues converged; otherwise not all converged.
%
%   EIGS(A,B) solves the generalized eigenvalue problem A*V == B*V*D. B must be
%   the same size as A. EIGS(A,[],...) indicates the standard eigenvalue problem
%   A*V == V*D.
%
%   EIGS(A,K) and EIGS(A,B,K) return the K largest magnitude eigenvalues.
%
%   EIGS(A,K,SIGMA) and EIGS(A,B,K,SIGMA) return K eigenvalues. If SIGMA is:
%      'LM' or 'SM' - Largest or Smallest Magnitude
%   For real symmetric problems, SIGMA may also be:
%      'LA' or 'SA' - Largest or Smallest Algebraic
%      'BE' - Both Ends, one more from high end if K is odd
%   For nonsymmetric or complex problems, SIGMA may also be:
%      'LR' or 'SR' - Largest or Smallest Real part
%      'LI' or 'SI' - Largest or Smallest Imaginary part
%   If SIGMA is a real or complex scalar including 0, EIGS finds the
%   eigenvalues closest to SIGMA.
%
%   EIGS(A,K,SIGMA,OPTS) and EIGS(A,B,K,SIGMA,OPTS) specify options:
%   OPTS.issym: symmetry of A or A-SIGMA*B represented by AFUN [{false} |
%   true]
%   OPTS.isreal: complexity of A or A-SIGMA*B represented by AFUN [false | {true}]
%   OPTS.tol: convergence: Ritz estimate residual <= tol*NORM(A) [scalar | {eps}]
%   OPTS.maxit: maximum number of iterations [integer | {300}]
%   OPTS.p: number of vectors after truncation: K+1<p<=N [integer | {2K}].
%           Full search space size is 2*p.
%   OPTS.v0: starting vector [N-by-1 vector | {randomly generated}]
%   OPTS.disp: diagnostic information display level [{0} | 1 | 2]
%   OPTS.cholB: B is actually its Cholesky factor CHOL(B) [{false} | true]
%   OPTS.permB: sparse B is actually CHOL(B(permB,permB)) [permB | {1:N}]
%   OPTS.nvlock: check convergence, deflate & lock Ritz vectors by 'nvlock' vectors at once [integer |{K}]. 
%   Use CHOL(B) instead of B when SIGMA is a string other than 'SM'.
%
%   EIGS(AFUN,N) and EIGS(AFUN,N,B) accept the function AFUN instead of the
%   matrix A. AFUN is a function handle and Y = AFUN(X) should return
%      A*X            if SIGMA is unspecified, or a string other than 'SM'
%      A\X            if SIGMA is 0 or 'SM'
%      (A-SIGMA*I)\X  if SIGMA is a nonzero scalar (standard problem)
%      (A-SIGMA*B)\X  if SIGMA is a nonzero scalar (generalized problem)
%   N is the size of A. The matrix A, A-SIGMA*I or A-SIGMA*B represented by
%   AFUN is assumed to be real and nonsymmetric unless specified otherwise
%   by OPTS.isreal and OPTS.issym. 
%
%   EIGS(AFUN,N,...) is equivalent to EIGS(A,...) for all previous syntaxes.
%
%   Example:
%      A = delsq(numgrid('C',15));  d1 = eigs(A,5,'SM');
%
%   Equivalently, if dnRk is the following one-line function:
%      %----------------------------%
%      function y = dnRk(x,R,k)
%      y = (delsq(numgrid(R,k))) \ x;
%      %----------------------------%
%
%      n = size(A,1);  opts.issym = 1;
%      d2 = eigs(@(x)dnRk(x,'C',15),n,5,'SM',opts);
%
%   See also EIG, SVDS, FUNCTION_HANDLE.

%   Copyright 2006-2022 Advanpix, LLC.


%%
% TODO:
%
% * CONVERGENCE SPEED:
% - Use convergence acceleration for Ritz vectors [19].
% - Chebyshev filtering & acceleration [15].
% - Use Cholesky (or LDLT) factorization for symmetric matrices.
% - Use B-inner product (v'Bv) in case of Hermitian generalized problem
%   which keeps the symmetry.
%
% * TECHNICAL SPEEDUP:
% -> (DONE) Implement Gram-Schmidt in C++.
% - Implement the whole thing in C++: 
%   A-la EIGS + ARPACK style - with data exchange using working buffers.
%   The only thing should be done in MATLAB language is application of
%   input operator to vector.
%   This is not of high priority since the heaviest parts are matrix
%   factorization/operator application, Gram-Schmidt, schur and ordschur.
%
% * USABILITY:
% - Handle warnings from ordschur, etc., e.g.:
%   'reordering failed because some eigenvalues are too close to swap.'
%
% * STABILITY:
% - Scaling in case of unsymmetric generalized problem.

% Brief overview of Krylov projection Methods [6]:
%
% Projection methods for eigenvalue problems are intended for computing a partial eigensolution,
% that is, given a square matrix A of order n, the objective is to compute a small number of
% eigenpairs, (lambda(i), x(i)), i = 1...k, k << n.
%
% The basic principle of projection methods is to find the best approximations to the eigenvectors 
% in a given subspace of small dimension. This is achieved with the so-called Rayleigh-Ritz procedure, described next.
%
% Given an n x m matrix V, with k <= m << n, whose columns v(i) constitute an orthonormal
% basis of a given subspace W, i.e. V'V = I and span{v(1), v(2), ..., v(m)} = W, the Rayleigh-Ritz
% procedure consists in computing H = V'AV (Rayleigh quotient), which is a square matrix of order m, and then
% solving the eigenvalue problem associated with it:
%
%                          Hy(i) = theta(i)y(i) 
%
% Then the approximate eigenpairs lambda(i) and x(i) of the original eigenvalue problem are
%
%                  lambda(i) = theta(i) and x(i) = Vy(i)
%
% which are called Ritz values and Ritz vectors, respectively. 
%
% Note that Ritz vectors belong to subspace W. 
% It can be shown that the Ritz values and vectors are the best possible 
% approximations in that subspace. Matrix H is the projection of A onto W, 
% hence the name of this class of methods.
%
% References:
%
% [1] G. Stewart. A Krylov–Schur Algorithm for Large Eigenproblems. SIAM J. Matrix Anal.
%     Appl., 23(3):601–614, 2001.
%
% [2] G. Stewart. Addendum to "A Krylov–Schur Algorithm for Large Eigenproblems". SIAM J.
%     Matrix Anal. Appl., 24(2):599–601, 2002.
%
% [3] D. Kressner. Numerical Methods for General and Structured Eigenvalue Problems , volume 46
%     of Lecture Notes in Computational Science and Engineering. Springer,
%     Berlin, 2005.
%
% [4] (ARPACK) R. Lehoucq, D. Sorensen, and C. Yang. ARPACK Users’ Guide, Solution of Large-Scale
%     Eigenvalue Problems by Implicitly Restarted Arnoldi Methods . Society for Industrial and Applied
%     Mathematics, Philadelphia, PA, 1998.
%
% [5] Y. Saad. Numerical Methods for Large Eigenvalue Problems, SIAM 2011.
%
% [6] (SLEPc) V. Hernandez, J. E. Roman, A. Tomas, V. Vidal. Krylov-Schur Methods in SLEPc. 
%     SLEPc Technical Report STR-7.
%
% [7] Z. Bujanovic. Krylov Type Methods for Large Scale Eigenvalue
%     Computations. Ph.D. thesis, University of Zagreb
% 
% [8] (ANASAZI) C. G. Baker, U. L. Hetmaniuk, R. B. Lehoucq, H. K. Thornquist.  Anasazi software for the 
%     numerical solution of large-scale eigenvalue problems, 
%     ACM Transactions on Mathematical Software, 
%     Volume 36 Issue 3, July 2009.
% 
% [9] (MultiParEig) A. Muhic, B. Plestenjak. On the quadratic two-parameter eigenvalue problem 
%     and its linearization, Linear Algebra Appl. 432 (2010) 2529-2542. 
%     https://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig
%
% [10] X. Ding. Krylov-Schur method, https://github.com/dingxiong/KrylovSchur
%
% [11] R. Radke. A Matlab implementation of the Implicitly Restarted Arnoldi 
%      Method for solving large-scale eigenvalue problems, (1996) Master’s Thesis, Rice University. 
%      http://hdl.handle.net/1911/14054.
% 
% [12](SLEPc) J. E. Roman, Practical Implementation of Harmonic Krylov-Schur
%     SLEPc Technical Report STR-9. http://slepc.upv.es/documentation/reports/str9.pdf
%
% [13] Parlett, B. N.: The symmetric Eigenvalue Problem, Englewood Cliffs, N.
%      J. Prentice-Hall (1980).
%
% [14] C. Hegedus, Reorthogonalization Methods Revisited.
%      http://www.uni-obuda.hu/journal/Hegedus_64.pdf
%
% [15] J. A. Scott, An Arnoldi Code for Computing Selected Eigenvalues of Sparse, Real, Unsymmetric Matrices.
%      Article?in?ACM Transactions on Mathematical Software 21(4):432-475,
%      December 1995.
%
% [16] R. Lehoucq and J. A. Scott, An evaluation of software for computing eigenvalues of sparse nonsymmetric matrices, 
%      Preprint MCS-P547-1195, Argonne National Laboratory.
%      https://goo.gl/F5IPLx
%
% [17] M. Bennani and T. Braconnier. Stopping criteria for eigensolvers. Technical report, 
%      CERFACS, November 1994.
%
% [18] F. Chatelin. Eigenvalues of Matrices. Wiley, 1993.
%
% [19] Sequence convergence acceleration techniques:
%      Wynn's epsilon algorithm, very nice and tested in SIAM100 tasks to
%      eccelerate convergence of eigenvalues. 
%
%      General information on acceleration/extrapolation:
%      - Hamming. Numerical Methods for Scientists and Engineers, p206-207. 
%      - https://en.wikipedia.org/wiki/Shanks_transformation
%      - Comprehensive review: https://arxiv.org/abs/math.NA/0306302v1
%
%      Acceleration/extrapolation for vector iterative processes:
%      - Avram Sidi. Vector extrapolation methods with applications to solution 
%        of large systems of equations and to PageRank computations:
%        http://www.cs.technion.ac.il/users/wwwb/cgi-bin/tr-get.cgi/2007/CS/CS-2007-09.revised.pdf
%
%      - Sadok. Vector extrapolation methods. Applications and numerical comparison
%        http://www.sciencedirect.com/science/article/pii/S0377042700003575
%
%      - Sadok. Nonlinear Schwarz iterations with Reduced Rank Extrapolation (RRE):
%        https://math.temple.edu/~szyld/reports/RRE.Schwarz.report.pdf
%
%      - p39-43: https://people.cs.kuleuven.be/~raf.vandebril/bari2007/talks/thu1.pdf
%
%    Copyright (c) 2006 - 2022 Advanpix LLC.

    if (nargout > 4)
        error(message('MATLAB:eigs:TooManyOutputs'));
    end
    
    % Platform dependent integer type
    if strfind(computer, '64')
        intconvert = @(arraytoconvert) int64(arraytoconvert);
        inttype = 'int64';
    else
        intconvert = @(arraytoconvert) int32(arraytoconvert);
        inttype = 'int32';
    end
    
    % We borrow parsing of arguments from original EIGS - code is quite
    % convoluted and ugly, but we use it for the moment.
    % Error check inputs and derive some information from them
    [A,Amatrix,isrealprob,issymA,n,B,class_t,k,eigs_sigma,whch, ...
        sigma,usertol,maxit,basis,info,verbose,cholB,permB,v0,useeig, ...
        afunNargs,style,mode,nvlock,Afactors,Bfactors] = ...
        checkInputs(varargin{:});

    % Now have enough information to do early return on cases EIGS does not
    % handle. For these cases, use the full EIG code.
    if useeig
        fullEig(nargout);
        return
    end
    
    precompLU = [];
    if ~isempty(Afactors)
        LA = Afactors.L;     Afactors.L = [];
        UA = Afactors.U;     Afactors.U = [];
        pp = Afactors.pp;   Afactors.pp = [];
        qq = Afactors.qq;   Afactors.qq = [];
        dgAsB = Afactors.dgAsB;    Afactors.dgAsB = [];
        precompLU = Afactors.precompLU; Afactors.precompLU = [];
        clear Afactors;
    end

    precompLUB = [];
    if ~isempty(Bfactors)
        BisHpd = Bfactors.BisHpd;
        if BisHpd
            RB = Bfactors.RB;       Bfactors.RB = [];
            RBT = Bfactors.RBT;     Bfactors.RBT = [];
            permB = Bfactors.permB; Bfactors.permB = [];
        else
            LB = Bfactors.LB;       Bfactors.LB = [];
            UB = Bfactors.UB;       Bfactors.UB = [];
            ppB = Bfactors.ppB;     Bfactors.ppB = [];
            qqB = Bfactors.qqB;     Bfactors.qqB = [];
            dgB = Bfactors.dgB;     Bfactors.dgB = [];
            precompLUB = Bfactors.precompLUB; Bfactors.precompLUB = [];            
        end
        clear Bfactors;
    end
    
    % Modes & specral transformations:
    % LM, LI, LR, SI, SR, SA, LA - no transformation
    % SM - invert transformation
    % LM + shift-and-invert if sigma is numeric scalar
    if isempty(B)
        if mode ~= 3
            % OP = A
            applyOP = @(v)Amtimes(v);
        else
            % OP = (A-\sigma*I)^{-1}
            applyOP = @(v)AminusSigmaBsolve(v);
        end
    else
        if mode ~= 3
            if BisHpd == true
                % OP = L^{-1}AL^{-T} (B = LL^T)
                applyOP = @(v)RBTsolve(Amtimes(RBsolve(v)));
            else
                % OP = U^{-1}L^{-1}A (B = LU)
                applyOP = @(v)Bsolve(Amtimes(v));
            end
        else
            % OP = (A-\sigma*B)^{-1}B
            applyOP = @(v)AminusSigmaBsolve(Bmtimes(v));
        end
    end
    
    applyM = @(v) v;
    if strcmp(style,'G') && ~isempty(B) && BisHpd == true
        applyM = @(v) Bmtimes(v);
    end

    basis = min(basis,max(1,ceil(n/4)));
    b = basis;   % number of vectors kept after restart (we call it a basis)             
    m = 2*basis; % full set of vectors in search subspace, 'basis' is number of vectors we leave after truncation.

    if verbose
        fprintf('\n   Krylov-Schur Iterations...\n');
    end
    
    Q = zeros(n, m+1, class_t); 
    H = zeros(m+1, m, class_t);
    Q(:,1) = v0 / norm(v0);
    p = 1;                     % indicates number of already converged vectors (= p-1).
    
    % Compute initial Krylov subspace
    [Q, H, stopAlgorithm] = expandKrylov(applyOP, Q, H, 1, b); % adds vectors 2:b+1, returns Q(n, b+1) and H(b+1,b)

    if stopAlgorithm
        % Failure to build an orthogonal subspace
        error('Unable to build an orthogonal subspace because the problem is ill-conditioned. For generalized problems, B may have low rank.');
    end
        
    iterations = 0;
    while iterations < maxit && p <= k

        iterations = iterations+1;
        
        % Expand Krylov subspace using Arnoldi process.
        [Q, H, stopAlgorithm] = expandKrylov(applyOP, Q, H, b+1, m); % adds vectors b+2:m+1, returns Q(n, m+1) and H(m+1,m)
        
        if stopAlgorithm
            % Failure to build an orthogonal subspace
            error('Unable to build an orthogonal subspace because the problem is ill-conditioned. For generalized problems, B may have low rank.');
        end
        
        % Compute its Krylov-Schur decomposition and Ritz values.
        [U, T] = schur(H(p:m, p:m));  % produces real quasi-triangluar Schur form if H is real.
        ritz   = ordeig(T);        
        
        % Sort Ritz values according to mode
        if strcmpi(whch,'LM') || strcmpi(whch,'SM')
            [~, order] = sort(abs(ritz), 'descend');            
        elseif strcmpi(whch,'LA') || strcmpi(whch,'LR')
            [~, order] = sort(real(ritz),'descend');            
        elseif strcmpi(whch,'SA') || strcmpi(whch,'SR')
            [~, order] = sort(real(ritz),'ascend');
        elseif strcmpi(whch,'LI')
            [~, order] = sort(imag(ritz),'descend');            
        elseif strcmpi(whch,'SI')
            [~, order] = sort(imag(ritz),'ascend');
        end;

        % Prepare clusters markers for reordering. 
        % Ritz values with the highest cluster marker go to left top corner.
        nritz = m-p+1;  % total number of Ritz values at hand
        clusters = zeros(nritz,1);
        clusters(order) = nritz:-1:1;
        
        % Make sure conjugate pairs in real Schur form fall in the same (and unique) cluster
        if isreal(T)
            for j = 1:nritz-1
                if T(j+1,j) ~= 0 
                    clusters(j+1) = clusters(j);
                end;
            end
        end

        % Re-order and update Krylov-Schur according to sorted order of Ritz values    
        [U, T] = ordschur(U, T, clusters);

        Q(:,     p:m) = Q(:, p:m) * U;        
        H(p:m,   p:m) = T;
        H(1:p-1, p:m) = H(1:p-1, p:m) * U;
        H(m+1,   p:m) = H(  m+1, p:m) * U;
        
        % Truncate Krylov-Schur decomposition
        % Make sure conjugate Ritz are not splitted on truncation.
        b = basis + (H(basis+1,basis) ~= 0);
        
        Q = [Q(  :, 1:b)  Q(  :, m+1)];
        H = [H(1:b, 1:b); H(m+1, 1:b)];

        %% Check convergence, do the deflation and locking.
        %
        % Deflation & locking reduces number of iterations and beneficial 
        % for computing independent vectors for multiple eigenvalue,
        % see [1], section 5. 
        %
        % There is a chance that some Ritz vectors will be locked earlier than wanted ones appeared.
        % Wanted Ritz vectors might require some number of iterations to appear in basis.
        % Locking might already lock "unwanted" vectors with small Ritz
        % estimate long before wanted vectors appear in the basis and start
        % to converge. This is especially troublesome for non-normal
        % operators, where small Ritz estimate doesn't guarantee
        % high accuracy of approximation for original eigenpair. 
        % 
        % Similar findings were reported in papers and other libraries. ARPACK uses locking of eigenvectors 
        % but it runs additional verification steps since incorrect vectors can be locked instead of correct ones.
        % As extra protection, it checks convergence of set of vectors.
        % Same was reported by Scott as for EB13 solver [15] - unwanted vectors might be locked early 
        % and spoil all the fun. 
        %
        % TODO: 
        % - More intelligent locking to include wanted Ritz pairs which
        % converged late. Once new vector converged, we need to check whether it 
        % should be inserted before some already locked Ritz pairs.
        % 
        % Convergence detection and stopping criteria in iterative eigensolvers 
        % are not trivial questions (especially in unsymmetric case).
        % We collect important facts in notes nelow.
        %
        
        %% I. Arnoldi factorization (ARPACK, EB13 [15], etc.)
        %
        % Having the Arnoldi factorization of length m:
        %
        %     A*Q_m = Q_m*H_m + f_m*e_m',  f = Q(:,m+1) - last computed vector. 
        %                                  f_m*e_m'     - its last component.
        %                                  H_m          - upper Hessenberg matrix.
        %
        % and Ritz pair {theta(i),y(i)} of Rayleigh quotient H_m:
        %
        %     H_m*y(i) = theta(i)*y(i), ||y(i)|| = 1
        %
        % Then approximation accuracy of the eigenpairs {theta(i), Q_m*y(i)} 
        % of original problem can be estimated as:
        %
        %    ||A*x(i)-theta(i)*x(i)|| = 
        %                             = ||A*Q_m*y(i)-theta(i)*Q_m*y(i)|| = 
        %                             = ||(A*Q_m-Q_m*H_m)*y(i)|| = ||f_m*e_m'*y(i)||
        %                             = ||f_m||*|e_m'*y(i)|, beta = ||f_m||
        %
        % Quantity ||f_m||*|e_m'*y(i)| is called Ritz estimate and it can
        % be used to detect convergence of Ritz eigenpair {theta(i), y(i)}.
        % In particular ARPACK terminates computations when all Ritz values
        % satisfy (see [15] + ARPACK source code):
        % 
        %            ||f_m||*|e_m'*y(i)| < |theta(i)| * tol
        %
        % In the Hermitian case, this estimate on the residual can be 
        % turned into a precise statement about the accuracy 
        % of the Ritz value theta(i) as an approximation to the eigenvalue of
        % A that is nearest to theta(i). 
        %
        % However, an analogous statement in the non-Hermitian case 
        % is not possible without further information concerning non-normality 
        % and defectiveness. If A is highly non-normal, small Ritz estimate 
        % doesn't necessarily mean that actual eigenvector residual is small too [17], [18]. 
        % See also [15-16] for good overviews of the issues related to stopping criteria.
        %
        % As the algorithm converges, Ritz estimate become poor estimates 
        % of the actual error in each Ritz pair [11]. Thus, if guaranteed
        % accuracy is required, it is advisable to use Ritz estimate only 
        % at the begining of the process and switch to direct error evaluation 
        % when Ritz estimate becomes small. Direct error evaluation
        %
        %             ||A*x(i)-theta(i)*x(i)||, x(i) = Q*y(i)
        %
        % is expensive process but ensures that the user gets the desired
        % eigenpairs to the requested tolerance.
        %
        %% I. Krylov-Schur factorization
        %
        % Krylov generalizes Arnoldi decomposition into: 
        %
        %     A*V_m = V_m*B_m + v_{m+1}*b_{m+1}' 
        % 
        % B_m is not restricted to be upper Hessenberg and b_{m+1} is an
        % arbitrary vector. Special (and very useful) case is Krylov-Schur
        % decomposition, in which B_m is in real Schur form (see [6] for
        % nice schemes).
        %
        % Analogously to Arnoldi, having Ritz pair {theta(i),y(i)} of
        % Rayleigh quotient B_m 
        %
        %     B_m*y(i) = theta(i)*y(i), ||y(i)|| = 1
        %
        % we can estimate approximation accuracy of the eigenpairs {theta(i), V_m*y(i)} 
        % of original problem as:
        %
        %    ||A*x(i)-theta(i)*x(i)|| = 
        %                             = ||A*V_m*y(i)-theta(i)*V_m*y(i)|| = 
        %                             = ||(A*V_m-V_m*B_m)*y(i)|| = 
        %                             = ||v_{m+1}*b_{m+1}'*y(i)||.
        %
        % Having ||v_{m+1}|| = ||y(i)|| = 1, we finally arrive to:
        %
        %    ||A*x(i)-theta(i)*x(i)|| = ||b_{m+1}||.
        %
        % Thus, in case of Krylov-Schur convergence check is simplier.
        %
        
        
        %% Check convergence, do deflation & locking if needed.
        %
        % Convergence detection in Krylov-Schur is outlined in details 
        % above.
        %
        % We check residual for every Ritz value (last row of H = b_{m+1})
        % as convergence indicator. Please note, that if nlock < k,
        % then some of the "unwanted" Ritz values might buble up 
        % and converged much earlier than wanted ones. That is why
        % in locking we need to check if sorting criteria is not
        % broken (e.g. newly found Ritz value is larger than already 
        % locked one in 'LM' mode).

        r = 1:k+(H(k+1,k) ~= 0);
        D = ordeig(H(r,r));
        
        safetol = norm(H(r,r),'fro') * numeric_t('eps',class_t); % safeguard tolerance
        
        if verbose
            % Compute residual for each Ritz value before deflation
            % Will be needed for detailed report later.
            eachRitzResidual = zeros(length(D),1,class_t);

            w = p; % save value of 'p' before deflation & locking
            u = 1; % compute residual for all Ritz pairs we have at hand.
            while u <= k
                conjugateRitzValues = (H(u+1,u) ~= 0);
                t = u + conjugateRitzValues;
                eachRitzResidual(u) = norm(H(b+1,u:t));
                if conjugateRitzValues
                    eachRitzResidual(u+1) = eachRitzResidual(u);
                end
                u = t+1;
            end
        end;

        %% Deflation & Locking
        % Locking reduces number of iterations.
        % In 50 tests of finding 8-eigenpairs of sparse random 1000x1000
        % matrix with non-zeros density of 0.01, locking (nlock = 1) needs 
        % less iterations by 30% compared to no locking (nlock = 8).
        
        % TODO:
        % There is no guarantee that Ritz vectors will appear in order we
        % need (e.g. from largest to smallest in case of 'L*').
        % Out-of-order vectors might buble up and converge long before
        % wanted vectors appear in the basis. To prevent this, we have to do simple
        % check and reset the locking position to preserve sorting order
        % of Ritz values. It would be more intelligent to do the
        % deflation & locking before Schur re-ordering, taking into account
        % all these effects. This task is for the future.
        
        converged = true;
        while p <= k && converged
            t = min(p-1+nvlock,k);     % last index of next block to check
            t = t + (H(t+1,t) ~= 0);   % include conjugate Ritz values if any

            residual = norm(H(b+1,p:t));

            converged = residual <= max(numeric_t('eps',class_t), max(usertol*norm(D(p:t)),safetol));
            if converged
                H(b+1,p:t) = 0;        % deflate converged vector(s) by zeroing their residuals
                p = t+1;               % and lock them.
            end;
        end
        
        if verbose 
            fprintf('\nIteration %d: Ritz values with residual, convergence tolerance and status:\n', iterations);
            
            % Display true Ritz values, after shift-and-invert transformation. 
            if mode ~= 3
                d = D;
            else
                if ~ischar(eigs_sigma)
                    d = eigs_sigma + 1./D;
                else
                    d = 1./D;
                end
            end
            
            % number of digits after comma
            if strcmpi(class_t,'mp')
                precision = mp.Digits();
            else
                precision = 5;
            end                
            width = precision+5; % total number of characters for number
            format = ['%' num2str(width) '.' num2str(precision) 'g'];
            schar  = ['-', '+', '+'];
            
            u = 1;
            while u <= k
                conjugateRitzValues = (H(u+1,u) ~= 0);
                t = u + conjugateRitzValues;

                tolerance = max(numeric_t('eps^(6/7)',class_t), max(usertol*norm(D(u:t)),safetol));

                if u < w
                    status = 'converged/locked';
                elseif u < p
                    status = 'just converged';                                
                else
                    status = '-';                            
                end;

                if isreal(d),  fprintf([format '%15.5e\t%15.5e\t(%s)\n'], d(u), eachRitzResidual(u), tolerance, status);
                else,          fprintf([format ' %c' format 'i' '%15.5e\t%15.5e\t(%s)\n'], real(d(u)), schar(sign(imag(d(u)))+2), abs(imag(d(u))), eachRitzResidual(u), tolerance, status);                    
                end
                
                if conjugateRitzValues
                    if isreal(d),  fprintf([format '%15.5e\t%15.5e\t(%s)\n'], d(u+1), eachRitzResidual(u), tolerance, status);
                    else,          fprintf([format ' %c' format 'i' '%15.5e\t%15.5e\t(%s)\n'], real(d(u+1)), schar(sign(imag(d(u+1)))+2), abs(imag(d(u+1))), eachRitzResidual(u), tolerance, status);                    
                    end
                end;
                
                u = t+1;
            end
        end;
    end  % main loop
    
    if ~isempty(precompLU)
        delete(precompLU);
    end;
    
    if ~isempty(precompLUB)
        delete(precompLUB);
    end;
        
    % Compute Ritz eigenpair of original matrix A from Krylov-Schur decomposition
    % using Rayleigh-Ritz procedure:
    r = 1:k+(H(k+1,k) ~= 0);
    if nargout <= 1
        d = ordeig(H(r,r));
        if mode ~= 3
            varargout{1} = d(1:k);
        else
            if ~ischar(eigs_sigma)
                varargout{1} = eigs_sigma + 1./d(1:k);            
            else
                varargout{1} = 1./d(1:k);                            
            end
        end
    else
        [V,D] = eig(H(r,r));
        V = Q(:,r)*V;
        
        varargout{1} = V(:,1:k);
        if mode ~= 3
            varargout{2} = D(1:k,1:k);
        else
            if ~ischar(eigs_sigma)
                varargout{2} = diag(eigs_sigma + 1./diag(D(1:k,1:k)));            
            else
                varargout{2} = diag(1./diag(D(1:k,1:k)));                            
            end
        end
    end

    flag = double(p <= k);
    if (flag && nargout < 3)
        if (p-1 == 0)
            warning('No eigenvalues converged');
        else
            warning('%d of the %d requested eigenvalues converged. Specify convergence flag output to determine the values that did not converge.', p-1, k);
        end
    end

    if nargout >= 3
        varargout{3} = flag;
    end

    if nargout >= 4
        varargout{4} = iterations;
    end

%% Nested functions
% Some of the functions are borrowed from original EIGS. 

% checkInputs error checks the inputs to EIGS and also derives some
%   variables from them:
% A may be a matrix or a function applying OP.
% Amatrix is true if A is a matrix, false if A is a function.
% isrealprob is true if all of A, B and sigma are real, false otherwise.
% issymA is true if A is symmetric, false otherwise.
% n is the size of (square) A and B.
% B is [] for the standard problem. Otherwise it may be one of B, CHOL(B)
%   or CHOL(B(permB,permB)).
% classAB is single if either A or B is single, otherwise double.
% k is the number of eigenvalues to be computed.
% eigs_sigma is the value for sigma passed in by the user, 'LM' if it was
%   unspecified. eigs_sigma may be either a string or a scalar value.
% whch is the ARPACK string corresponding to eigs_sigma and mode.
% sigma is the ARPACK scalar corresponding to eigs_sigma and mode.
% tol is the convergence tolerance.
% maxit is the maximum number of iterations.
% p is the number of Lanczos vectors.
% info is the start value, initialized to 1 or 0 to indicate whether to use
% resid as the start vector or not.
% eigs_display is true if Ritz values should be displayed, false otherwise.
% cholB is true if CHOL(B) was passed in instead of B, false otherwise.
% permB may be [], otherwise it is the permutation in CHOL(B(permB,permB)).
% resid is the start vector if specified and info=1, otherwise all zero.
% useeig is true if we need to use EIG instead of ARPACK, otherwise false.
% afunNargs is the range of EIGS' varargin that are to be passed as
%   trailing parameters to the function as in afun(X,P1,P2,...).
    function [A,Amatrix,isrealprob,issymA,n,B,classAB,k, ...
            eigs_sigma,whch,sigma,tol,maxit,p,info,eigs_display,cholB,...
            permB,resid,useeig,afunNargs,style,mode,nvlock,Afactors,Bfactors] = checkInputs(varargin)
        % Process inputs and do error-checking
        
        % #1 Argument
        % Process the input A or the inputs AFUN and N
        % Start to derive some qualities (real, symmetric) about the problem
        if isfloat(varargin{1})
            A = varargin{1};
            Amatrix = true;
        else
            % By checking the function A with fcnchk, we can now use direct
            % function evaluation on the result, without resorting to feval
            [A,notFunc] = fcnchk(varargin{1});
            Amatrix = false;
            if ~isempty(notFunc)
                error(message('MATLAB:eigs:NonDoubleOrFunction'));
            end
        end
        
        isrealprob = true;
        issymA = false;
        if Amatrix
            isrealprob = isreal(A);
            issymA = ishermitian(A);
            [m,n] = size(A);
            if (m ~= n)
                error(message('MATLAB:eigs:NonSquareMatrixOrFunction'));
            end
        else
            n = varargin{2};
            if ~isscalar(n) || ~isreal(n) || (n<0) || ~isfinite(n)
                error(message('MATLAB:eigs:NonPosIntSize'));
            end
            if issparse(n)
                n = full(n);
            end
            if (round(n) ~= n)
                warning(message('MATLAB:eigs:RoundNonIntSize'));
                n = round(n);
            end
        end
        
        % #2 Argument        
        % Process the input B and derive the class of the problem.
        % Is B present in the eigs call or not?
        Bpresent = true;
        if (nargin < (3-Amatrix))
            B = [];
            Bpresent = false;
        else
            % Is the next input B or K?
            B = varargin{3-Amatrix};
            if ~isempty(B) % allow eigs(A,[],k,sigma,opts);
                if isscalar(B)
                    if n ~= 1
                        % this input is really K and B is not specified
                        B = [];
                        Bpresent = false;
                    else
                        % This input could be B or K.
                        % If A is scalar, then the only valid value for k is 1.
                        % So if this input is scalar, let it be B, namely
                        % eigs(4,2,...) assumes A=4, B=2, NOT A=4, k=2
                        if ~isnumeric(B)
                            error(message('MATLAB:eigs:BsizeMismatchA'));
                        end
                        % Unless, of course, the scalar is 1, in which case
                        % assume that it is meant to be K.
                        if (B == 1) && ((Amatrix && nargin <= 4) || ...
                                (~Amatrix && nargin <= 5))
                            B = [];
                            Bpresent = false;
                        elseif ~isfloat(B)
                            error(message('MATLAB:eigs:BsizeMismatchA'));
                        end
                    end
                else
                    % B is a not a scalar.
                    if ~isfloat(B) || ~isequal(size(B),[n,n])
                        error(message('MATLAB:eigs:BsizeMismatchA'));
                    end
                    isrealprob = isrealprob && isreal(B);
                end
            end
        end
        % We can only handle homogeneous inputs
        if Amatrix
            classAB = 'mp'; % superiorfloat(A,B);
            A = cast(A,classAB);
            B = cast(B,classAB);
        else
            if ~isempty(B)
                classAB = 'mp';% class(B);
            else
                classAB = 'mp';
            end
        end
        
        % argOffset tells us where to get the eigs inputs K, SIGMA and OPTS.
        % If A is really the function afun, then it also helps us find the
        % trailing parameters in eigs(afun,n,[B],k,sigma,opts,P1,P2,...)
        % Values of argOffset:
        %  0: Amatrix is false and Bpresent is true:
        %     eigs(afun,n,B,k,sigma,opts,P1,P2,...)
        %  1: Amatrix and Bpresent are both true, or both false
        %     eigs(A,B,k,sigma,opts)
        %     eigs(afun,n,k,sigma,opts,P1,P2,...)
        %  2: Amatrix is true and Bpresent is false:
        %     eigs(A,k,sigma,opts)
        argOffset = Amatrix + ~Bpresent;
        
        if Amatrix && ((nargin - Bpresent)>4)
            error(message('MATLAB:eigs:TooManyInputs'));
        end
        
        % Process the input K.
        if (nargin < (4-argOffset))
            k = min(n,6);
        else
            k = varargin{4-argOffset};
            if ~isnumeric(k) || ~isscalar(k) || ~isreal(k) || (k>n) || ...
                    (k<0) || ~isfinite(k)
                if isnumeric(k) && isscalar(k)
                    error(message('MATLAB:eigs:NonIntegerEigQtyDetail', n, num2str(k)));
                elseif ischar(k)
                    error(message('MATLAB:eigs:NonIntegerEigQtyDetail', n, ['''' k '''']));
                elseif isstruct(k)
                    error(message('MATLAB:eigs:NonIntegerEigQtyStruct', n));
                else
                    error(message('MATLAB:eigs:NonIntegerEigQty', n));
                end
            end
            if issparse(k)
                k = full(k);
            end
            if (round(k) ~= k)
                warning(message('MATLAB:eigs:RoundNonIntegerEigQty'));
                k = round(k);
            end
        end
        
        % Process the input SIGMA and derive ARPACK values whch and sigma.
        % eigs_sigma is the value documented in the help as "SIGMA" that is
        % passed in to EIGS. eigs_sigma may be either a scalar, including 0,
        % or a string, including 'SM'.
        % In ARPACK, eigs_sigma corresponds to two variables:
        % 1.  which, called "whch" to avoid conflict with MATLAB's function
        % 2.  sigma
        % whch is always a string. sigma is always a scalar.
        % Valid combinations are shown below. Note eigs_sigma = 0/'SM' has
        % the same sigma/whch values as eigs_sigma='LM' (default) so these
        % must be distinguished by the mode.
        % eigs_sigma = 'SM' or 0 => sigma = 0, whch = 'LM' (mode=3)
        % eigs_sigma is a string not 'SM' => sigma = 0, whch = eigs_sigma (mode=1)
        % eigs_sigma is any scalar => sigma = eigs_sigma, whch = 'LM'
        % (mode=1)
        if (nargin < (5-argOffset))
            % default: eigs 'LM' => ARPACK which='LM', sigma=0
            eigs_sigma = 'LM';
            whch = 'LM';
            sigma = 0;
        else
            eigs_sigma = varargin{5-argOffset};
            if ischar(eigs_sigma)
                % eigs(string) => ARPACK which=string, sigma=0
                if ~isequal(size(eigs_sigma),[1,2])
                    error(message('MATLAB:eigs:EigenvalueRangeNotValid'));
                end
                eigs_sigma = upper(eigs_sigma);
                if strcmp(eigs_sigma,'SM')
                    % eigs('SM') => ARPACK which='LM', sigma=0
                    whch = 'LM';
                else
                    % eigs(string), where string~='SM' => ARPACK which=string, sigma=0
                    whch = eigs_sigma;
                end
                sigma = zeros(1,classAB);
            else
                % eigs(scalar) => ARPACK which='LM', sigma=scalar
                if ~isfloat(eigs_sigma) || ~isscalar(eigs_sigma)
                    error(message('MATLAB:eigs:EigenvalueShiftNonScalar'));
                end
                sigma = eigs_sigma;
                if issparse(sigma)
                    sigma = full(sigma);
                end
                sigma = cast(sigma,classAB);
                isrealprob = isrealprob && isreal(sigma);
                whch = 'LM';
            end
        end
        
        % Process the input OPTS and derive some ARPACK values.
        % ARPACK's minimum tolerance is eps/2 ([S/D]LAMCH's EPS)
        tol = numeric_t('eps^(6/7)',classAB);
        maxit = [];
        p = [];
        % Always use resid as the start vector, whether it is OPTS.v0 or
        % randomly generated within eigs.  We default resid to empty here.
        % If the user does not initialize it, we provide a random residual
        % below.
        info = intconvert(1);
        resid = [];
        nvlock = [];
        eigs_display = 0;
        fptype = [];
        cholB = false; % do we have B or its Cholesky factor?
        permB = []; % if cholB, is it chol(B), or chol(B(permB,permB))?
        if (nargin >= (6-argOffset))
            opts = varargin{6-argOffset};
            if ~isa(opts,'struct')
                error(message('MATLAB:eigs:OptionsNotStructure'));
            end
            if isfield(opts,'issym') && ~Amatrix
                issymA = opts.issym;
                if (issymA ~= false) && (issymA ~= true)
                    error(message('MATLAB:eigs:InvalidOptsIssym'));
                end
            end
            if isfield(opts,'isreal') && ~Amatrix
                if (opts.isreal ~= false) && (opts.isreal ~= true)
                    error(message('MATLAB:eigs:InvalidOptsIsreal'));
                end
                isrealprob = isrealprob && opts.isreal;
            end
            if ~isempty(B) && (isfield(opts,'cholB') || isfield(opts,'permB'))
                if isfield(opts,'cholB')
                    cholB = opts.cholB;
                    if (cholB ~= false) && (cholB ~= true)
                        error(message('MATLAB:eigs:InvalidOptsCholB'));
                    end
                    if isfield(opts,'permB')
                        if issparse(B) && cholB
                            permB = opts.permB;
                            if ~isvector(permB) || ~isequal(sort(permB(:)),(1:n)')
                                error(message('MATLAB:eigs:InvalidOptsPermB'));
                            end
                        else
                            warning(message('MATLAB:eigs:IgnoredOptionPermB'));
                        end
                    end
                end
            end
            if isfield(opts,'tol')
                if ~isfloat(tol) || ~isscalar(opts.tol) || ~isreal(opts.tol) || (opts.tol<=0)
                    error(message('MATLAB:eigs:InvalidOptsTol'));
                end
                tol = cast(full(opts.tol),classAB);
            end
            
            if isfield(opts,'p')
                p = opts.p;
                if ~isnumeric(p) || ~isscalar(p) || ~isreal(p) || (p<=0) || (p>n)
                    error(message('MATLAB:eigs:InvalidOptsP'));
                end
                if issparse(p)
                    p = full(p);
                end
                if (round(p) ~= p)
                    warning(message('MATLAB:eigs:NonIntegerVecQty'));
                    p = round(p);
                end
            end

            if isfield(opts,'nvlock')
                nvlock = min(opts.nvlock,k);
                if ~isnumeric(nvlock) || ~isscalar(nvlock) || ~isreal(nvlock) || (nvlock<=0) || (nvlock>n)
                    error(message('MATLAB:eigs:InvalidOptsP'));
                end
                if issparse(nvlock)
                    nvlock = full(nvlock);
                end
                if (round(nvlock) ~= nvlock)
                    warning(message('MATLAB:eigs:NonIntegerVecQty'));
                    nvlock = round(nvlock);
                end
            end

            if isfield(opts,'fptype')
                fptype = opts.fptype;
                if ~ischar(fptype)
                    error(message('OPTS.fptype should be a string'));
                end
            end
            
            if isfield(opts,'maxit')
                maxit = opts.maxit;
                if ~isnumeric(maxit) || ~isscalar(maxit) || ~isreal(maxit) ...
                        || (maxit<=0) || ~isfinite(maxit)
                    error(message('MATLAB:eigs:OptsMaxitNotPosInt'));
                end
                if issparse(maxit)
                    maxit = full(maxit);
                end
                if (round(maxit) ~= maxit)
                    warning(message('MATLAB:eigs:NonIntegerIterationQty'));
                    maxit = round(maxit);
                end
            end
            
            if isfield(opts,'v0')
                if ~isfloat(opts.v0) || ~isequal(size(opts.v0),[n,1])
                    error(message('MATLAB:eigs:WrongSizeOptsV0'));
                end
                
                if isrealprob && ~isreal(opts.v0)
                    error(message('MATLAB:eigs:NotRealOptsV0'));
                end
                resid(1:n,1) = cast(full(opts.v0),classAB);                
            end
            
            if isfield(opts,'disp')
                eigs_display = opts.disp;
                if ~isnumeric(eigs_display) || ~isscalar(eigs_display) || ...
                        ~isreal(eigs_display) || (eigs_display<0)
                    error(message('MATLAB:eigs:NonIntegerDiagnosticLevel'));
                end
                if (round(eigs_display) ~= eigs_display)
                    warning(message('MATLAB:eigs:RoundNonIntDiagnosticLevel'));
                    eigs_display = round(eigs_display);
                end
            end
        end
        
        if (isempty(resid))
            resid = randn(n,1,classAB);
        end
        
        if eigs_display
            fprintf('\n Multiprecision Krylov-Schur parameters:\n');
            fprintf('   Number of wanted eigenvalues: %d\n',k);
            fprintf('   User-supplied tolerance: %e\n',tol);
        end
        
        afunNargs = zeros(1,0);
        if ~Amatrix
            % The trailing parameters for afun start at varargin{7-argOffset}
            % in eigs(afun,n,[B],k,sigma,opts,P1,P2,...). If there are no
            % trailing parameters in eigs, then afunNargs is a 1-by-0 empty
            % and no trailing parameters are passed to afun(x)
            afunNargs = 7-argOffset:nargin;
        end
        
        % Now that OPTS has been processed, do final error checking and
        % assign ARPACK variables
        
        % Extra check on input B
        if ~isempty(B)
            % B must be symmetric (Hermitian) positive (semi-)definite
            if cholB
                if ~isequal(triu(B),B)
                    error(message('MATLAB:eigs:BNotChol'));
                end
            end
        end
        
        style = 'S';
        if strcmp(eigs_sigma,'SM') || isscalar(eigs_sigma) || ~isempty(B)
            style = 'G';
        end
        
        % **** DO WE NEED SUCH SCALING AT ALL?
        
%         if ~isempty(B)
%             scaleB = norm(B,'fro')/sqrt(n);
%             scaleB = 2^(floor(log2(scaleB+1)));
%             B = B/scaleB;
%             if cholB
%                 scaleB = scaleB^2;
%             end
%             if isscalar(eigs_sigma)
%                 sigma = scaleB*eigs_sigma;
%             end
%         end
        
        
        if strcmp(eigs_sigma,'SM') || ~ischar(eigs_sigma)
            % eigs(A,B,k,scalarSigma) or eigs(A,B,k,'SM'), B may be []
            % Note: sigma must be real for [s,d]saupd and [s,d]naupd
            % If sigma is complex, even if A and B are both real, we use
            % [c,z]naupd.
            % This means that mode=3 in [s,d]naupd, which has
            % OP = real(inv(A - sigma*M)*M) and B = M
            % reduces to the same OP as [s,d]saupd and [c,z]naupd.
            % A*x = lambda*M*x, M symmetric (positive) semi-definite
            % => OP = inv(A - sigma*M)*M and B = M
            % => shift-and-invert mode
            mode = 3;
        elseif strcmp(style,'S')
            % eigs(A,k,stringSigma) or eigs(A,[],k,stringSigma),
            % stringSigma~='SM'
            % A*x = lambda*x
            % => OP = A and B = I
            mode = 1;
        else
            % eigs(A,B,k,stringSigma), stringSigma~='SM'
            % A*x = lambda*B*x
            % => OP = inv(B)*A and use standard inner product.
            mode = 1;
        end
        
        BisHpd = false;
        % Enable HPD only if CHOL(B) was supplied explicitly.
        % At the moment we have no CHOL nor LDLT for sparse matrices.
        if cholB %|| (~isempty(B) && ishermitian(B)) % <-- enable when it CHOL/LDLT will be availible
            % The reordering permutation permB is [] unless B is sparse
            [RB,RBT,permB,BisHpd] = CHOLfactorB;
            if mode == 3 && ~cholB
                RB = [];    RBT = [];   permB = [];
            end
        end
        
        qqB = [];
        if BisHpd == false && (mode == 1 || mode == 2)
            [LB,UB,ppB,qqB,dgB,precompLUB] = LUfactorB;
        end
        
        Bfactors = [];
        if ~isempty(B)
            if BisHpd == true
                Bfactors = struct('RB',RB,'RBT',RBT,'permB',permB,'BisHpd',BisHpd);
            elseif (mode == 1 || mode == 2)
                Bfactors = struct('LB',LB,'UB',UB,'ppB',ppB,'qqB',qqB,'dgB',dgB,'BisHpd',BisHpd,'precompLUB',precompLUB);
            end
        end
        
        Afactors = [];
        qq = [];
        if (mode == 3) && Amatrix % need lu(A-sigma*B)
            % The reordering permutation permAsB is [] unless A-sigma*B is sparse
            
            if eigs_display
                fprintf('   Compute decomposition of A...\n');
            end
        
            [LA,UA,pp,qq,dgAsB,precompLU] = LUfactorAminusSigmaB;
            Afactors = struct('L',LA,'U',UA,'pp',pp,'qq',qq,'dgAsB',dgAsB,'precompLU',precompLU);
        end % if (mode == 3) && Amatrix
        
        % under these conditions, OP must be unsymmetric
        % note that OP = inv(A-\sigma*B)*B IS symmetric if A is symmetric
        % and B-inner product is used!
        if ~isempty(B) && (BisHpd == false || (strcmp(style,'S') && mode == 3))
            issymA = false;
        end
        % Extra check on input K
        % We fall back on using the full EIG code if K is too large.
        useeig = false;
        if (k == 0)
            useeig = true;
        end
        if isrealprob && issymA
            if (k > n-1)
                if (n >= 6)
                    warning(message('MATLAB:eigs:TooManyRequestedEigsForRealSym'));
                end
                useeig = true;
            end
        else
            if (k > n-2)
                if (n >= 7)
                    warning(message('MATLAB:eigs:TooManyRequestedEigsForComplexNonsym'));
                end
                useeig = true;
            end
        end
        
        % Extra check on input SIGMA
        if issymA && ~isreal(sigma)
            warning(message('MATLAB:eigs:ComplexShiftForHermitianProblem'));
        end
        
        if isrealprob && issymA
            if strcmp(whch,'LR')
                whch = 'LA';
                warning(message('MATLAB:eigs:SigmaChangedToLA'));
            end
            if strcmp(whch,'SR')
                whch = 'SA';
                warning(message('MATLAB:eigs:SigmaChangedToSA'));
            end
            if ~ismember(whch,{'LM', 'SM', 'LA', 'SA', 'BE'})
                error(message('MATLAB:eigs:EigenvalueRangeNotValid'));
            end
        else
            if strcmp(whch,'BE')
                warning(message('MATLAB:eigs:SigmaChangedToLM'));
                whch = 'LM';
            end
            if ~ismember(whch,{'LM', 'SM', 'LR', 'SR', 'LI', 'SI'})
                error(message('MATLAB:eigs:EigenvalueRangeNotValid'));
            end
        end
        
        % The remainder of the error checking does not apply for the large
        % values of K that force us to use full EIG instead of ARPACK.
       
        % Extra check on input OPTS.p
        if isempty(p)
            if isrealprob && ~issymA
                p = min(max(2*k+1,20),n);
            else
                p = min(max(2*k,20),n);
            end
        else
            if isrealprob && issymA
                if (p <= k)
                    error(message('MATLAB:eigs:InvalidOptsPforRealSymProb'));
                end
            else
                if (p <= k+1)
                    error(message('MATLAB:eigs:InvalidOptsPforComplexOrNonSymProb'));
                end
            end
        end

        if p >= n || k == 0
            % Since we are going to build the whole subspace anyway, just use
            % full eig. Ignore starting vector and InnerOpts.
            useeig = true;
        end            
        
        if useeig
            return
        end
        
        if isempty(nvlock)
            nvlock = k; % Check convergence (deflate & lock) all wanted Ritz values at once.
        end
        
        % Extra check on input OPTS.maxit
        if isempty(maxit)
            maxit = max(300,ceil(2*n/max(p,1)));
        end
        
        if eigs_display
            fprintf('   Search subspace size from %d to %d\n',p,2*p); 
            fprintf('   Maximum number of iterations: %d\n',maxit);
        end
        
        % Nested functions in checkInputs
        %-------------------------------------------------------------------------%
        function [RB,RBT,ppB,BisHpd] = CHOLfactorB
            % permB may be [] (from checkInputs) if the problem is not sparse
            % or if it was not passed in as opts.permB
            ppB = permB;
            if cholB
                % CHOL(B) was passed in as B
                RB = B;
                RBT = B';
                BisHpd = true;
            else
                % CHOL(B) was not passed into EIGS
                if ~isempty(B)
                    % Algorithm requires CHOL(B) to be computed
                    if issparse(B)
                        [RB,idxB,ppB] = chol(B,'vector');
                    else
                        [RB,idxB] = chol(B);
                    end
                    if mode == 3
                        RB = [];    ppB = [];
                    end
                    RBT = RB';
                    if (idxB == 0)
                        BisHpd = true;
                    elseif mode == 3 && isreal(B)
                        % LDL decomposition is only for 'SM' eigenvalues of the
                        % pair (A,B) where B is Hermitian positive
                        % semi-definite; in this case, as ARPACK users' guide
                        % suggests, one should still use B-(semi)inner product
                        [~,DB,~] = ldl(B,'vector');
                        BisHpd = checkTridiagForHSD(diag(DB), diag(DB,1));
                    else
                        BisHpd = false;
                    end
                end
            end
            if ~isempty(B) && issparse(B)
                tmp = speye(length(B));
                ppB = tmp(ppB,1:length(B));
            end
        end % CHOLfactorB
        %-------------------------------------------------------------------------%
        function [LB,UB,ppB,qqB,dgB,precompLUB] = LUfactorB
            % LU factor B, including a reordering perm if it is sparse
            precompLUB = [];
            if issparse(B)
                if isa(B,'double')
                    [LB,UB,ppB,qqB,dgB] = lu(B);
                elseif isa(B,'mp')
                    precompLUB = precomputeLU(B);
                    LB = []; UB = []; ppB = []; qqB = []; dgB = [];                
                end
            else
                [LB,UB,ppB] = lu(B,'vector');
                qqB = []; dgB = [];
            end
            
            if isempty(precompLUB)
                % Warn if lu(B) is ill-conditioned
                dUB = diag(UB);
                if any(dUB == 0) || any(diag(LB) == 0)
                    error(message('MATLAB:eigs:SingularB'));
                end
                rcondestUB = full(min(abs(dUB)) / max(abs(dUB)));
                if (rcondestUB < numeric_t('eps',classAB))
                    warning(message('MATLAB:eigs:IllConditionedB', sprintf('%f',rcondestUB)));
                end
            end
        end % LUfactorB
        
        %-------------------------------------------------------------------------%
        function [L,U,pp,qq,dgAsB,precompLU] = LUfactorAminusSigmaB
            % Called only if A is matrix
            % LU factor A-sigma*B, including a reordering perm if it is sparse
            if sigma == 0
                AsB = A;
            elseif isempty(B)
                if issparse(A)
                    AsB = A - sigma * speye(n);
                else
                    AsB = A - sigma * eye(n,classAB);
                end
            else
                if cholB
                    if issparse(B)
                        AsB = A - sigma * checkInputBmtimes(speye(n));
                    else
                        AsB = A - sigma * checkInputBmtimes(eye(n,classAB));
                    end
                else
                    AsB = A - sigma * B;
                end
            end
            
            precompLU = [];
            
            if issparse(AsB)
                if isa(AsB,'double')
                    [L,U,pp,qq,dgAsB] = lu(AsB);
                elseif isa(AsB,'mp')
                    % At the moment, our LU is not very advanced,
                    % without sophisticated P, Q and R. So that 
                    % it frequently gives warning about singular matrix.
                    precompLU = precomputeLU(AsB);
                    L = []; U = []; pp = []; qq = []; dgAsB = [];                
                end
            else
                [L,U,pp] = lu(AsB,'vector');
                qq = []; dgAsB = [];
            end
            
            if isempty(precompLU) % <-- remove aftewards
                % DO WE HAVE TO DO THAT? LU ALREADY CHECKS SINGULARITY.
                % Warn if lu(A-sigma*B) is ill-conditioned
                % => sigma is close to an exact eigenvalue of (A,B)
                dU = diag(U);
                if any(dU == 0) || any(diag(L) == 0)
                    error(message('MATLAB:eigs:AminusBSingular'));
                end
                rcondestU = full(min(abs(dU)) / max(abs(dU)));
                if (rcondestU < numeric_t('eps',classAB))
                    warning(message('MATLAB:eigs:SigmaNearExactEig',sprintf('%f',rcondestU)));
                end
            end
        end % LUfactorAminusSigmaB
        
        %-------------------------------------------------------------------------%
        function v = checkInputBmtimes(u)
            % Matrix-vector multiply v = B*u
            if cholB % use B's cholesky factor and its transpose
                if ~isempty(permB)
                    v = permB'*(RBT * (RB * (permB*u)));
                else
                    v = RBT * (RB * u);
                end
            else
                v = B * u;
            end
        end
        
    end % checkInputs
    
    %% fullEig
    function fullEig(nOutputs)
        % Use EIG(FULL(A)) or EIG(FULL(A),FULL(B)) instead of ARPACK
        if ~isempty(B)
            B = Bmtimes(eye(n,class_t));
            B = full(B); % full(B*scaleB);
        end
        if isfloat(A)
            if issparse(A)
                A = full(A);
            end
        else
            % A is specified by a function.
            % Form the matrix A by applying the function.
            if ischar(eigs_sigma) && ~strcmp(eigs_sigma,'SM')
                % A is a function multiplying A*x
                AA = eye(n,class_t);
                for i = 1:n
                    AA(:,i) = A(AA(:,i),varargin{afunNargs});
                end
                A = AA;
            else
                if (isfloat(eigs_sigma) && eigs_sigma == 0) || strcmp(eigs_sigma,'SM')
                    % A is a function solving A\x
                    invA = eye(n,class_t);
                    for i = 1:n
                        invA(:,i) = A(invA(:,i),varargin{afunNargs});
                    end
                    A = eye(n,class_t) / invA;
                else
                    % A is a function solving (A-sigma*B)\x
                    % B may be [], indicating the identity matrix
                    % U = (A-sigma*B)\sigma*B
                    % => (A-sigma*B)*U = sigma*B
                    % => A*U = sigma*B(U + eye(n))
                    % => A = sigma*B(U + eye(n)) / U
                    if isempty(B)
                        sB = eigs_sigma*eye(n,class_t);
                    else
                        sB = eigs_sigma*B;
                    end
                    U = zeros(n,n,class_t);
                    for i = 1:n
                        U(:,i) = A(sB(:,i),varargin{afunNargs});
                    end
                    A = sB*(U+eye(n,class_t)) / U;
                end
            end
        end
        
%         if isempty(B)
%             eigInputs = {A};
%         else
%             eigInputs = {A,B};
%         end
        
        % Check that values are finite (infs and NaNs cause eig to fail)
        if ~all(isfinite(A),'all') || ~all(isfinite(B),'all')
            error(message('MATLAB:eigs:VeryBadCondition'))
        end

        % Now with full floating point matrices A and B, use EIG:
        if nOutputs <= 1
            if isempty(B)
                d = eig(A);
            else
                d = eig(A, B);
            end
        else
            if isempty(B)
                [V, d] = eig(A, 'vector');
            else
                [V, d] = eig(A, B, 'vector');
            end
        end
        
        % Grab the eigenvalues we want, based on sigma
        if ischar(eigs_sigma)
            if ishermitian(A) && isreal(d) && ismember(eigs_sigma, {'largestimag', 'smallestimag', 'bothendsimag','LI','SI'})
                eigs_sigma = 'largestabs';
            end
            ind = whichEigenvalues(d, eigs_sigma);
        else
            % sigma is a scalar
            [~,ind] = sort(abs(d-eigs_sigma));
        end

        ind = ind(1:k);

        if isequal(eigs_sigma, 'bothendsreal')
            % Choose eigenvalues alternating between top and bottom, but then sort
            % them in ascending direction
            [~, ind2] = sort(real(d(ind)), 'ascend');
            ind = ind(ind2);
        end
        
        if (nOutputs <= 1)
            varargout{1} = d(ind);
        else
            varargout{1} = V(:,ind);
            varargout{2} = diag(d(ind));
            if (nOutputs == 3)
                % flag indicates "convergence"
                varargout{3} = 0;
            end
        end
        
    end % FULLEIG

    function ind = whichEigenvalues(d, method)

    switch method
        case {'largestabs','LM'}
            [~, ind] = sort(abs(d), 'descend');
        case {'smallestabs','SM'}
            [~,ind] = sort(abs(d), 'ascend');
        case {'largestreal','LR'}
            [~, ind] = sort(real(d), 'descend');
        case {'smallestreal','SR'}
            [~, ind] = sort(real(d), 'ascend');
        case {'largestimag','LI'}
            [~, ind] = sort(imag(d), 'descend');
        case {'smallestimag','SI'}
            [~, ind] = sort(imag(d), 'ascend');
        case {'bothendsreal','BE'}
            [~, ind] = sort(real(d), 'descend');
            ind2 = [ind, flip(ind)]';
            ind2 = ind2(:);
            ind = ind2(1:size(d,1));
        case 'bothendsimag'  % todo
            [~, ind] = sort(imag(d), 'descend');
            ind2 = [ind, flip(ind)]';
            ind2 = ind2(:);
            ind = ind2(1:size(d,1));
        case 'smallestimagabs' % todo
            [~,ind] = sort(abs(imag(d)), 'ascend');
    end

    end

%% Operators
    function v = Amtimes(u)
        % Matrix-vector multiply v = A*u
        if Amatrix
            v = A * u;
        else % A is a function
            v = A(u,varargin{afunNargs});
            if isrealprob && ~isreal(v)
                error(message('MATLAB:eigs:complexFunction'));
            end
        end
    end

    function v = Bmtimes(u)
        % Matrix-vector multiply v = B*u
        if cholB % use B's cholesky factor and its transpose
            if ~isempty(permB)
                v = permB'*(RBT * (RB * (permB*u)));
            else
                v = RBT * (RB * u);
            end
        else
            v = B * u;
        end
    end

%-------------------------------------------------------------------------%
    function v = RBsolve(u)
        % Solve v = RB\u for v
        if issparse(B)
            if ~isempty(permB)
                v = permB'*(RB \ u);
            else
                v = RB \ u;
            end
        else
            RBopts.UT = true;
            v = linsolve(RB,u,RBopts);
        end
    end

%-------------------------------------------------------------------------%
    function v = RBTsolve(u)
        % Solve v = RB'\u for v
        if issparse(B)
            if ~isempty(permB)
                v = RBT \ (permB*u);
            else
                v = RBT \ u;
            end
        else
            RBTopts.LT = true;
            v = linsolve(RBT,u,RBTopts);
        end
    end

%-------------------------------------------------------------------------%
    function v = AminusSigmaBsolve(u)
        % Solve v = (A-sigma*B)\u for v
        if Amatrix
            if ~isempty(dgAsB)
                % use LU reordering permAsB
                 v = qq*(UA \ (LA \ (pp*(dgAsB\u))));
            else
                if ~isempty(precompLU)
                    v = precompLU\u;
                else
                    v = UA \ (LA \ u(pp));
                end;
            end
        else % A is a function
            v = A(u,varargin{afunNargs});
            if isrealprob && ~isreal(v)
                error(message('MATLAB:eigs:complexFunction'));
            end
        end
    end % AminusSigmaBsolve
%-------------------------------------------------------------------------%
    function v = Bsolve(u)
        % Solve v = (A-sigma*B)\u for v
        if ~isempty(dgB)
            % use LU reordering permAsB
            v = qqB*(UB \ (LB \ (ppB*(dgB\u))));
        else
            if ~isempty(precompLUB)
                v = precompLUB\u;
            else
                v = UB \ (LB \ u(ppB));                
            end;
        end
    end % AminusSigmaBsolve

end % mpeigs

function [projr, w] = projectInto(V, r, j)
   w = V(:,1:j)'*r;
   projr = V(:,1:j)*w;
end 

function [r, normRes, stopAlgorithm, w] = robustReorthogonalize(V, r, index, wIn)

    normr0 = norm(r,2);

    if nargin < 4
        wIn = zeros(index,1,'mp');
    end
    w = wIn;

    stopAlgorithm = false;
    
    % Reorthogonalize:
    [projr, dw] = projectInto(V, r, index);
    r = r - projr;
    w = w + dw;
    normRes = norm(r,2);

    numReorths = 1;
    while normRes <= mp('1/sqrt(2)')*normr0 && numReorths < 5
        [projr, dw] = projectInto(V, r, index);
        r = r - projr;
        w = w + dw;
        normr0 = normRes;
        normRes = norm(r,2);
        numReorths = numReorths + 1;
    end

    if normRes <= mp('1/sqrt(2)')*normr0
        % Cannot Reorthogonalize, invariant subspace found.
        normRes = 0;
        w = wIn;

        % Try a random restart
        stopAlgorithm = true;

        for restart=1:3
            % Do a random restart: Will try at most three times
            %r = randn(randStr, size(r,1), 1);
            r = randn(size(r,1), 1);

            % Orthogonalize r
            r = r - projectInto(V, r, index);
            rMr = sqrt(abs(r'*r));
            r = r / rMr;

            % Re-orthogonalize if necessary
            stopAlgorithm = true;
            for reorth=1:5

                % Check orthogonality
                Mr = r;
                [projr, VMr] = projectInto(V, Mr, index);
                rMr = sqrt(abs(r'*Mr));

                %if abs(rMr - 1) <= 1e-10 && all(abs(VMr) <= 1e-10)
                if abs(rMr - 1) <= mp('eps*1e6') && all(abs(VMr) <= mp('eps*1e6'))                
                    stopAlgorithm = false;
                    break;
                end

                % Re-orthogonalize
                r = r - projr;
                r = r / norm(r,2);
            end

            if ~stopAlgorithm
                break;
            end
        end
    else
        r = r/normRes;
    end
end

function [Q, H, stopAlgorithm] = expandKrylov(A, Q, H, sk, ek)
% Expands Krylov subspace.
%
%  The function contruct the sk+1, sk+2, ..., ek_th column of Q.
%  A * Q(:, 1:ek) = Q(:, 1:ek+1) * H
% Parameters:
%   sk       start index
%   ek       end index
% Return:
%   Q        the expanded orthonormal matrix with dimension [n x ek+1]
%   H        dimension [ek+1 x ek], the upper [ek x ek] block is Hessenberg

% TODO:
%   - (Speed-up) Convert projectInto to C++
%
    if 1
       % Heavily optimized and implemented in C++:
       mgs = modifiedGramSchmidt(Q,H,sk,ek); % initialize internal structures

       stopAlgorithm = false;
       for j = sk:ek
            if isfloat(A), v = A * mgs.lastVector();
            else,          v = A(  mgs.lastVector()); end % apply linear operator to last vector in a basis       
            
            % Orthogonolize new vector and update basis.
            if mgs.addVector(v) == false
                % If orthogonal vector cannot be constructed (invariant subspace was already found) 
                % then terminate algorithm.
                stopAlgorithm = true;   
                return;
            end
       end

       [Q,H] = mgs.factors();      % extract factors
       delete(mgs);                % free resources to avoid memory leaks
       
    else    
        for j = sk : ek
            if isfloat(A)
                r = A * Q(:,j);
            else
                r = A(Q(:,j));
            end
                     
            [projr, w] = projectInto(Q, r, j);
            r = r - projr;
            
            % Reorthogonalize  
            [r, beta, stopAlgorithm, w] = robustReorthogonalize(Q, r, j, w);

            if stopAlgorithm
                return;
            end

            % Save data
            H(1:j, j) = w;
            H(j+1, j) = beta;
            Q(:,j+1)  = r;            
        end
    end
end

function object = numeric_t(expression, class_t)
    narginchk(2,2);
    if (strcmpi(class_t,'mp')), object = mp(expression);
    else
        if isnumeric(expression)
            object = expression;
        else
            object = eval(expression);
        end;
    end;
end  % numeric_t
