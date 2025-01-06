function r = magic(n, classname)
    r = magicImpl(n);
    
    if nargin < 2, classname = 'double'; end    
    if strcmpi(classname,'mp') || mp.OverrideDoubleBasicArrays(), r = mp(r); end
end

function M = magicImpl(n)
%MAGIC  Magic square.
%   MAGIC(N) is an N-by-N matrix constructed from the integers
%   1 through N^2 with equal row, column, and diagonal sums.
%   Produces valid magic squares for all N > 0 except N = 2.

%   Copyright 1984-2011 The MathWorks, Inc. 
%   Copyright 2008-2021 Advanpix LLC.

    n = floor(real(double(n(1))));

    if mod(n,2) == 1
       % Odd order 
       M = oddOrderMagicSquare(n);
    elseif mod(n,4) == 0
       % Doubly even order.
       % Doubly even order.
       J = fix(mod(1:n,4)/2);
       K = bsxfun(@eq,J',J);
       M = bsxfun(@plus,(1:n:(n*n))',0:n-1);
       M(K) = n*n+1 - M(K);
    else
       % Singly even order.
       p = n/2;   %p is odd.
       M = oddOrderMagicSquare(p);
       M = [M M+2*p^2; M+3*p^2 M+p^2];
       if n == 2
          return
       end
       i = (1:p)';
       k = (n-2)/4;
       j = [1:k (n-k+2):n];
       M([i; i+p],j) = M([i+p; i],j);
       i = k+1;
       j = [1 i];
       M([i; i+p],j) = M([i+p; i],j);
    end
end

function M = oddOrderMagicSquare(n)
p = 1:n;
M = n*mod(bsxfun(@plus,p',p-(n+3)/2),n) + mod(bsxfun(@plus,p',2*p-2),n) + 1;
end