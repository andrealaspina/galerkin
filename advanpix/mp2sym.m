function r = mp2sym(s)
%MP2SYM Converts multiprecision floating point object to sym object
%
%    Example:
%
%         A = mp(rand(2,2));
%         B = mp2sym(A)
% 
%         B =
%         [ 0.65547789017755664353614974970696,  0.70604608801960877517700509997667]
%         [ 0.17118668781156176628144294227241, 0.031832846377420676020619794144295]
%
%    See also mp.Digits

%    Copyright (c) 2006 - 2022 Advanpix LLC. 
    
    r = sym(zeros(size(s)));
    for n=1:numel(s), r(n) = sym(strrep(num2str(s(n)),'i','*i')); end;
end
