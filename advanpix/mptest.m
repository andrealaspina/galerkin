function mptest()
% MPTEST - unit tests for Multiprecision Computing Toolbox
% Copyright (c) 2008 - 2022 Advanpix LLC. 
% http://advanpix.com  

    % Setup working precision to 34 decimal digits (quadruple precision)
    mp.Digits(34);
    mp.GuardDigits(0);

    % In most cases we test compatibility with MATLAB's built-in functions.
    % Since MATLAB doesn't support high precision we have to level with its capabilities 
    % which is around 7-10 (max 15 ) correct digits, so we use it as tolerance.
    tolerance = mp('1e-7');

    % In other cases when MCT introduces new functions (for which MATLAB doesn't have counterparts)
    % we test their correctness by mathematical properties, special values and relations with already tested functions. 
    
    % Setup size of test matrices
    Rows = 5;
    Cols = 10;
    
    % Repetitions where appropriate
    N = 10;
    
    % Every test is separated by line of '*' symbols
    
    % **************************************************************************
    % Digits()     
    fprintf('%-20s:','mp.Digits()')    
    
    digits = mp.Digits();

    % Default precision control
    mp.Digits(500);
    pi = mp('pi');
    d = mp.Digits(pi)-500;    
        if d ~= 0, fprintf(' 0 <- fail \n'), end;
     
    % Per-number precision control
    pi = mp('pi',1000);
    d  = mp.Digits(pi)-1000;    
        if d ~= 0, fprintf(' 1 <- fail \n'), end;
    
    mp.Digits(digits);
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % GuardDigits()     
    fprintf('%-20s:','mp.GuardDigits()')    
    
    digits = mp.GuardDigits();

    mp.GuardDigits(50);
    d = mp.GuardDigits()-50;    
        if d ~= 0, fprintf(' 0 <- fail \n'), end;

    mp.GuardDigits(digits);
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % double()
    fprintf('%-20s:','double()')    

    X = rand(Rows,Cols)-0.5;

    d = X-double(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = X-double(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % int8()
    fprintf('%-20s:','int8()')    

    X = double(intmax('int8')) * (rand(Rows,Cols)-0.5);

    d = int8(X)-int8(mp(X));    
        if norm(double(d),1)>tolerance, fprintf(' 0 <- fail \n'), end;
  
    fprintf('\t<- success \n')

    % **************************************************************************
    % uint8()
    fprintf('%-20s:','uint8()')    

    X = double(intmax('uint8')) * (rand(Rows,Cols)-0.5);

    d = uint8(X)-uint8(mp(X));    
        if norm(double(d),1)>tolerance, fprintf(' 0 <- fail \n'), end;
  
    fprintf('\t<- success \n')

    % **************************************************************************
    % int16()
    fprintf('%-20s:','int16()')    

    X = double(intmax('int16')) * (rand(Rows,Cols)-0.5);

    d = int16(X)-int16(mp(X));    
        if norm(double(d),1)>tolerance, fprintf(' 0 <- fail \n'), end;
  
    fprintf('\t<- success \n')

    % **************************************************************************
    % uint16()
    fprintf('%-20s:','uint16()')    

    X = double(intmax('uint16')) * (rand(Rows,Cols)-0.5);

    d = uint16(X)-uint16(mp(X));    
        if norm(double(d),1)>tolerance, fprintf(' 0 <- fail \n'), end;
  
    fprintf('\t<- success \n')

    % **************************************************************************
    % int32()
    fprintf('%-20s:','int32()')    

    X = double(intmax('int32')) * (rand(Rows,Cols)-0.5);

    d = int32(X)-int32(mp(X));    
        if norm(double(d),1)>tolerance, fprintf(' 0 <- fail \n'), end;
  
    fprintf('\t<- success \n')

    % **************************************************************************
    % uint32()
    fprintf('%-20s:','uint32()')    

    X = double(intmax('uint32')) * (rand(Rows,Cols)-0.5);

    d = uint32(X)-uint32(mp(X));    
        if norm(double(d),1)>tolerance, fprintf(' 0 <- fail \n'), end;
  
    fprintf('\t<- success \n')

    % **************************************************************************
    % int64()
    fprintf('%-20s:','int64()')    

    X = double(intmax('int64')) * (rand(Rows,Cols)-0.5);

    d = double(int64(X))-double(int64(mp(X)));    
        if norm(double(d),1)>tolerance, fprintf(' 0 <- fail \n'), end;
  
    fprintf('\t<- success \n')

    % **************************************************************************
    % uint64()
    fprintf('%-20s:','uint64()')    

    X = double(intmax('uint64')) * (rand(Rows,Cols)-0.5);

    d = double(uint64(X))-double(uint64(mp(X)));    
        if norm(double(d),1)>tolerance, fprintf(' 0 <- fail \n'), end;
  
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % colon()
    fprintf('%-20s:','colon()')    

    X = mp('-1.01'):mp('1e-3'):mp('1.01');
    Y = -1.01 : 1e-3 : 1.01;
    
    d = Y-double(X);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = Y-X;    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

    X = mp('-100'):mp('100');
    Y = -100 : 100;
    
    d = Y-double(X);    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    d = Y-X;    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;

    X = mp('-100'):mp('35'):mp('100');
    Y = -100 :35: 100;
    
    d = Y-double(X);    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;

    d = Y-X;    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    %% **************************************************************************    
    fprintf('%-20s:','plus()')        
    
    %% plus() - real:

    % scalar - scalar
    for n=1:N
        X = rand(1)-0.5;
        Y = rand(1)-0.5;
        
        d = (X+Y)-(mp(X)+mp(Y));    
            if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

        d = (X+Y)-double(mp(X)+mp(Y));    
            if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    end;
    
    % matrix - scalar
    X = rand(Rows,Cols)-0.5;
    Y = rand(1)-0.5;

    d = (X+Y)-(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    d = (X+Y)-double(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
        
    % scalar - matrix
    X = rand(1)-0.5;
    Y = rand(Rows,Cols)-0.5;

    d = (X+Y)-(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;

    d = (X+Y)-double(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;
    
    % matrix - matrix
    X = rand(Rows,Cols)-0.5;
    Y = rand(Rows,Cols)-0.5;

    d = (X+Y)-(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;

    d = (X+Y)-double(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;
        
    %% plus() - complex:
    
    % scalar - scalar
    for n=1:N
        X = rand(1)-0.5+(rand(1)-0.5)*1i;
        Y = rand(1)-0.5+(rand(1)-0.5)*1i;
        
        d = (X+Y)-(mp(X)+mp(Y));    
            if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), end;

        d = (X+Y)-double(mp(X)+mp(Y));    
            if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), end;
    end;
    
    % matrix - scalar
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(1)-0.5+(rand(1)-0.5)*1i;

    d = (X+Y)-(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;

    d = (X+Y)-double(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;
        
    % scalar - matrix
    X = rand(1)-0.5+(rand(1)-0.5)*1i;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = (X+Y)-(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;

    d = (X+Y)-double(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;
    
    % matrix - matrix
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = (X+Y)-(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), end;

    d = (X+Y)-double(mp(X)+mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    fprintf('%-20s:','minus()')    
    
    %% minus() - real
    
    % scalar - scalar
    for n=1:N
        X = rand(1)-0.5;
        Y = rand(1)-0.5;
        
        d = (X-Y)-(mp(X)-mp(Y));    
            if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

        d = (X-Y)-double(mp(X)-mp(Y));    
            if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    end;
    
    % matrix - scalar
    X = rand(Rows,Cols)-0.5;
    Y = rand(1)-0.5;

    d = (X-Y)-(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    d = (X-Y)-double(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
        
    % scalar - matrix
    X = rand(1)-0.5;
    Y = rand(Rows,Cols)-0.5;

    d = (X-Y)-(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;

    d = (X-Y)-double(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;
    
    % matrix - matrix
    X = rand(Rows,Cols)-0.5;
    Y = rand(Rows,Cols)-0.5;

    d = (X-Y)-(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;

    d = (X-Y)-double(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;
        
    %% minus() - complex:
    
    % scalar - scalar
    for n=1:N
        X = rand(1)-0.5+(rand(1)-0.5)*1i;
        Y = rand(1)-0.5+(rand(1)-0.5)*1i;
        
        d = (X-Y)-(mp(X)-mp(Y));    
            if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), end;

        d = (X-Y)-double(mp(X)-mp(Y));    
            if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), end;
    end;
    
    % matrix - scalar
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(1)-0.5+(rand(1)-0.5)*1i;

    d = (X-Y)-(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;

    d = (X-Y)-double(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;
        
    % scalar - matrix
    X = rand(1)-0.5+(rand(1)-0.5)*1i;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = (X-Y)-(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;

    d = (X-Y)-double(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;
    
    % matrix - matrix
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = (X-Y)-(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), end;

    d = (X-Y)-double(mp(X)-mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % times()
    fprintf('%-20s:','times()')    

    X = rand(Rows,Cols)-rand(Rows,Cols);
    Y = rand(Rows,Cols)-rand(Rows,Cols);

    d = (X.*Y)-(mp(X).*mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = (X.*Y)-double(mp(X).*mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = (X.*Y)-(mp(X).*mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    d = (X.*Y)-double(mp(X).*mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % mtimes()
    fprintf('%-20s:','mtimes()')    
    
    X = rand(Rows,Cols);
    Y = rand(Cols,Rows);

    d = (X*Y)-(mp(X)*mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(Cols,Rows)-0.5+(rand(Cols,Rows)-0.5)*1i;

    d = (X*Y)-(mp(X)*mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

    d = (X*Y)-double(mp(X)*mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % rdivide()
    fprintf('%-20s:','rdivide()')    

    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = (X./Y)-(mp(X)./mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = (X./Y)-double(mp(X)./mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % ldivide()
    fprintf('%-20s:','ldivide()')    

    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = (X.\Y)-(mp(X).\mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = (X.\Y)-double(mp(X).\mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % mldivide() - solution of Ax = b - Dense cases
    fprintf('%-20s:','mldivide()')    

    % General square problem
    A = mp(rand(Rows));
    b = mp(rand(Rows,1));
    
    [L,U,P] = lu(A);
    
    x = U\(L\(P*b));
    y = A\b;
    
    d = x - y;
         if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
    
    error = norm(A*x-b,1) / (norm(A,1) * norm(b, 1));
         if error > tolerance, fprintf(' 1 <- fail \n'), end;
         
    A = mp(rand(Rows,Rows))+ mp(rand(Rows,Rows)*1i);
    b = mp(rand(Rows,1))   + mp(rand(Rows,1)*1i);
    
    [L,U,P] = lu(A);
   
    x = U\(L\(P*b));
    y = A\b;
    
    d = x - y;
         if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
         
    error = norm(A*x-b,1) / (norm(A,1) * norm(b, 1));
         if error > tolerance, fprintf(' 3 <- fail \n'), end;
    
    % Hermitian positive definite problem
    A = mp(rand(Rows));
    b = mp(rand(Rows,1));
    A = A*A';
    
    R = chol(A);
    
    x = R\(R'\b);
    y = A\b;
    
    d = x - y;
         if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;
         
    error = norm(A*x-b,1) / (norm(A,1) * norm(b, 1));
         if error > tolerance, fprintf(' 5 <- fail \n'), end;
         
    A = mp(rand(Rows,Rows))+ mp(rand(Rows,Rows)*1i);
    b = mp(rand(Rows,1))   + mp(rand(Rows,1)*1i);
    A = A*A';
    
    R = chol(A);
    
    x = R\(R'\b);
    y = A\b;
    
    d = x - y;
         if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;
         
    error = norm(A*x-b,1) / (norm(A,1) * norm(b, 1));
         if error > tolerance, fprintf(' 7 <- fail \n'), end;
    
    % under-determined system
    A = mp(rand(Rows, 2*Rows));
    b = mp(rand(Rows, 3*Rows));
    
    x = A\b;
    error = norm(A*x-b,1) / (norm(A,1) * norm(b, 1));
         if error > tolerance, fprintf(' 8 <- fail \n'), end;
        
    % over-determined system        
    A = mp(rand(2*Rows, Rows));
    b = mp(rand(2*Rows, 3*Rows));

    x = A\b;
    y = pinv(A)*b; 
    d = x-y;       % compare with SVD-based least squares solution
        if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), end;

    % under-determined system
    A = mp(rand(Rows, 2*Rows)*1i);
    b = mp(rand(Rows, 3*Rows)*1i);
    
    x = A\b;
    error = norm(A*x-b,1) / (norm(A,1) * norm(b, 1));
         if error > tolerance, fprintf(' 10 <- fail \n'), end;
        
    % over-determined system        
    A = mp(rand(2*Rows, Rows)*1i);
    b = mp(rand(2*Rows, 3*Rows)*1i);

    x = A\b;
    y = pinv(A)*b;
    d = x-y;       % compare with SVD-based least squares solution
        if norm(d,1)>tolerance, fprintf(' 11 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % mrdivide()
    fprintf('%-20s:','mrdivide()')    

    A = mp(magic(Rows));
    B = A/2;

    d = (B/A)-(A'\B')';    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    A = mp(magic(Rows)+magic(Rows)*i);
    B = A/2;

    d = (B/A)-(A'\B')';    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % mpower()
    fprintf('%-20s:','mpower()')    

    X = rand();
    Y = rand();

    d = (mp(X)^mp(Y))-(X^Y);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    X = rand()-0.5+(rand()-0.5)*1i;
    Y = rand();

    d = (mp(X)^mp(Y))-(X^Y);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    X = rand();
    Y = rand()-0.5+(rand()-0.5)*1i;

    d = (mp(X)^mp(Y))-(X^Y);    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    X = rand()-0.5+(rand()-0.5)*1i;
    Y = rand()-0.5+(rand()-0.5)*1i;

    d = (mp(X)^mp(Y))-(X^Y);    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
    
    for n=0:10
        
        % general real matrix        
        X = mp(rand(Rows,Rows)-0.5);

        d = double(X)^n - X^n;
            if norm(d,1) > tolerance, fprintf(' 4 <- fail \n'), end;

        % general complex matrix   
        X = mp(rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i);

        d = double(X)^n - X^n;
            if norm(d,1) > tolerance, fprintf(' 5 <- fail \n'), end;
    end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % power()
    fprintf('%-20s:','power()')    

    N = Rows;
    X = rand(Rows,Cols)-0.5;
    Y = rand(Rows,Cols)-0.5;

    d = (mp(X).^mp(Y))-(X.^Y);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = (mp(X).^mp(Y))-(X.^Y);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    X = rand(Rows,Cols)-0.5;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = (mp(X).^mp(Y))-(X.^Y);    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = (mp(X).^mp(Y))-(X.^Y);    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % power()
    fprintf('%-20s:','realpow()')    

    X = rand();
    Y = rand();

    d = realpow(mp(X),mp(Y))-realpow(X,Y);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
        
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % transpose()
    fprintf('%-20s:','transpose()')    

    X = rand(Rows,Cols);

    d = mp(X).'-X.';    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    X = rand(Rows,Cols)+rand(Rows,Cols)*i;

    d = mp(X).'-X.';    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % ctranspose()
    fprintf('%-20s:','ctranspose()')    

    X = rand(Rows,Cols);

    d = mp(X)'-X';    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    X = rand(Rows,Cols)+rand(Rows,Cols)*i;

    d = mp(X)'-X';    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % uminus()     
    fprintf('%-20s:','uminus()')    

    X = rand(Rows,Cols)-0.5;
    
    d = uminus(X)-uminus(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = uminus(X)-uminus(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % uplus()     
    fprintf('%-20s:','uplus()')    

    X = rand(Rows,Cols)-0.5;
    
    d = uplus(X)-uplus(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = uplus(X)-uplus(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % sin()     
    fprintf('%-20s:','sin()')    

    X = rand(Rows,Cols)-0.5;
    
    d = sin(X)-sin(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = sin(X)-sin(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % cos()     
    fprintf('%-20s:','cos()')    

    X = rand(Rows,Cols)-0.5;
    
    d = cos(X)-cos(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = cos(X)-cos(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % tan()     
    fprintf('%-20s:','tan()')    

    X = rand(Rows,Cols)-0.5;
    
    d = tan(X)-tan(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = tan(X)-tan(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % sec()     
    fprintf('%-20s:','sec()')    

    X = rand(Rows,Cols)-0.5;
    
    d = sec(X)-sec(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = sec(X)-sec(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % csc()     
    fprintf('%-20s:','csc()')    

    X = rand(Rows,Cols)-0.5;
    
    d = csc(X)-csc(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = csc(X)-csc(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % cot()     
    fprintf('%-20s:','cot()')    

    X = rand(Rows,Cols)-0.5;
    
    d = cot(X)-cot(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = cot(X)-cot(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % acos()     
    fprintf('%-20s:','acos()')    

    X = rand(Rows,Cols)-0.5;
    
    d = acos(X)-acos(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = acos(X)-acos(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % asin()     
    fprintf('%-20s:','asin()')    

    X = rand(Rows,Cols)-0.5;
    
    d = asin(X)-asin(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = asin(X)-asin(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % atan()     
    fprintf('%-20s:','atan()')    

    X = rand(Rows,Cols)-0.5;
    
    d = atan(X)-atan(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = atan(X)-atan(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % acot()     
    fprintf('%-20s:','acot()')    

    X = rand(Rows,Cols)-0.5;
    
    d = acot(X)-acot(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = acot(X)-acot(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % atan2()
    fprintf('%-20s:','atan2()')    

    X = rand(Rows,Cols)-0.5;
    Y = rand(Rows,Cols)-0.5;

    d = atan2(mp(Y),mp(X))-atan2(Y,X);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % hypot()
    fprintf('%-20s:','hypot()')    

    X = rand(Rows,Cols)-0.5;
    Y = rand(Rows,Cols)-0.5;

    d = hypot(mp(Y),mp(X))-hypot(Y,X);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(Rows,Cols)-0.5;

    d = hypot(mp(Y),mp(X))-hypot(Y,X);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

    X = rand(Rows,Cols)-0.5;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = hypot(mp(Y),mp(X))-hypot(Y,X);    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    Y = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;

    d = hypot(mp(Y),mp(X))-hypot(Y,X);    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % cosh()     
    fprintf('%-20s:','cosh()')    

    X = rand(Rows,Cols)-0.5;
    
    d = cosh(X)-cosh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = cosh(X)-cosh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % sinh()     
    fprintf('%-20s:','sinh()')    

    X = rand(Rows,Cols)-0.5;
    
    d = sinh(X)-sinh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = sinh(X)-sinh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % tanh()     
    fprintf('%-20s:','tanh()')    

    X = rand(Rows,Cols)-0.5;
    
    d = tanh(X)-tanh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = tanh(X)-tanh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % sech()     
    fprintf('%-20s:','sech()')    

    X = rand(Rows,Cols)-0.5;
    
    d = sech(X)-sech(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = sech(X)-sech(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % csch()     
    fprintf('%-20s:','csch()')    

    X = rand(Rows,Cols)-0.5;
    
    d = csch(X)-csch(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = csch(X)-csch(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % coth()     
    fprintf('%-20s:','coth()')    

    X = rand(Rows,Cols)-0.5;
    
    d = coth(X)-coth(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = coth(X)-coth(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % acosh()     
    fprintf('%-20s:','acosh()')    

    X = rand(Rows,Cols)-0.5;
    
    d = acosh(X)-acosh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = acosh(X)-acosh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % asinh()     
    fprintf('%-20s:','asinh()')    

    X = rand(Rows,Cols)-0.5;
    
    d = asinh(X)-asinh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = asinh(X)-asinh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % atanh()     
    fprintf('%-20s:','atanh()')    

    X = rand(Rows,Cols)-0.5;
    
    d = atanh(X)-atanh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = atanh(X)-atanh(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % acoth()
    fprintf('%-20s:','acoth()')    

    X = rand(Rows,Cols) + 1;
    
    d = acoth(X)-acoth(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = acoth(X)-acoth(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % asech()     
    fprintf('%-20s:','asech()')    

    X = rand(Rows,Cols);
    
    d = asech(X)-asech(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = asech(X)-asech(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % acsch()     
    fprintf('%-20s:','acsch()')    

    X = rand(Rows,Cols)-0.5;
    
    d = acsch(X)-acsch(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = acsch(X)-acsch(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % exp()     
    fprintf('%-20s:','exp()')    

    X = rand(Rows,Cols)-0.5;
    
    d = exp(X)-exp(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = exp(X)-exp(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % expm1()     
    fprintf('%-20s:','expm1()')    

    X = rand(Rows,Cols)-0.5;
    
    d = expm1(X)-expm1(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = expm1(X)-expm1(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % log()     
    fprintf('%-20s:','log()')    

    X = rand(Rows,Cols)-0.5;
    
    d = log(X)-log(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = log(X)-log(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % log10()     
    fprintf('%-20s:','log10()')    

    X = rand(Rows,Cols)-0.5;
    
    d = log10(X)-log10(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = log10(X)-log10(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % log1p()     
    fprintf('%-20s:','log1p()')    

    X = rand(Rows,Cols)-0.5;
    
    d = log1p(X)-log1p(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = log1p(X)-log1p(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % log2()     
    fprintf('%-20s:','log2()')    

    X = rand(Rows,Cols)-0.5;
    
    d = log2(X)-log2(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = log2(X)-log2(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % nextpow2()
    % MATLAB 2008b has bug in nextpow2, so use another approach.
    fprintf('%-20s:','nextpow2()')    

    X = rand(Rows,Cols)-0.5;
    
    d = ceil(log2(abs(X)))-nextpow2(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = ceil(log2(abs(X)))-nextpow2(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % pow2()     
    fprintf('%-20s:','pow2()')    

    X = rand(Rows,Cols)-0.5;
    
    d = pow2(X)-pow2(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = pow2(X)-pow2(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % sqrt()     
    fprintf('%-20s:','sqrt()')    

    X = rand(Rows,Cols)-0.5;
    
    d = sqrt(X)-sqrt(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    
    d = sqrt(X)-sqrt(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % reallog()     
    fprintf('%-20s:','reallog()')    

    X = rand(Rows,Cols);
    
    d = reallog(X)-reallog(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % realsqrt()     
    fprintf('%-20s:','realsqrt()')    

    X = rand(Rows,Cols);
    
    d = realsqrt(X)-realsqrt(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % nthroot()     
    fprintf('%-20s:','nthroot()')    

    X = rand(Rows,Cols);
    Y = rand(Rows,Cols);    
    
    d = nthroot(X,Y)-nthroot(mp(X),mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % pow2(F,E)     
    fprintf('%-20s:','pow2(F,E)')    

    F = rand(Rows,Cols);
    E = rand(Rows,Cols);    
    
    d = pow2(F,E)-pow2(mp(F),mp(E));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % min(),max()     
    fprintf('%-20s:','min(), max()')    

    x = rand(Rows,Cols)-rand(Rows,Cols);
    y = rand(Rows,Cols)-rand(Rows,Cols);

    X = mp(x);
    Y = mp(y);

    [c,J] = min(x);
    [C,I] = min(X);    
    if norm(C - c,1) > tolerance, fprintf(' 0 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 1 <- fail \n'), end;

    [c,J] = min(x,[],1);
    [C,I] = min(X,[],1);    
    if norm(C - c,1) > tolerance, fprintf(' 2 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 3 <- fail \n'), end;

    [c,J] = min(x,[],2);
    [C,I] = min(X,[],2);    
    if norm(C - c,1) > tolerance, fprintf(' 4 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 5 <- fail \n'), end;
   
    [c,J] = max(x);
    [C,I] = max(X);    
    if norm(C - c,1) > tolerance, fprintf(' 6 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 7 <- fail \n'), end;
        
    [c,J] = max(x,[],1);
    [C,I] = max(X,[],1);    
    if norm(C - c,1) > tolerance, fprintf(' 8 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 9 <- fail \n'), end;

    [c,J] = max(x,[],2);
    [C,I] = max(X,[],2);    
    if norm(C - c,1) > tolerance, fprintf(' 10 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 11 <- fail \n'), end;

    d = min(x,y)-min(X,Y);    
        if norm(d,1) > tolerance, fprintf(' 12 <- fail \n'), end;

    d = max(x,y)-max(X,Y);    
        if norm(d,1) > tolerance, fprintf(' 13 <- fail \n'), end;
    
    x = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    y = rand(Rows,Cols)+rand(Rows,Cols)*1i;

    X = mp(x);
    Y = mp(y);

    [c,J] = min(x);
    [C,I] = min(X);    
    if norm(C - c,1) > tolerance, fprintf(' 14 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 15 <- fail \n'), end;

    [c,J] = min(x,[],1);
    [C,I] = min(X,[],1);    
    if norm(C - c,1) > tolerance, fprintf(' 16 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 17 <- fail \n'), end;

    [c,J] = min(x,[],2);
    [C,I] = min(X,[],2);    
    if norm(C - c,1) > tolerance, fprintf(' 18 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 19 <- fail \n'), end;
   
    [c,J] = max(x);
    [C,I] = max(X);    
    if norm(C - c,1) > tolerance, fprintf(' 20 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 21 <- fail \n'), end;
        
    [c,J] = max(x,[],1);
    [C,I] = max(X,[],1);    
    if norm(C - c,1) > tolerance, fprintf(' 22 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 23 <- fail \n'), end;

    [c,J] = max(x,[],2);
    [C,I] = max(X,[],2);    
    if norm(C - c,1) > tolerance, fprintf(' 24 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 25 <- fail \n'), end;

    d = min(x,y)-min(X,Y);    
        if norm(d,1) > tolerance, fprintf(' 26 <- fail \n'), end;

    d = max(x,y)-max(X,Y);    
        if norm(d,1) > tolerance, fprintf(' 27 <- fail \n'), end;

    x = rand(Rows,Cols);
    y = rand(Rows,Cols)+rand(Rows,Cols)*1i;

    X = mp(x);
    Y = mp(y);

    [c,J] = min(x);
    [C,I] = min(X);    
    if norm(C - c,1) > tolerance, fprintf(' 28 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 29 <- fail \n'), end;

    [c,J] = min(x,[],1);
    [C,I] = min(X,[],1);    
    if norm(C - c,1) > tolerance, fprintf(' 30 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 31 <- fail \n'), end;

    [c,J] = min(x,[],2);
    [C,I] = min(X,[],2);    
    if norm(C - c,1) > tolerance, fprintf(' 32 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 33 <- fail \n'), end;
   
    [c,J] = max(x);
    [C,I] = max(X);    
    if norm(C - c,1) > tolerance, fprintf(' 34 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 35 <- fail \n'), end;
        
    [c,J] = max(x,[],1);
    [C,I] = max(X,[],1);    
    if norm(C - c,1) > tolerance, fprintf(' 36 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 37 <- fail \n'), end;

    [c,J] = max(x,[],2);
    [C,I] = max(X,[],2);    
    if norm(C - c,1) > tolerance, fprintf(' 38 <- fail \n'), end;
    if norm(I - J,1) > tolerance, fprintf(' 39 <- fail \n'), end;

    d = min(x,y)-min(X,Y);    
        if norm(d,1) > tolerance, fprintf(' 40 <- fail \n'), end;

    d = max(x,y)-max(X,Y);    
        if norm(d,1) > tolerance, fprintf(' 41 <- fail \n'), end;

    d = numel(max(X))-numel(max(x));    
        if norm(d,1)>tolerance, fprintf(' 43 <- fail \n'), end;

	d = numel(min(X))-numel(min(x));    
        if norm(d,1)>tolerance, fprintf(' 44 <- fail \n'), end;
		
	x = rand(1,Cols);
	X = mp(x);
	
    d = numel(max(X))-numel(max(x));    
        if norm(d,1)>tolerance, fprintf(' 45 <- fail \n'), end;
  
    d = numel(min(X))-numel(min(x));    
        if norm(d,1)>tolerance, fprintf(' 46 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % prod(matrix)     
    fprintf('%-20s:','prod(matrix)')    

    X = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    
    d = prod(X)-prod(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    d = prod(X,2)-prod(mp(X),2);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % prod(vector)     
    fprintf('%-20s:','prod(vector)')    

    X = rand(1,Cols);
    d = prod(X)-prod(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,1);
    d = prod(X)-prod(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % sum(matrix)     
    fprintf('%-20s:','sum(matrix)')    

    X = rand(Rows,Cols);
    
    d = sum(X)-sum(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    d = sum(X,2)-sum(mp(X),2);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % sum(vector)     
    fprintf('%-20s:','sum(vector)')    

    X = rand(1,Cols);
    d = sum(X)-sum(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,1);
    d = sum(X)-sum(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % cumsum(matrix)     
    fprintf('%-20s:','cumsum(matrix)')    

    X = rand(Rows,Cols);
    
    d = cumsum(X)-cumsum(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    d = cumsum(X,2)-cumsum(mp(X),2);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

	X = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    
    d = cumsum(X)-cumsum(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
   
    d = cumsum(X,2)-cumsum(mp(X),2);    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % cumsum(vector)     
    fprintf('%-20s:','cumsum(vector)')    

    X = rand(1,Cols);
    d = cumsum(X)-cumsum(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,1);
    d = cumsum(X)-cumsum(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

	X = rand(1,Cols)+rand(1,Cols)*1i;
    d = cumsum(X)-cumsum(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
   
    X = rand(Rows,1)+rand(Rows,1)*1i;
    d = cumsum(X)-cumsum(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % cumprod(matrix)     
    fprintf('%-20s:','cumprod(matrix)')    

    X = rand(Rows,Cols);
    
    d = cumprod(X)-cumprod(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    d = cumprod(X,2)-cumprod(mp(X),2);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

	X = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    
    d = cumprod(X)-cumprod(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
   
    d = cumprod(X,2)-cumprod(mp(X),2);    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % cumprod(vector)     
    fprintf('%-20s:','cumprod(vector)')    

    X = rand(1,Cols);
    d = cumprod(X)-cumprod(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,1);
    d = cumprod(X)-cumprod(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

	X = rand(1,Cols)+rand(1,Cols)*1i;
    d = cumprod(X)-cumprod(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
   
    X = rand(Rows,1)+rand(Rows,1)*1i;
    d = cumprod(X)-cumprod(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
    
    fprintf('\t<- success \n')
	
	% **************************************************************************
    % dot(matrix)     
    fprintf('%-20s:','dot()')    

    X = rand(Rows,Cols);
	Y = rand(Rows,Cols);
    
    d = dot(X,Y)-dot(mp(X),mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
   
    d = dot(X,Y,2)-dot(mp(X),mp(Y),2);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

	X = rand(Rows,Cols)+rand(Rows,Cols)*1i;
	Y = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    
    d = dot(X,Y)-dot(mp(X),mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
   
    d = dot(X,Y,2)-dot(mp(X),mp(Y),2);    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;

    % **************************************************************************
    % dot(vector)     
    X = rand(1,Cols);
    Y = rand(1,Cols);	
	
    d = dot(X,Y)-dot(mp(X),mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;
   
    d = dot(X,Y,2)-dot(mp(X),mp(Y),2);    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;
   
    X = rand(Rows,1);
    Y = rand(Rows,1);	
	
    d = dot(X,Y)-dot(mp(X),mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;
   
    d = dot(X,Y,2)-dot(mp(X),mp(Y),2);    
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;

	X = rand(1,Cols)+rand(1,Cols)*1i;
	Y = rand(1,Cols)+rand(1,Cols)*1i;
	
    d = dot(X,Y)-dot(mp(X),mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), end;
   
    d = dot(X,Y,2)-dot(mp(X),mp(Y),2);    
        if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), end;
   
    X = rand(Rows,1)+rand(Rows,1)*1i;
    Y = rand(Rows,1)+rand(Rows,1)*1i;	
	
    d = dot(X,Y)-dot(mp(X),mp(Y));    
        if norm(d,1)>tolerance, fprintf(' 10 <- fail \n'), end;
   
    d = dot(X,Y,2)-dot(mp(X),mp(Y),2);    
        if norm(d,1)>tolerance, fprintf(' 11 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % svd()     
    fprintf('%-20s:','svd()')    

    X = rand(Rows,Rows)-0.5;
    [U,S,V] = svd(mp(X));
    d = X - U*S*V';
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
  
    X = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;
    [U,S,V] = svd(mp(X));
    d = X - U*S*V';
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

    X = rand(Rows,2*Rows)-0.5;
    [U,S,V] = svd(mp(X));
    d = X - U*S*V';
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
  
    X = rand(Rows,2*Rows)-0.5+(rand(Rows,2*Rows)-0.5)*1i;
    [U,S,V] = svd(mp(X));
    d = X - U*S*V';
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;

    X = rand(2*Rows,Rows)-0.5;
    [U,S,V] = svd(mp(X));
    d = X - U*S*V';
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;
  
    X = rand(2*Rows,Rows)-0.5+(rand(2*Rows,Rows)-0.5)*1i;
    [U,S,V] = svd(mp(X));
    d = X - U*S*V';
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;
        
    X = rand(Rows,2*Rows)-0.5;
    [U,S,V] = svd(mp(X),'econ');
    d = X - U*S*V';
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;
  
    X = rand(Rows,2*Rows)-0.5+(rand(Rows,2*Rows)-0.5)*1i;
    [U,S,V] = svd(mp(X),'econ');
    d = X - U*S*V';
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;

    X = rand(2*Rows,Rows)-0.5;
    [U,S,V] = svd(mp(X),'econ');
    d = X - U*S*V';
        if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), end;
  
    X = rand(2*Rows,Rows)-0.5+(rand(2*Rows,Rows)-0.5)*1i;
    [U,S,V] = svd(mp(X),'econ');
    d = X - U*S*V';
        if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % qr()     
    fprintf('%-20s:','qr()')    

    % Real & Full decomposition
    X = rand(Rows,Cols);

    [Q,R] = qr(mp(X));

    d = X - Q*R;
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = Q'*Q - mp(eye(Rows));
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

    d = Q'*X - R;
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
        
    [Q,R,E] = qr(mp(X));

    d = X*E - Q*R;
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;

    % Real & Economic decomposition
    X = rand(2*Rows,Cols);

    [Q,R] = qr(mp(X),0);

    d = X - Q*R;
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;

    d = Q'*Q - mp(eye(2*Rows));
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;

    d = Q'*X - R;
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;
        
    % Real & Full & Rectangular decomposition
    X = rand(2*Rows,Cols);

    [Q,R] = qr(mp(X));

    d = X - Q*R;
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;

    d = Q'*Q - mp(eye(2*Rows));
        if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), end;

    d = Q'*X - R;
        if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), end;
        
    % Complex & Full decomposition
    X = rand(Rows,Cols)+rand(Rows,Cols)*1i;

    [Q,R] = qr(mp(X));

    d = X - Q*R;
        if norm(d,1)>tolerance, fprintf(' 10 <- fail \n'), end;

    d = Q'*Q - mp(eye(Rows));
        if norm(d,1)>tolerance, fprintf(' 11 <- fail \n'), end;

    d = Q'*X - R;
        if norm(d,1)>tolerance, fprintf(' 12 <- fail \n'), end;
        
    [Q,R,E] = qr(mp(X));

    d = X*E - Q*R;
        if norm(d,1)>tolerance, fprintf(' 13 <- fail \n'), end;
        
    % Complex & Economic decomposition
    X = rand(2*Rows,Cols)+rand(2*Rows,Cols)*1i;

    [Q,R] = qr(mp(X),0);

    d = X - Q*R;
        if norm(d,1)>tolerance, fprintf(' 14 <- fail \n'), end;

    d = Q'*Q - mp(eye(2*Rows));
        if norm(d,1)>tolerance, fprintf(' 15 <- fail \n'), end;

    d = Q'*X - R;
        if norm(d,1)>tolerance, fprintf(' 16 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % lu(square)     
    fprintf('%-20s:','lu(square)')    
    
    X = rand(Rows,Rows);
    
    [L,U] = lu(mp(X));
    d = X-L*U;
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    [L,U,P] = lu(mp(X));
    d = P*X - L*U;
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    X = rand(Rows,Rows)+rand(Rows,Rows)*1i;
   
    [L,U] = lu(mp(X));
    d = X-L*U;
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
        
    [L,U,P] = lu(mp(X));
    d = P*X - L*U;
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
        
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % lu(rect)     
    fprintf('%-20s:','lu(rect)')    

    X = rand(Rows,2*Rows);
    
    [L,U] = lu(mp(X));
    d = X-L*U;
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    [L,U,P] = lu(mp(X));
    d = P*X - L*U;
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

    X = rand(Rows,2*Rows)+rand(Rows,2*Rows)*1i;
   
    [L,U] = lu(mp(X));
    d = X-L*U;
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
        
    [L,U,P] = lu(mp(X));
    d = P*X - L*U;
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
   
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % pinv()     
    fprintf('%-20s:','pinv()')    

    X = rand(Rows,Cols);
    Y = pinv(mp(X));
 
    d = X - X*Y*X;
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = Y - Y*X*Y;
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    d = ctranspose(X*Y) - X*Y;
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    d = ctranspose(Y*X) - Y*X;
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;

    X = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    Y = pinv(mp(X));
 
    d = X - X*Y*X;
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;

    d = Y - Y*X*Y;
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;
        
    d = ctranspose(X*Y) - X*Y;
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;

    d = ctranspose(Y*X) - Y*X;
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;
   
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % null()     
    fprintf('%-20s:','null()')    

    X = rand(Rows,Cols);
    Y = null(mp(X));
 
    d = X*Y;
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = Y'*Y;
        if norm(d-mp(eye(size(d,2))),1)>tolerance, fprintf(' 1 <- fail \n'), end;

    X = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    Y = null(mp(X));
 
    d = X*Y;
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    d = Y'*Y;
        if norm(d-mp(eye(size(d,2))),1)>tolerance, fprintf(' 3 <- fail \n'), end;
   
    fprintf('\t<- success \n')
    
    % **************************************************************************    
    % balance() 
    fprintf('%-20s:','balance()')
    
    X = rand(Rows,Rows);
    idx = randi([1,Rows*Rows],min(3,Rows),1);
    X(idx) = 1000;
    
    [T,B] = balance(mp(X));    
    d = norm(B - T\X*T,1);
        if d>tolerance, fprintf(' 0 <- fail \n'), end;
        
    [T,B] = balance(mp(X),'noperm');    
    d = norm(B - T\X*T,1);
        if d>tolerance, fprintf(' 1 <- fail \n'), end;
    
    [S,P,B] = balance(mp(X));    
    T(:,P) = diag(S);
    B(P,P) = diag(1./S)*X*diag(S);
 
    d = norm(B - T\X*T,1);        
        if d>tolerance, fprintf(' 2 <- fail \n'), end;

    X = rand(Rows,Rows)+rand(Rows,Rows)*1i;
    idx = randi([1,Rows*Rows],min(3,Rows),1);
    X(idx) = 1000;
    
    [T,B] = balance(mp(X));    
    d = norm(B - T\X*T,1);
        if d>tolerance, fprintf(' 3 <- fail \n'), end;

    [T,B] = balance(mp(X),'noperm');    
    d = norm(B - T\X*T,1);
        if d>tolerance, fprintf(' 4 <- fail \n'), end;
        
    [S,P,B] = balance(mp(X));    
    T(:,P) = diag(S);
    B(P,P) = diag(1./S)*X*diag(S);
 
    d = norm(B - T\X*T,1);        
        if d>tolerance, fprintf(' 5 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % eig()     
    fprintf('%-20s:','eig()')    

    % eig(A) - nonsymmetric
    X = mp(rand(Rows,Rows));
    [V,D] = eig(X);
    d = X*V - V*D;
         if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
         
    X = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);
    [V,D] = eig(X);
    d = X*V - V*D;
         if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

    % eig(A) - symmetric/hermitian
    X = mp(rand(Rows,Rows));        
    X = X+X';
    [V,D] = eig(X);
    d = X*V - V*D;
         if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    X = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);
    X = X+X';    
    [V,D] = eig(X);
    d = X*V - V*D;
         if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;

    % eig(A) - symmetric tridiagonal
    X = mp(full(gallery('tridiag',Rows,-1,2,-1)));        
    [V,D] = eig(X);
    d = X*V - V*D;
         if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;
         
    [V,D,W] = eig(X);
    d = X*V - V*D;
         if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;
    d = W'*X - D*W';
         if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;

    X = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);
    X = X+X';    
    [V,D] = eig(X);
    d = X*V - V*D;
         if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;
         
    % eig(A) - both eigenvectors
    X = mp(rand(Rows,Rows));
    [V,D,W] = eig(X);
    d = X*V - V*D;
         if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), end;
    d = W'*X - D*W';
         if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), end;

    X = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);
    [V,D,W] = eig(X);
    d = X*V - V*D;
         if norm(d,1)>tolerance, fprintf(' 10 <- fail \n'), end;
    d = W'*X - D*W';
         if norm(d,1)>tolerance, fprintf(' 11 <- fail \n'), end;

    X = mp(rand(Rows,Rows));        
    X = X+X';
    [V,D,W] = eig(X);
    d = X*V - V*D;
         if norm(d,1)>tolerance, fprintf(' 12 <- fail \n'), end;
    d = W'*X - D*W';
         if norm(d,1)>tolerance, fprintf(' 13 <- fail \n'), end;

    X = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);
    X = X+X';    
    [V,D,W] = eig(X);
    d = X*V - V*D;
         if norm(d,1)>tolerance, fprintf(' 14 <- fail \n'), end;
    d = W'*X - D*W';
         if norm(d,1)>tolerance, fprintf(' 15 <- fail \n'), end;
        
    % eig(A,B) - real nonsymmetric
    X = mp(rand(Rows,Rows));
    Y = mp(rand(Rows,Rows));

    [V,D,W] = eig(X,Y);
    d = X*V - Y*V*D;
         if norm(d,1)>tolerance, fprintf(' 16 <- fail \n'), end;
    d = W'*X - D*W'*Y;
         if norm(d,1)>tolerance, fprintf(' 17 <- fail \n'), end;

    [V,D,W] = eig(X,Y,'chol');
    d = X*V - Y*V*D;
         if norm(d,1)>tolerance, fprintf(' 18 <- fail \n'), end;
    d = W'*X - D*W'*Y;
         if norm(d,1)>tolerance, fprintf(' 19 <- fail \n'), end;

    [V,D,W] = eig(X,Y,'qz');
    d = X*V - Y*V*D;
         if norm(d,1)>tolerance, fprintf(' 20 <- fail \n'), end;
    d = W'*X - D*W'*Y;
         if norm(d,1)>tolerance, fprintf(' 21 <- fail \n'), end;

    % eig(A,B) - symmetric
    A = X+X'; 
    B = Y*Y';
    
    [V,D,W] = eig(A,B);
    d = A*V - B*V*D;
         if norm(d,1) > tolerance, fprintf(' 22 <- fail \n'), end;
    d = W'*A - D*W'*B;
         if norm(d,1)>tolerance, fprintf(' 23 <- fail \n'), end;

    [V,D,W] = eig(A,B,'chol');
    d = A*V - B*V*D;
         if norm(d,1) > tolerance, fprintf(' 24 <- fail \n'), end;
    d = W'*A - D*W'*B;
         if norm(d,1)>tolerance, fprintf(' 25 <- fail \n'), end;

    [V,D,W] = eig(A,B,'qz');
    d = A*V - B*V*D;
         if norm(d,1) > tolerance, fprintf(' 26 <- fail \n'), end;
    d = W'*A - D*W'*B;
         if norm(d,1)>tolerance, fprintf(' 27 <- fail \n'), end;

    % eig(A,B) - complex nonsymmetric             
    X = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);
    Y = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);    
       
    [V,D,W] = eig(X,Y);
    d = X*V - Y*V*D;
         if norm(d,1)>tolerance, fprintf(' 28 <- fail \n'), end;
    d = W'*X - D*W'*Y;
         if norm(d,1)>tolerance, fprintf(' 29 <- fail \n'), end;

    [V,D,W] = eig(X,Y,'chol');
    d = X*V - Y*V*D;
         if norm(d,1)>tolerance, fprintf(' 30 <- fail \n'), end;
    d = W'*X - D*W'*Y;
         if norm(d,1)>tolerance, fprintf(' 31 <- fail \n'), end;
         
    [V,D,W] = eig(X,Y,'qz');
    d = X*V - Y*V*D;
         if norm(d,1)>tolerance, fprintf(' 32 <- fail \n'), end;
    d = W'*X - D*W'*Y;
         if norm(d,1)>tolerance, fprintf(' 33 <- fail \n'), end;

    % eig(A,B) - Hermitian
    A = X+X'; 
    B = Y*Y';

    [V,D,W] = eig(A,B);
    d = A*V - B*V*D;
         if norm(d,1) > tolerance, fprintf(' 34 <- fail \n'), end;
    d = W'*A - D*W'*B;
         if norm(d,1)>tolerance, fprintf(' 35 <- fail \n'), end;
         
    [V,D,W] = eig(A,B,'chol');
    d = A*V - B*V*D;
         if norm(d,1) > tolerance, fprintf(' 36 <- fail \n'), end;
    d = W'*A - D*W'*B;
         if norm(d,1)>tolerance, fprintf(' 37 <- fail \n'), end;

    [V,D,W] = eig(A,B,'qz');
    d = A*V - B*V*D;
         if norm(d,1) > tolerance, fprintf(' 38 <- fail \n'), end;
    d = W'*A - D*W'*B;
         if norm(d,1)>tolerance, fprintf(' 39 <- fail \n'), end;
         
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % qz()     
        
    fprintf('%-20s:','qz()');            

    A = mp(rand(Rows,Rows));        
    B = mp(rand(Rows,Rows));
    s = (norm(A,1) * norm(B,1));
    
    [AA,BB,Q,Z,V,W] = qz(A,B); % complex

    a = diag(AA);
    b = diag(BB);
    l = a./b;

    d = Q*A*Z-AA;
      if norm(d,1)/s > tolerance, fprintf(' 1 <- fail \n'), end;

    d = Q*B*Z-BB;
      if norm(d,1)/s > tolerance, fprintf(' 2 <- fail \n'), end;

    d = A*V - B*V*diag(l);
      if norm(d,1)/s > tolerance, fprintf(' 3 <- fail \n'), end;

    d = W'*A - diag(l)*W'*B;
      if norm(d,1)/s > tolerance, fprintf(' 4 <- fail \n'), end;

    d = A*V*b - B*V*a;
      if norm(d,1)/s > tolerance, fprintf(' 5 <- fail \n'), end;

    [AA,BB,Q,Z,V,W] = qz(A,B,'real');

    d = Q*A*Z-AA;
      if norm(d,1)/s > tolerance, fprintf(' 6 <- fail \n'), end;

    d = Q*B*Z-BB;
      if norm(d,1)/s > tolerance, fprintf(' 7 <- fail \n'), end;
      
    A = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);        
    B = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);
    s = (norm(A,1) * norm(B,1));
    
    [AA,BB,Q,Z,V,W] = qz(A,B);

    a = diag(AA);
    b = diag(BB);
    l = a./b;

    d = Q*A*Z-AA;
      if norm(d,1)/s > tolerance, fprintf(' 8 <- fail \n'), end;

    d = Q*B*Z-BB;
      if norm(d,1)/s > tolerance, fprintf(' 9 <- fail \n'), end;

    d = A*V - B*V*diag(l);
      if norm(d,1)/s > tolerance, fprintf(' 10 <- fail \n'), end;

    d = W'*A - diag(l)*W'*B;
      if norm(d,1)/s > tolerance, fprintf(' 11 <- fail \n'), end;

    d = A*V*b - B*V*a;
      if norm(d,1)/s > tolerance, fprintf(' 12 <- fail \n'), end;

    fprintf('\t<- success \n')
        
    % **************************************************************************
    % hess()     
    fprintf('%-20s:','hess()')    
    
   
    X = mp(rand(Rows,Rows));
    [P,H] = hess(X);

    d = X - P*H*P';
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = P*P' - eye(size(X));
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

    X = mp(rand(Rows,Rows))+mp(rand(Rows,Rows))*1i;
    [P,H] = hess(X);

    d = X - P*H*P';
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    d = P*P' - eye(size(X));
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % chol()     
    fprintf('%-20s:','chol()')    

    X = rand(Rows,Rows);
    A = X'*X;

    R = chol(mp(A));
    L = chol(mp(A),'lower');
    
    d = A - R'*R;
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = A - L*L';
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

    X = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;
    A = X'*X;

    R = chol(mp(A));
    L = chol(mp(A),'lower');
    
    d = A - R'*R;
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    d = A - L*L';
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % schur()     
    fprintf('%-20s:','schur()')    

    A = rand(Rows,Rows);
    s = norm(A,1);
    
    [U,T] = schur(mp(A));
    
    d = A - U*T*U';
        if norm(d,1)/s >tolerance, fprintf(' 0 <- fail \n'), end;

    [U,T] = schur(mp(A),'complex');
    
    d = A - U*T*U';
        if norm(d,1)/s>tolerance, fprintf(' 1 <- fail \n'), end;

    A = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;

    [U,T] = schur(mp(A));
    
    d = A - U*T*U';
        if norm(d,1)/s>tolerance, fprintf(' 2 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % ordschur() 
    % We cannot compare ordschur to MATLAB's one 
    % since we do minimization of permutations
    fprintf('%-20s:','ordschur()')    
    
    A = mp(rand(Rows, Rows));
    [U,T] = schur(A,'complex');
    
    for n = 1:50
        clusters = randi(Rows-1, 1, Rows);
        [US, TS] = ordschur(U, T, clusters);

         d = A - US*TS*US';
            if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
    end

    fprintf('\t<- success \n')
    
    % **************************************************************************
    % rank()     
    fprintf('%-20s:','rank()')    

    X = rand(Rows,Cols);

    d = rank(X) - rank(mp(X));
        if norm(d,1) > tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Cols)+rand(Rows,Cols)*1i;

    d = rank(X) - rank(mp(X));
        if norm(d,1) > tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % trace()     
    fprintf('%-20s:','trace()')    

    X = rand(Rows,Rows);

    d = trace(X) - trace(mp(X));
        if norm(d,1) > tolerance, fprintf(' 0 <- fail \n'), end;
   
    X = rand(Rows,Rows)+rand(Rows,Rows)*1i;

    d = trace(X) - trace(mp(X));
        if norm(d,1) > tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % det()     
    fprintf('%-20s:','det()')    

    % general real matrix
    X = mp(rand(Rows,Rows));

    d = det(double(X)) - det(X);
        if norm(d,1) > tolerance, fprintf(' 0 <- fail \n'), end;
    
    % general complex matrix   
    X = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);

    d = det(double(X)) - det(X);
        if norm(d,1) > tolerance, fprintf(' 1 <- fail \n'), end;
        
    % symmetric real matrix
    X = mp(rand(Rows,Rows));
    X = X+X';    
    
    d = det(double(X)) - det(X);
        if norm(d,1) > tolerance, fprintf(' 2 <- fail \n'), end;
        
    % hermitian complex matrix   
    X = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);
    X = X+X';    
    
    d = det(double(X)) - det(X);
        if norm(d,1) > tolerance, fprintf(' 3 <- fail \n'), end;
        
    % symmetric positive definite real matrix
    X = mp(rand(Rows,Rows));
    X = X*X';    
    
    d = det(double(X)) - det(X);
        if norm(d,1) > tolerance, fprintf(' 4 <- fail \n'), end;
        
    % self-adjoint positive definite complex matrix   
    X = mp(rand(Rows,Rows)+rand(Rows,Rows)*1i);
    X = X*X';    
    
    d = det(double(X)) - det(X);
        if norm(d,1) > tolerance, fprintf(' 5 <- fail \n'), end;
        
    % upper triangular real matrix
    X = mp(triu(rand(Rows,Rows)));
   
    d = det(double(X)) - det(X);
        if norm(d,1) > tolerance, fprintf(' 6 <- fail \n'), end;
        
    % lower triangular real matrix
    X = mp(tril(rand(Rows,Rows)));
   
    d = det(double(X)) - det(X);
        if norm(d,1) > tolerance, fprintf(' 7 <- fail \n'), end;

    % upper triangular complex matrix
    X = mp(triu(rand(Rows,Rows)+rand(Rows,Rows)*1i));
   
    d = det(double(X)) - det(X);
        if norm(d,1) > tolerance, fprintf(' 8 <- fail \n'), end;
        
    % lower triangular complex matrix
    X = mp(tril(rand(Rows,Rows)+rand(Rows,Rows)*1i));
   
    d = det(double(X)) - det(X);
        if norm(d,1) > tolerance, fprintf(' 9 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % inv()     
    fprintf('%-20s:','inv()')    

    X = mp(magic(5));
    Y = inv(X);

    d = X - inv(Y);
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;    

    d = inv(X') - Y';
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;    

    X = mp(magic(5)*1i);
    Y = inv(X);

    d = X - inv(Y);
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;    

    d = inv(X') - Y';
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;    
       
    fprintf('\t<- success \n')

    % **************************************************************************
    % sort()     
    fprintf('%-20s:','sort(real)')    

    A = magic(10);
    X = mp(A);

    [B,IX] = sort(X);
    [b,ix] = sort(A);

    d = B - b;
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;    

    d = IX - ix;
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;    

    [B,IX] = sort(X,1,'descend');
    [b,ix] = sort(A,1,'descend');

    d = B - b;
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;    

    d = IX - ix;
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;    

    [B,IX] = sort(X,2);
    [b,ix] = sort(A,2);

    d = B - b;
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;    

    d = IX - ix;
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;    

    [B,IX] = sort(X,2,'descend');
    [b,ix] = sort(A,2,'descend');
    
    d = B - b;
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;    

    d = IX - ix;
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;    

%%
    A = cat(3, magic(3), magic(3)');
    A = cat(4, A, A, A);    
    X = mp(A);

    [B,IX] = sort(X,3);
    [b,ix] = sort(A,3);

    d = B(:) - b(:);
        if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), end;    

    d = IX(:) - ix(:);
        if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), end;    

    [B,IX] = sort(X,3,'descend');
    [b,ix] = sort(A,3,'descend');

    d = B(:) - b(:);
        if norm(d,1)>tolerance, fprintf(' 10 <- fail \n'), end;    

    d = IX(:) - ix(:);
        if norm(d,1)>tolerance, fprintf(' 11 <- fail \n'), end;    

    [B,IX] = sort(X,4);
    [b,ix] = sort(A,4);

    d = B(:) - b(:);
        if norm(d,1)>tolerance, fprintf(' 12 <- fail \n'), end;    

    d = IX(:) - ix(:);
        if norm(d,1)>tolerance, fprintf(' 13 <- fail \n'), end;    

    [B,IX] = sort(X,4,'descend');
    [b,ix] = sort(A,4,'descend');
    
    d = B(:) - b(:);
        if norm(d,1)>tolerance, fprintf(' 14 <- fail \n'), end;    

    d = IX(:) - ix(:);
        if norm(d,1)>tolerance, fprintf(' 15 <- fail \n'), end;    

%%

    fprintf('\t<- success \n')

    % **************************************************************************
    % sort()     
    fprintf('%-20s:','sort(complex)')    

    A = magic(10)+magic(10)*i;
    X = mp(A);

    [B,IX] = sort(X);
    [b,ix] = sort(A);

    d = B - b;
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;    

    d = IX - ix;
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;    

    [B,IX] = sort(X,1,'descend');
    [b,ix] = sort(A,1,'descend');

    d = B - b;
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;    

    d = IX - ix;
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;    

    [B,IX] = sort(X,2);
    [b,ix] = sort(A,2);

    d = B - b;
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;    

    d = IX - ix;
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;    

    [B,IX] = sort(X,2,'descend');
    [b,ix] = sort(A,2,'descend');
    
    d = B - b;
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;    

    d = IX - ix;
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;    
       
    fprintf('\t<- success \n')

    % **************************************************************************
    % find()     
    fprintf('%-20s:','find()')    

    A = magic(10);
    X = mp(A);

    [R,C,V] = find(X);
    [r,c,v] = find(A);

    d = R - r;
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;    
        
    d = C - c;
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;    
        
    d = V - v;
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;    

    [R,C,V] = find(X,5,'last');
    [r,c,v] = find(A,5,'last');

    d = R - r;
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;    
        
    d = C - c;
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;    
        
    d = V - v;
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;    

    IX = find(X);
    ix = find(A);

    d = IX - ix;
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), end;    

    IX = find(X,5,'last');
    ix = find(A,5,'last');

    d = IX - ix;
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), end;    
       
    fprintf('\t<- success \n')

    % **************************************************************************
    % lt(), le(), gt(), ge()     
    fprintf('%-20s:','<,<=,>,>=,==,~=')    

    X = rand(Rows,Cols);
    Y = rand(Rows,Cols);
    x = mp(X);
    y = mp(Y);
    
    d = (x > y)-(X > Y);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = (x >= y)-(X >= Y);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
   
    d = (x < y)-(X < Y);    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    d = (x <= y)-(X <= Y);    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;

    X(:,1) = Y(:,1);
    x(:,1) = y(:,1);
    
    d = (x == y)-(X == Y);    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;

    d = (x ~= y)-(X ~= Y);    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;

    X = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    Y = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    x = mp(X);
    y = mp(Y);
    
    d = (x > y)-(X > Y);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = (x >= y)-(X >= Y);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
   
    d = (x < y)-(X < Y);    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

    d = (x <= y)-(X <= Y);    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;

    X(:,1) = Y(:,1);
    x(:,1) = y(:,1);
    
    d = (x == y)-(X == Y);    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;

    d = (x ~= y)-(X ~= Y);    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % and, or, not, xor
    fprintf('%-20s:','and,or,not,xor')    

    X = rand(Rows,Cols);
    Y = rand(Rows,Cols);
    X(:,1) = 0;
    X(1,:) = 0;
    Y(:,1) = 0;
    Y(1,:) = 0;

    x = mp(X);
    y = mp(Y);

    d = and(X,Y)-and(x,y);    
        if norm(d,1)>tolerance, fprintf(' and <- fail \n'), end;

    d = or(X,Y)-or(x,y);    
        if norm(d,1)>tolerance, fprintf(' or <- fail \n'), end;
        
    d = xor(X,Y)-xor(x,y);    
        if norm(d,1)>tolerance, fprintf(' xor <- fail \n'), end;

    d = not(X)-not(x);    
        if norm(d,1)>tolerance, fprintf(' not <- fail \n'), end;
        
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % all()
    fprintf('%-20s:','all()')    

    X = rand(Rows,Cols)-0.5;
    X(1,1) = 0;
    X(Rows-1,Cols-1) = 0;
    
    d = all(X)-all(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = all(X,2)-all(mp(X),2);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % any()
    fprintf('%-20s:','any()')    

    X = rand(Rows,Cols)-0.5;
    X(:,1) = 0;
    X(1,:) = 0;
    
    d = any(X)-any(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

    d = any(X,2)-any(mp(X),2);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % isinf   
    fprintf('%-20s:','isinf()')    

    X = magic(5);
    X(1,1) = inf+1i*inf;
    X(2,2) = nan+1i*nan;
    
    Y = mp(X);
    
    d = isinf(X) - isinf(Y);
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;    
      
    fprintf('\t<- success \n')

    % **************************************************************************
    % isnan
    fprintf('%-20s:','isnan()')    

    X = magic(5);
    X(1,1) = inf+1i*inf;
    X(2,2) = nan+1i*nan;
    
    Y = mp(X);

    d = isnan(X) - isnan(Y);
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;    
       
    fprintf('\t<- success \n')

    % **************************************************************************
    % isfinite   
    fprintf('%-20s:','isfinite()')    

    X = magic(5);
    X(1,1) = inf+1i*inf;
    X(2,2) = nan+1i*nan;
    
    Y = mp(X);

    d = isfinite(X) - isfinite(Y);
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;    
       
    fprintf('\t<- success \n')

    % **************************************************************************
    % isreal   
    fprintf('%-20s:','isreal()')    

    X = rand(Rows,Cols)+1i*rand(Rows,Cols);
    Y = mp(X);
    
    d = isreal(X) - isreal(Y);
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;    

    X = rand(Rows,Cols);
    Y = mp(X);
    
    d = isreal(X) - isreal(Y);
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;    
      
    fprintf('\t<- success \n')

    % **************************************************************************
    % abs()
    fprintf('%-20s:','abs()')    

    x = rand(Rows,Cols)-rand(Rows,Cols);
    X = mp(x);
    if norm(abs(X)-abs(x),1)>tolerance, fprintf(' 0 <- fail \n'), end;

    x = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    X = mp(x);
    if norm(abs(X)-abs(x),1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % sign()
    fprintf('%-20s:','sign()')    

    x = rand(Rows,Cols)-rand(Rows,Cols);
    X = mp(x);
    if norm(sign(X)-sign(x),1)>tolerance, fprintf(' 0 <- fail \n'), end;

    x = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    X = mp(x);
    if norm(sign(X)-sign(x),1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % conj()
    fprintf('%-20s:','conj()')    

    x = rand(Rows,Cols)-rand(Rows,Cols);
    X = mp(x);
    if norm(conj(X)-conj(x),1)>tolerance, fprintf(' 0 <- fail \n'), end;

    x = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    X = mp(x);
    if norm(conj(X)-conj(x),1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % angle()
    fprintf('%-20s:','angle()')    

    x = rand(Rows,Cols)-rand(Rows,Cols);
    X = mp(x);
    if norm(angle(X)-angle(x),1)>tolerance, fprintf(' 0 <- fail \n'), end;

    x = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    X = mp(x);
    if norm(angle(X)-angle(x),1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % imag()
    fprintf('%-20s:','imag()')    

    x = rand(Rows,Cols)-rand(Rows,Cols);
    X = mp(x);
    if norm(imag(X)-imag(x),1)>tolerance, fprintf(' 0 <- fail \n'), end;

    x = rand(Rows,Cols)-0.5+(rand(Rows,Cols)-0.5)*1i;
    X = mp(x);
    if norm(imag(X)-imag(x),1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % real()
    fprintf('%-20s:','real()')    

    x = rand(Rows,Cols)-rand(Rows,Cols);
    X = mp(x);
    if norm(real(X)-real(x),1)>tolerance, fprintf(' 0 <- fail \n'), end;

    x = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    X = mp(x);
    if norm(real(X)-real(x),1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % complex()
    fprintf('%-20s:','complex()')    

    x = rand(Rows,Cols)-0.5;
    y = rand(Rows,Cols)-0.5;    
    X = mp(x);
    Y = mp(y);    
    if norm(complex(X,Y)-complex(x,y),1)>tolerance, fprintf(' 0 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % ceil()
    fprintf('%-20s:','ceil()')    

    x = rand(Rows,Cols)-rand(Rows,Cols);
    X = mp(x);
    if norm(ceil(X)-ceil(x),1)>tolerance, fprintf(' 0 <- fail \n'), end;

    x = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    X = mp(x);
    if norm(ceil(X)-ceil(x),1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % fix()
    fprintf('%-20s:','fix()')    

    x = rand(Rows,Cols)-rand(Rows,Cols);
    X = mp(x);
    if norm(fix(X)-fix(x),1)>tolerance, fprintf(' 0 <- fail \n'), end;

    x = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    X = mp(x);
    if norm(fix(X)-fix(x),1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % floor()
    fprintf('%-20s:','floor()')    

    x = rand(Rows,Cols)-rand(Rows,Cols);
    X = mp(x);
    if norm(floor(X)-floor(x),1)>tolerance, fprintf(' 0 <- fail \n'), end;

    x = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    X = mp(x);
    if norm(floor(X)-floor(x),1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % idivide()
    fprintf('%-20s:','idivide()')    

    Y = randi([-100 100],Rows,Cols,'int32');
    X = randi([-10 10],Rows,Cols,'int32');
    
    % remove zeros - for compatibility with MATLAB
    Y(Y == 0) = int32(1);
    
    if norm(idivide(X,Y)         -idivide(mp(X),mp(Y)),1) > tolerance, fprintf(' 0 <- fail \n'), end;
    if norm(idivide(X,Y,'fix')   - idivide(mp(X),mp(Y),'fix'),  1) > tolerance, fprintf(' 1 <- fail \n'), end;
    if norm(idivide(X,Y,'round') - idivide(mp(X),mp(Y),'round'),1) > tolerance, fprintf(' 2 <- fail \n'), end;
    if norm(idivide(X,Y,'floor') - idivide(mp(X),mp(Y),'floor'),1) > tolerance, fprintf(' 3 <- fail \n'), end;
    if norm(idivide(X,Y,'ceil')  - idivide(mp(X),mp(Y),'ceil'), 1) > tolerance, fprintf(' 4 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % round()
    fprintf('%-20s:','round()')    

    x = rand(Rows,Cols)-rand(Rows,Cols);
    X = mp(x);
    if norm(round(X)-round(x),1)>tolerance, fprintf(' 0 <- fail \n'), end;

    x = rand(Rows,Cols)+rand(Rows,Cols)*1i;
    X = mp(x);
    if norm(round(X)-round(x),1)>tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % rem()
    fprintf('%-20s:','rem()')    

    X = rand(Rows,Cols)-0.5;
    Y = rand(Rows,Cols)-0.5;    

    if norm(rem(X,Y)-rem(mp(X),mp(Y)),1) > tolerance, fprintf(' 0 <- fail \n'), end;

    X = magic(Rows)-rand(Rows,Rows);
    Y = magic(Rows)-rand(Rows,Rows);    

    if norm(rem(X,Y)-rem(mp(X),mp(Y)),1) > tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % mod()
    fprintf('%-20s:','mod()')    

    X = rand(Rows,Cols)-0.5;
    Y = rand(Rows,Cols)-0.5;    

    if norm(mod(X,Y)-mod(mp(X),mp(Y)),1) > tolerance, fprintf(' 0 <- fail \n'), end;

    X = magic(Rows)-rand(Rows,Rows);
    Y = magic(Rows)-rand(Rows,Rows);    

    if norm(mod(X,Y)-mod(mp(X),mp(Y)),1) > tolerance, fprintf(' 1 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % Matrix tril()     
    fprintf('%-20s:','tril(matrix)')    

    X = rand(Rows,Cols);
    
    d = tril(X)-tril(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 0 <- fail \n'), end;
    
    d = tril(X,int32((Rows-1)/2))-tril(mp(X),int32((Rows-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 1 <- fail \n'), end;

    d = tril(X,int32((Cols-1)/2))-tril(mp(X),int32((Cols-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 2 <- fail \n'), end;

    d = tril(X,-int32((Rows-1)/2))-tril(mp(X),-int32((Rows-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 3 <- fail \n'), end;

    d = tril(X,-int32((Cols-1)/2))-tril(mp(X),-int32((Cols-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 4 <- fail \n'), end;

    fprintf('\t<- success \n')
    
    % **************************************************************************
    % Vector tril()     
    fprintf('%-20s:','tril(vector)')    

    X = rand(1,Cols);
    
    d = tril(X)-tril(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 0 <- fail \n'), end;
    
    d = tril(X,int32((Cols-1)/2))-tril(mp(X),int32((Cols-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 1 <- fail \n'), end;

    d = tril(X,0)-tril(mp(X),0);
        if max(max(abs(double(d))))>tolerance, fprintf(' 2 <- fail \n'), end;

    d = tril(X,-int32((Cols-1)/2))-tril(mp(X),-int32((Cols-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 3 <- fail \n'), end;
    
    X = rand(Rows,1);
    d = tril(X,-int32((Rows-1)/2))-tril(mp(X),-int32((Rows-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 4 <- fail \n'), end;

    d = tril(X,0)-tril(mp(X),0);
        if max(max(abs(double(d))))>tolerance, fprintf(' 5 <- fail \n'), end;

    d = tril(X,int32((Rows-1)/2))-tril(mp(X),int32((Rows-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 6 <- fail \n'), end;

    fprintf('\t<- success \n')

    % **************************************************************************
    % Matrix triu()     
    fprintf('%-20s:','triu(matrix)')    

    X = rand(Rows,Cols);
    
    d = triu(X)-triu(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 0 <- fail \n'), end;
    
    d = triu(X,int32((Rows-1)/2))-triu(mp(X),int32((Rows-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 1 <- fail \n'), end;

    d = triu(X,int32((Cols-1)/2))-triu(mp(X),int32((Cols-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 2 <- fail \n'), end;

    d = triu(X,-int32((Rows-1)/2))-triu(mp(X),-int32((Rows-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 3 <- fail \n'), end;

    d = triu(X,-int32((Cols-1)/2))-triu(mp(X),-int32((Cols-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 4 <- fail \n'), end;

    fprintf('\t<- success \n')
    
    % **************************************************************************
    % Vector triu()     
    fprintf('%-20s:','triu(vector)')    

    X = rand(1,Cols);
    
    d = triu(X)-triu(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 0 <- fail \n'), end;
    
    d = triu(X,int32((Cols-1)/2))-triu(mp(X),int32((Cols-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 1 <- fail \n'), end;

    d = triu(X,0)-triu(mp(X),0);
        if max(max(abs(double(d))))>tolerance, fprintf(' 2 <- fail \n'), end;

    d = triu(X,-int32((Cols-1)/2))-triu(mp(X),-int32((Cols-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 3 <- fail \n'), end;
    
    X = rand(Rows,1);
    d = triu(X,-int32((Rows-1)/2))-triu(mp(X),-int32((Rows-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 4 <- fail \n'), end;

    d = triu(X,0)-triu(mp(X),0);
        if max(max(abs(double(d))))>tolerance, fprintf(' 5 <- fail \n'), end;

    d = triu(X,int32((Rows-1)/2))-triu(mp(X),int32((Rows-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 6 <- fail \n'), end;

    fprintf('\t<- success \n')

    % **************************************************************************
    % Matrix diag()     
    fprintf('%-20s:','diag(matrix)')    

    X = rand(Rows,Cols);
    
    d = diag(X)-diag(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 0 <- fail \n'), end;
    
    d = diag(X,int32((Cols-1)/2))-diag(mp(X),int32((Cols-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 1 <- fail \n'), end;

    d = diag(X,-int32((Rows-1)/2))-diag(mp(X),-int32((Rows-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 2 <- fail \n'), end;

    fprintf('\t<- success \n')

    % **************************************************************************
    % Vector diag()     
    fprintf('%-20s:','diag(vector)')    

    X = rand(1,Cols);
    
    d = diag(X)-diag(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 0 <- fail \n'), end;
    
    d = diag(X,int32((Rows-1)/2))-diag(mp(X),int32((Rows-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 1 <- fail \n'), end;

    d = diag(X,int32((Cols-1)/2))-diag(mp(X),int32((Cols-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 2 <- fail \n'), end;

    d = diag(X,-int32((Rows-1)/2))-diag(mp(X),-int32((Rows-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 3 <- fail \n'), end;

    d = diag(X,-int32((Cols-1)/2))-diag(mp(X),-int32((Cols-1)/2));
        if max(max(abs(double(d))))>tolerance, fprintf(' 4 <- fail \n'), end;

    fprintf('\t<- success \n')
    
    % **************************************************************************    
    % Matrix norm() 
    fprintf('%-20s:','norm(matrix)')
    
    X = rand(Rows,Cols);
    
    d = norm(X)-norm(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 0 <- fail \n'), end;
    
    d = norm(X,1)-norm(mp(X),1);
        if max(max(abs(double(d))))>tolerance, fprintf(' 1 <- fail \n'), end;

    d = norm(X,Inf)-norm(mp(X),Inf);
        if max(max(abs(double(d))))>tolerance, fprintf(' 2 <- fail \n'), end;

    d = norm(X,'fro')-norm(mp(X),'fro');
        if max(max(abs(double(d))))>tolerance, fprintf(' 3 <- fail \n'), end;

    fprintf('\t<- success \n')
    
    % **************************************************************************
    % Vector norm()     
    fprintf('%-20s:','norm(vector)')    

    X = rand(1,Cols);
    
    d = norm(X)-norm(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 0 <- fail \n'), end;
    
    d = norm(X,1)-norm(mp(X),1);
        if max(max(abs(double(d))))>tolerance, fprintf(' 1 <- fail \n'), end;

    d = norm(X,Inf)-norm(mp(X),Inf);
        if max(max(abs(double(d))))>tolerance, fprintf(' 2 <- fail \n'), end;

    d = norm(X,-Inf)-norm(mp(X),-Inf);
        if max(max(abs(double(d))))>tolerance, fprintf(' 3 <- fail \n'), end;

    d = norm(X,'fro')-norm(mp(X),'fro');
        if max(max(abs(double(d))))>tolerance, fprintf(' 4 <- fail \n'), end;

    d = norm(X,3.14)-norm(mp(X),3.14);
        if max(max(abs(double(d))))>tolerance, fprintf(' 5 <- fail \n'), end;

    fprintf('\t<- success \n')
    
    % **************************************************************************    
    % cond() 
    fprintf('%-20s:','cond()')
    
    X = rand(Rows,Cols);
    
    d = cond(X)-cond(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 0 <- fail \n'), end;
        
    X = rand(Rows,Rows);    
    
    d = cond(X,1)-cond(mp(X),1);
        if max(max(abs(double(d))))>tolerance, fprintf(' 1 <- fail \n'), end;

    d = cond(X,Inf)-cond(mp(X),Inf);
        if max(max(abs(double(d))))>tolerance, fprintf(' 2 <- fail \n'), end;

    d = cond(X,'fro')-cond(mp(X),'fro');
        if max(max(abs(double(d))))>tolerance, fprintf(' 3 <- fail \n'), end;

    X = rand(Rows,Cols)*1i;
    
    d = cond(X)-cond(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 4 <- fail \n'), end;

    X = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;        

    d = cond(X,1)-cond(mp(X),1);
        if max(max(abs(double(d))))>tolerance, fprintf(' 5 <- fail \n'), end;

    d = cond(X,Inf)-cond(mp(X),Inf);
        if max(max(abs(double(d))))>tolerance, fprintf(' 6 <- fail \n'), end;

    d = cond(X,'fro')-cond(mp(X),'fro');
        if max(max(abs(double(d))))>tolerance, fprintf(' 7 <- fail \n'), end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************    
    % rcond() 
    fprintf('%-20s:','rcond()')
    
    X = rand(Rows,Rows);
    d = rcond(X)-rcond(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 0 <- fail \n'), end;
        
    X = rand(Rows,Rows)*1i;
    d = rcond(X)-rcond(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 1 <- fail \n'), end;

    X = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;        
    d = rcond(X)-rcond(mp(X));    
        if max(max(abs(double(d))))>tolerance, fprintf(' 2 <- fail \n'), end;
    
    fprintf('\t<- success \n')

    % **************************************************************************
    % factorial()     
    fprintf('%-20s:','factorial()')    

    X = magic(3);
    
    d = factorial(X)-factorial(mp(X));    
        if norm(d,1) > tolerance, fprintf(' 0 <- fail \n'), end;
        
    fprintf('\t<- success \n')    

    % **************************************************************************
    % mean()     
    fprintf('%-20s:','mean()')    

    X = rand(Rows,Cols);
    
    d = mean(X) - mean(mp(X));    
        if norm(d,1) > tolerance, fprintf(' 0 <- fail \n'), end;

	d = mean(X,2) - mean(mp(X),2);    
        if norm(d,1) > tolerance, fprintf(' 1 <- fail \n'), end;
		
    X = rand(Rows,Cols) + rand(Rows,Cols)*1i;     
	
    d = mean(X) - mean(mp(X));    
        if norm(d,1) > tolerance, fprintf(' 2 <- fail \n'), end;

	d = mean(X,2) - mean(mp(X),2);    
        if norm(d,1) > tolerance, fprintf(' 3 <- fail \n'), end;

	X = rand(1,Cols);
    
    d = mean(X) - mean(mp(X));    
        if norm(d,1) > tolerance, fprintf(' 4 <- fail \n'), end;

	d = mean(X,1) - mean(mp(X),1);    
        if norm(d,1) > tolerance, fprintf(' 5 <- fail \n'), end;

	d = mean(X,2) - mean(mp(X),2);    
        if norm(d,1) > tolerance, fprintf(' 6 <- fail \n'), end;

	X = rand(Rows,1);
    
    d = mean(X) - mean(mp(X));    
        if norm(d,1) > tolerance, fprintf(' 7 <- fail \n'), end;

	d = mean(X,1) - mean(mp(X),1);    
        if norm(d,1) > tolerance, fprintf(' 8 <- fail \n'), end;

	d = mean(X,2) - mean(mp(X),2);    
        if norm(d,1) > tolerance, fprintf(' 9 <- fail \n'), end;
		
	X = rand(1,Cols)*1i;
    
    d = mean(X) - mean(mp(X));    
        if norm(d,1) > tolerance, fprintf(' 10 <- fail \n'), end;

	d = mean(X,1) - mean(mp(X),1);    
        if norm(d,1) > tolerance, fprintf(' 11 <- fail \n'), end;

	d = mean(X,2) - mean(mp(X),2);    
        if norm(d,1) > tolerance, fprintf(' 12 <- fail \n'), end;

	X = rand(Rows,1)*1i;
    
    d = mean(X) - mean(mp(X));    
        if norm(d,1) > tolerance, fprintf(' 13 <- fail \n'), end;

	d = mean(X,1) - mean(mp(X),1);    
        if norm(d,1) > tolerance, fprintf(' 14 <- fail \n'), end;

	d = mean(X,2) - mean(mp(X),2);    
        if norm(d,1) > tolerance, fprintf(' 15 <- fail \n'), end;
	
    fprintf('\t<- success \n')    
    
    % **************************************************************************
    % std()     
    fprintf('%-20s:','std()')    

    X = rand(Rows,Cols);
    
    d = std(X) - std(mp(X));    
        if norm(d,1) > tolerance, fprintf(' 0 <- fail \n'), end;

	d = std(X,1) - std(mp(X),1);    
        if norm(d,1) > tolerance, fprintf(' 1 <- fail \n'), end;

    d = std(X,1,2) - std(mp(X),1,2);    
        if norm(d,1) > tolerance, fprintf(' 2 <- fail \n'), end;
		
    X = rand(Rows,Cols) + rand(Rows,Cols)*1i;     
	
    d = std(X) - std(mp(X));    
        if norm(d,1) > tolerance, fprintf(' 3 <- fail \n'), end;

	d = std(X,1) - std(mp(X),1);    
        if norm(d,1) > tolerance, fprintf(' 4 <- fail \n'), end;

    d = std(X,1,2) - std(mp(X),1,2);    
        if norm(d,1) > tolerance, fprintf(' 5 <- fail \n'), end;
	
    fprintf('\t<- success \n')    

    
    % **************************************************************************
    % erf()     
    fprintf('%-20s:','erf()')    

    X = rand(Rows,Cols)-1/2;
    
    d = erf(X)-erf(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % erfc()     
    fprintf('%-20s:','erfc()')    

    X = rand(Rows,Cols)-1/2;
    
    d = erfc(X)-erfc(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % erfi()     
    fprintf('%-20s:','erfi()')    

    X = mp(rand(Rows,Cols))-1/2;
    
    d = erfi(X)+i*erf(i*X);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % fresnels()     
    fprintf('%-20s:','fresnels()')    

	% Real
    x = mp(rand(Rows,Cols)-0.5);
	
	d = fresnels(-x) + fresnels(x);   % = 0
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
		
	d = fresnels(i*x) + i*fresnels(x);   % = 0
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

	% Complex
    z = mp(rand(Rows,Cols)-0.5) + mp(rand(Rows,Cols)-0.5)*i;
    
    d = -i*fresnels(z) - fresnels(i*z);    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % fresnelc()     
    fprintf('%-20s:','fresnelc()')    

	% Real
    x = mp(rand(Rows,Cols)-0.5);

	d = fresnelc(-x) + fresnelc(x);   		% = 0
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	d = fresnelc(i*x) - i*fresnelc(x);   	% = 0
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
	
	% Complex
    z = mp(rand(Rows,Cols)-0.5) + mp(rand(Rows,Cols)-0.5)*i;
    
    d = fresnelc(i*z) - i*fresnelc(z);    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
        
	ksi = mp('sqrt(pi)*(1-i)/2')*z;
	d = fresnelc(z) + i*fresnels(z) - (erf(ksi)*mp('(1+i)/2'));
		if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
	
    fprintf('\t<- success \n')

    % **************************************************************************
    % gammaln()     
    fprintf('%-20s:','gammaln()')    

    X = rand(Rows,Cols);
    
    d = gammaln(X)-gammaln(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % gamma()     
    fprintf('%-20s:','gamma()')    

    X = rand(Rows,Cols);
    
    d = gamma(X)-gamma(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % gammainc()     
    fprintf('%-20s:','gammainc()')    

    X = rand(Rows,Cols);
	A = rand(Rows,Cols);
    
    d = gammainc(X,A)-gammainc(mp(X),mp(A));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	d = gammainc(X,A,'upper')-gammainc(mp(X),mp(A),'upper');    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % psi() - polygamma (digamma)    
    fprintf('%-20s:','psi()')    

    X = 100 * rand(Rows,Cols);
    
	% Nonnegative real arguments - compare with Matlab
    d = psi(X)-psi(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
		
	% Negative real arguments - compare with Matlab 
	% Matlab doesn't support x < 0 - so we use reflection formula
	d = (psi(X) + 1./X + pi * cot(pi*X)) - psi(-mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
	
	% Special Values
	x = mp('1'); % psi(1) = -0.5772 =  Euler constant
	
    d = - psi(x) - mp('euler');
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;
	
	x = mp('-1/2'); % psi(-1/2) = 2-0.5772..-2*ln(2), Euler constant
	
    d = psi(x) - mp('2 - euler - 2 * ln2'); % Use native evaluator
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;

	d = psi(x) - (2 - mp('euler') - 2 * log(mp('2'))); 
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;

	% Complex argument of large magnitude
	x = mp('(1-i)*1e+6'); 
	y = mp('14.162083898244246758816564786751940196311025665619-0.78539841339748997628232751248654238771546298476441*i');

	d = psi(x) - y; 
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), end;
	
    fprintf('\t<- success \n')

    % **************************************************************************
    % zeta()
    fprintf('%-20s:','zeta()')    
   
    X = -10:1:10;
    Y = [0, -0.00757575757575758, 0, 0.00416666666666667, 0,  -0.00396825396825397,  0,  0.00833333333333333,   0, -0.0833333333333333,  -0.5,   Inf,  1.64493406684823,  1.20205690315959,   1.08232323371114,  1.03692775514337,  1.01734306198445,  1.00834927738192, 1.00407735619794,      1.00200839282608,  1.00099457512782];
    d = Y-zeta(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % expint() - Exponential Integral
    if true % chebfun overloads the expint and breaks the code!!!
    fprintf('%-20s:','expint()')    

    X = rand(Rows,Cols);
    
    d = expint(X)-expint(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	X = rand(Rows,Cols)*1i;
    
    d = expint(X)-expint(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    fprintf('\t<- success \n')
    end
    
    % **************************************************************************
    % eint() - Exponential Integral - principal value
    fprintf('%-20s:','eint()')    

    X = rand(Rows,Cols);
    
	% 5.1.7
    d = expint(-mp(X))+eint(mp(X))+mp('i*pi');    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % logint() - Logarithmic Integral
    fprintf('%-20s:','logint()')    

    X = 1+rand(Rows,Cols);
    
	% 5.1.3
    d = logint(mp(X))-eint(log(mp(X)));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % cosint() - Cosine Integral
    fprintf('%-20s:','cosint()')    

    X = mp(rand(Rows,Cols));
    
	% 5.2.23
    d = cosint(X)+(expint(X*1i)+expint(-X*1i))/2;    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	X = mp(rand(Rows,Cols))*1i;
    
    d = cosint(X)+(expint(X*1i)+expint(-X*1i))/2;    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    fprintf('\t<- success \n')
	
    % **************************************************************************
    % sinint() - Sine Integral
    fprintf('%-20s:','sinint()')    

    X = mp(rand(Rows,Cols));
    
	% 5.2.21
    d = sinint(X)-(expint(X*1i)-expint(-X*1i))/(2*1i)-mp('pi/2');    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	X = mp(rand(Rows,Cols))*1i;
    
    d = sinint(X)-(expint(X*1i)-expint(-X*1i))/(2*1i)-mp('pi/2');    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % besselj()
    fprintf('%-20s:','besselj()')    

    nu = mp(rand(Rows,Cols));
	z  = mp(rand(Rows,Cols));
    
    d = besselj(nu,z)-besselj(double(nu),double(z));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	d = besselj(nu,z,1)-besselj(double(nu),double(z),1);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    fprintf('\t<- success \n')
	
    % **************************************************************************
    % bessely()
    fprintf('%-20s:','bessely()')    

    nu = mp(rand(Rows,Cols));
	z  = mp(rand(Rows,Cols));
    
    d = bessely(nu,z)-bessely(double(nu),double(z));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	d = bessely(nu,z,1)-bessely(double(nu),double(z),1);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % besseli()
    fprintf('%-20s:','besseli()')    

    nu = mp(rand(Rows,Cols));
	z  = mp(rand(Rows,Cols));
    
    d = besseli(nu,z)-besseli(double(nu),double(z));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	d = besseli(nu,z,1)-besseli(double(nu),double(z),1);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % besselk()
    fprintf('%-20s:','besselk()')    

    nu = mp(rand(Rows,Cols));
	z  = mp(rand(Rows,Cols));

    d = besselk(nu,z)-besselk(double(nu),double(z));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	d = besselk(nu,z,1)-besselk(double(nu),double(z),1);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % besselh()
    fprintf('%-20s:','besselh()')    

    nu = mp(rand(Rows,Cols));
	z  = mp(rand(Rows,Cols));

    d = besselh(nu,z)-besselh(double(nu),double(z));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	d = besselh(nu,2,z)-besselh(double(nu),2,double(z));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

	d = besselh(nu,1,z,1)-besselh(double(nu),1,double(z),1);    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

	d = besselh(nu,2,z,1)-besselh(double(nu),2,double(z),1);    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
        
    fprintf('\t<- success \n')

    % **************************************************************************
    % hypergeom()
    fprintf('%-20s:','hypergeom()')    

	z = mp(rand(Rows,Cols)-1/2) + i*mp(rand(Rows,Cols)-1/2);

	% 0F0 = exp(z)
    d = hypergeom([],[],z)-exp(z);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	% 1F0 = (1-z)^-a
	a = mp('1-i');
    d = hypergeom(a,[],z)-(1-z).^(-a);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;

	% 0F1, besselj(nu,z) = 0F1(;nu+a,-z*z/4)*(z/2)^nu/gamma(nu+1)
	nu = mp('1-i'); % note we support complex arguments in hypergeom & gamma & bessel
    d = (z./2).^nu .* hypergeom([],nu + 1,-z.*z./4)./gamma(nu+1)-besselj(nu,z);    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), end;

	% 1F1, erfc
	d = (1 - 2 * z .* hypergeom(1, 1.5, z.*z)./(mp('sqrt(pi)')*exp(z.*z))) - erfc(z); 
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), end;
		
	% 1F1, gammainc
	a = mp('1-i');
	d = (z.^a).*hypergeom(1,1 + a,z)./(a.*exp(z).*gamma(a)) - gammainc(z,a); 
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), end;
		
    % add your special cases 
	% ...
    fprintf('\t<- success \n')

    % **************************************************************************
    % kummerM() - same as hypergeom
    fprintf('%-20s:','kummerM()')    

	a = mp(rand(Rows,Cols)-1/2) + i*mp(rand(Rows,Cols)-1/2);	
	b = mp(rand(Rows,Cols)-1/2) + i*mp(rand(Rows,Cols)-1/2);		
	z = mp(rand(Rows,Cols)-1/2) + i*mp(rand(Rows,Cols)-1/2);

	% Kummer transformation
    d = kummerM(a,b,z)-exp(z).*kummerM(b-a,b,-z);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;
		
    fprintf('\t<- success \n')

    % **************************************************************************
    % kummerU()
    fprintf('%-20s:','kummerU()')    

	a = mp(rand(Rows,Cols)-1/2); %+ i*mp(rand(Rows,Cols)-1/2);	
	b = mp(rand(Rows,Cols)-1/2); %+ i*mp(rand(Rows,Cols)-1/2);		
	z = mp(rand(Rows,Cols)-1/2); %+ i*mp(rand(Rows,Cols)-1/2);

	% Transformations
    d = kummerU(a,b,z) - z.^(1-b).*kummerU(a - b + 1,2 - b,z);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), end;

	% Expression using kummerM
	U = gamma(1 - b) .* kummerM(a, b, z) ./ gamma(a - b + 1) + z.^(1 - b) .* gamma(b - 1) .* kummerM(1 + a - b, 2 - b, z) ./ gamma(a);	
	
	d = U - kummerU(a,b,z); 
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), end;
		
    fprintf('\t<- success \n')

    % **************************************************************************
    % expm()     
    fprintf('%-20s:','expm()')    

    X = rand(Rows,Rows)-0.5;
    
    d = expm(X)-expm(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1), X, end;
   
    X = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;
    
    d = expm(X)-expm(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), X, end;

    fprintf('\t<- success \n')
    
    % **************************************************************************
    % logm()     
    fprintf('%-20s:','logm()')    

    X = rand(Rows,Rows)-0.5;
    X = X*X';
    
    Y = logm(mp(X)); 

    d = X - expm(double(Y));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1), X, end;
    
    d = X - expm(Y);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), X, end;

    X = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;
    Y = logm(mp(X)); 
    
    d = X - expm(double(Y));    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), norm(d,1), X, end;
    
    d = X - expm(Y);    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), norm(d,1), X, end;
   
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % sqrtm()     
    fprintf('%-20s:','sqrtm()')    

    X = rand(Rows,Rows)-0.5;
    Y = sqrtm(mp(X));
    
    d = Y * Y - X;    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1), X, end;

    X = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;
    Y = sqrtm(mp(X));
    
    d = Y * Y - X;    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), X, end;
    
    fprintf('\t<- success \n')
    
    % **************************************************************************
    % sinm()     
    fprintf('%-20s:','sinm()')    

    X = rand(Rows,Rows)-0.5;
    
    d = funm(X,@sin)-sinm(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1), X, end;
   
    X = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;
    
    d = funm(X,@sin)-sinm(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), X, end;

    fprintf('\t<- success \n')

    % **************************************************************************
    % cosm()     
    fprintf('%-20s:','cosm()')    

    X = rand(Rows,Rows)-0.5;
    
    d = funm(X,@cos)-cosm(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1), X, end;
   
    X = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;
    
    d = funm(X,@cos)-cosm(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), X, end;

    fprintf('\t<- success \n')
   
    % **************************************************************************
    % sinhm()     
    fprintf('%-20s:','sinhm()')    

    X = rand(Rows,Rows)-0.5;
    
    d = funm(X,@sinh)-sinhm(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1), X, end;
   
    X = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;
    
    d = funm(X,@sinh)-sinhm(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), X, end;

    fprintf('\t<- success \n')
    
    % **************************************************************************
    % coshm()     
    fprintf('%-20s:','coshm()')    

    X = rand(Rows,Rows)-0.5;
    
    d = funm(X,@cosh)-coshm(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1), X, end;
   
    X = rand(Rows,Rows)-0.5+(rand(Rows,Rows)-0.5)*1i;
    
    d = funm(X,@cosh)-coshm(mp(X));    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), X, end;

    fprintf('\t<- success \n')

    % **************************************************************************
    % fft     
    fprintf('%-20s:','fft')    

    N = 32;
    X = rand(N,N)-0.5;
    Y = mp(X);
    
    d = fft(X)-fft(Y);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1), X, end;

    d = fft(X,N/2)-fft(Y,N/2);    % Truncation
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), X, end;

    d = fft(X,2*N)-fft(Y,2*N);    % Zero padding
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), norm(d,1), X, end;
        
    d = fft(X,N-1)-fft(Y,N-1);    % Non-even/power of 2
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), norm(d,1), X, end;

    d = fft(X,[],2)-fft(Y,[],2);    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), norm(d,1), X, end;

    d = fft(X,N/2,2)-fft(Y,N/2,2);    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), norm(d,1), X, end;
        
    d = fft(X,2*N,2)-fft(Y,2*N,2);    
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), norm(d,1), X, end;
        
    d = fft(X,N-1,2)-fft(Y,N-1,2);    
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), norm(d,1), X, end;
   
    d = Y-ifft(fft(Y));
        if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), norm(d,1), X, end;
   
    X = rand(N,N)-0.5+(rand(N,N)-0.5)*1i;
    Y = mp(X);
    
    d = fft(X)-fft(Y);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1), X, end;

    d = fft(X,N/2)-fft(Y,N/2);    
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), X, end;

    d = fft(X,2*N)-fft(Y,2*N);    
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), norm(d,1), X, end;
        
    d = fft(X,N-1)-fft(Y,N-1);    
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), norm(d,1), X, end;

    d = fft(X,[],2)-fft(Y,[],2);    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), norm(d,1), X, end;

    d = fft(X,N/2,2)-fft(Y,N/2,2);    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), norm(d,1), X, end;
        
    d = fft(X,2*N,2)-fft(Y,2*N,2);    
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), norm(d,1), X, end;
        
    d = fft(X,N-1,2)-fft(Y,N-1,2);    
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), norm(d,1), X, end;
   
    d = Y-ifft(fft(Y));
        if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), norm(d,1), X, end;
        
    fprintf('\t<- success \n')
      
    % **************************************************************************
    % ifft     
    fprintf('%-20s:','ifft')    

    N = 16;
    A = rand(N,N)-0.5;
    X = fft(A);
    Y = fft(mp(A));
    
    d = ifft(X)-ifft(Y);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,N/2)-ifft(Y,N/2);    % Truncation
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,2*N)-ifft(Y,2*N);    % Zero padding
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,N-1)-ifft(Y,N-1);    % Non-even/power of 2
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,[],2)-ifft(Y,[],2);    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,N/2,2)-ifft(Y,N/2,2);    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,2*N,2)-ifft(Y,2*N,2);    
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,N-1,2)-ifft(Y,N-1,2);    
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), norm(d,1), X; end;
   
    d = ifft(X,'symmetric')-ifft(Y,'symmetric');
        if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), norm(d,1), X; end;
   
    d = ifft(X,N/2,'symmetric')-ifft(Y,N/2,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,2*N,'symmetric')-ifft(Y,2*N,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 10 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,N-1,'symmetric')-ifft(Y,N-1,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 11 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,[],2,'symmetric')-ifft(Y,[],2,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 12 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,N/2,2,'symmetric')-ifft(Y,N/2,2,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 13 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,2*N,2,'symmetric')-ifft(Y,2*N,2,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 14 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,N-1,2,'symmetric')-ifft(Y,N-1,2,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 15 <- fail \n'), norm(d,1), X; end;
       
    A = rand(N,N)-0.5+(rand(N,N)-0.5)*1i;
    X = fft(A);
    Y = fft(mp(A));
    
    d = ifft(X)-ifft(Y);    
        if norm(d,1)>tolerance, fprintf(' 16 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,N/2)-ifft(Y,N/2);    
        if norm(d,1)>tolerance, fprintf(' 17 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,2*N)-ifft(Y,2*N);    
        if norm(d,1)>tolerance, fprintf(' 18 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,N-1)-ifft(Y,N-1);    
        if norm(d,1)>tolerance, fprintf(' 19 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,[],2)-ifft(Y,[],2);    
        if norm(d,1)>tolerance, fprintf(' 20 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,N/2,2)-ifft(Y,N/2,2);    
        if norm(d,1)>tolerance, fprintf(' 21 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,2*N,2)-ifft(Y,2*N,2);    
        if norm(d,1)>tolerance, fprintf(' 22 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,N-1,2)-ifft(Y,N-1,2);    
        if norm(d,1)>tolerance, fprintf(' 23 <- fail \n'), norm(d,1), X; end;
   
    d = ifft(X,'symmetric')-ifft(Y,'symmetric');
        if norm(d,1)>tolerance, fprintf(' 24 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,N/2,'symmetric')-ifft(Y,N/2,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 25 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,2*N,'symmetric')-ifft(Y,2*N,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 26 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,N-1,'symmetric')-ifft(Y,N-1,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 27 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,[],2,'symmetric')-ifft(Y,[],2,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 28 <- fail \n'), norm(d,1), X; end;

    d = ifft(X,N/2,2,'symmetric')-ifft(Y,N/2,2,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 29 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,2*N,2,'symmetric')-ifft(Y,2*N,2,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 30 <- fail \n'), norm(d,1), X; end;
        
    d = ifft(X,N-1,2,'symmetric')-ifft(Y,N-1,2,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 31 <- fail \n'), norm(d,1), X; end;
        
    fprintf('\t<- success \n')
   
    % **************************************************************************
    % fft2     
    fprintf('%-20s:','fft2')    

    N = 32;
    X = rand(N,N)-0.5;
    Y = mp(X);
    
    d = fft2(X) - fft2(Y);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1), X; end;

    d = fft2(X, N/2, N/2) - fft2(Y, N/2, N/2);    % Truncation
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), X; end;

    d = fft2(X, 2*N, 2*N) - fft2(Y, 2*N, 2*N);    % Zero padding
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), norm(d,1), X; end;
        
    d = fft2(X, N-1, N-1) - fft2(Y, N-1, N-1);    % Non-even/power of 2
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), norm(d,1), X; end;

    d = Y-ifft2(fft2(Y));
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), norm(d,1), X; end;
   
    X = rand(N,N)-0.5+(rand(N,N)-0.5)*1i;
    Y = mp(X);
    
    d = fft2(X) - fft2(Y);    
        if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), norm(d,1), X; end;

    d = fft2(X, N/2, N/2) - fft2(Y, N/2, N/2);    % Truncation
        if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), norm(d,1), X; end;

    d = fft2(X, 2*N, 2*N) - fft2(Y, 2*N, 2*N);    % Zero padding
        if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), norm(d,1), X; end;
        
    d = fft2(X, N-1, N-1) - fft2(Y, N-1, N-1);    % Non-even/power of 2
        if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), norm(d,1), X; end;

    d = Y-ifft2(fft2(Y));
        if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), norm(d,1), X; end;
        
    fprintf('\t<- success \n')
   
    % **************************************************************************
    % ifft2     
    fprintf('%-20s:','ifft2')    

    N = 16;
    A = rand(N,N)-0.5;
    X = fft2(A);
    Y = fft2(mp(A));
    
    d = ifft2(X)-ifft2(Y);    
        if norm(d,1)>tolerance, fprintf(' 0 <- fail \n'), norm(d,1),  end;

    d = ifft2(X, N/2, N/2)-ifft2(Y, N/2, N/2);    % Truncation
        if norm(d,1)>tolerance, fprintf(' 1 <- fail \n'), norm(d,1), end;

    d = ifft2(X, 2*N, 2*N)-ifft2(Y, 2*N, 2*N);    % Zero padding
        if norm(d,1)>tolerance, fprintf(' 2 <- fail \n'), norm(d,1), end;
        
    d = ifft2(X, N-1, N-1)-ifft2(Y, N-1, N-1);    % Non-even/power of 2
        if norm(d,1)>tolerance, fprintf(' 3 <- fail \n'), norm(d,1), end;

    d = ifft2(X,'symmetric')-ifft2(Y,'symmetric');    
        if norm(d,1)>tolerance, fprintf(' 4 <- fail \n'), norm(d,1), end;

    %d = ifft2(X, N/2, N/2,'symmetric')-ifft2(Y, N/2, N/2,'symmetric');    % Truncation
    %    if norm(d,1)>tolerance, fprintf(' 5 <- fail \n'), norm(d,1), end;

    %d = ifft2(X, 2*N, 2*N,'symmetric')-ifft2(Y, 2*N, 2*N,'symmetric');    % Zero padding
    %    if norm(d,1)>tolerance, fprintf(' 6 <- fail \n'), norm(d,1), end;
        
    %d = ifft2(X, N-1, N-1,'symmetric')-ifft2(Y, N-1, N-1,'symmetric');    % Non-even/power of 2
    %    if norm(d,1)>tolerance, fprintf(' 7 <- fail \n'), norm(d,1), end;
   
    A = rand(N,N)-0.5+(rand(N,N)-0.5)*1i;
    X = fft2(A);
    Y = fft2(mp(A));
    
    d = ifft2(X)-ifft2(Y);    
        if norm(d,1)>tolerance, fprintf(' 8 <- fail \n'), norm(d,1), end;

    d = ifft2(X, N/2, N/2)-ifft2(Y, N/2, N/2);    % Truncation
        if norm(d,1)>tolerance, fprintf(' 9 <- fail \n'), norm(d,1), end;

    d = ifft2(X, 2*N, 2*N)-ifft2(Y, 2*N, 2*N);    % Zero padding
        if norm(d,1)>tolerance, fprintf(' 10 <- fail \n'), norm(d,1), end;
        
    d = ifft2(X, N-1, N-1)-ifft2(Y, N-1, N-1);    % Non-even/power of 2
        if norm(d,1)>tolerance, fprintf(' 11 <- fail \n'), norm(d,1), end;

    %d = ifft2(X,'symmetric')-ifft2(Y,'symmetric');    
    %    if norm(d,1)>tolerance, fprintf(' 12 <- fail \n'), norm(d,1), end;

    %d = ifft2(X, N/2, N/2,'symmetric')-ifft2(Y, N/2, N/2,'symmetric');    % Truncation
    %    if norm(d,1)>tolerance, fprintf(' 13 <- fail \n'), norm(d,1), end;

    %d = ifft2(X, 2*N, 2*N,'symmetric')-ifft2(Y, 2*N, 2*N,'symmetric');    % Zero padding
    %    if norm(d,1)>tolerance, fprintf(' 14 <- fail \n'), norm(d,1), end;
        
    %d = ifft2(X, N-1, N-1,'symmetric')-ifft2(Y, N-1, N-1,'symmetric');    % Non-even/power of 2
    %    if norm(d,1)>tolerance, fprintf(' 15 <- fail \n'), norm(d,1), end;
        
    fprintf('\t<- success \n')
    
end

