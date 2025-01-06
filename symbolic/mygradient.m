function [gradF] = mygradient(F,var)
% Calculate symbolic gradient of scalar/vector/matrix function

% Get number of rows and columns
[nrow,ncol]=size(F);

% Get number of function variables
nfvar=length(symvar(F));

% Force variables to be a vector if a matrix is given
isvarmatrix=false;
sizevar=size(var);
if sizevar(1)>1 && sizevar(2)>1
  trasvar=transpose(var);
  var=trasvar(1:end);
  isvarmatrix=true;
end

% Get number of variables
nvar=length(var);

% Don't allow more components than number of function variables
if nrow>nfvar
  error('Number of rows greater than number of function variables not allowed!');
end
if ncol>nfvar
  error('Number of columns greater than number of function variables not allowed!');
end

% Check if scalar, vector or matrix
if nrow==1 && ncol==1
  whatis='scalar';
else
  if nrow==1 || ncol==1
    whatis='vector';
  else
    whatis='matrix';
  end
end

% Check if functional variables
funvar=setdiff(var,symvar(var));
isfunvar=not(isempty(funvar));

% Force vertical vectors
if strcmp(whatis,'vector') && ncol>1
  nrow=ncol;
  ncol=1;
  F=transpose(F);
end

% Calculate gradient
switch whatis
  case 'scalar'
    gradF=sym(zeros(1,nvar));
    for n=1:nvar
      if isfunvar
        gradF(n)=functionalDerivative(F,var(n));
      else
        gradF(n)=diff(F,var(n));
      end
    end
    if isvarmatrix
      gradF=transpose(reshape(gradF,[sqrt(nvar),sqrt(nvar)]));
    end
  case 'vector'
    gradF=sym(zeros(nrow,nvar));
    for i=1:nrow
      for n=1:nvar
        if isfunvar
          gradF(i,n)=functionalDerivative(F(i),var(n));
        else
          gradF(i,n)=diff(F(i),var(n));
        end
      end
    end
  case 'matrix'
    gradF=sym(zeros(nrow,ncol,nvar));
    for i=1:nrow
      for j=1:ncol
        for n=1:nvar
          if isfunvar
            gradF(i,j,n)=functionalDerivative(F(i,j),var(n));
          else
            gradF(i,j,n)=diff(F(i,j),var(n));
          end
        end
      end
    end
end

end