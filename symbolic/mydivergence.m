function [divF] = mydivergence(F,var)
% Calculate symbolic divergence of scalar/vector/matrix function

% Get number of rows and columns
[nrow,ncol]=size(F);

% Get number of function variables
nfvar=length(symvar(F));

% Force variables to be a vector if a matrix is given
sizevar=size(var);
if sizevar(1)>1 && sizevar(2)>1
  trasvar=transpose(var);
  var=trasvar(1:end);
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

% Calculate divergence
switch whatis
  case 'scalar'
    if 1==nvar
      if isfunvar
        divF=functionalDerivative(F,var(1));
      else
        divF=diff(F,var(1));
      end
    else
      error('Number of variables greater than one!');
    end
  case 'vector'
    if nrow==nvar
      divF=sym(zeros(1,1));
      for n=1:nvar
        if isfunvar
          divF=divF+functionalDerivative(F(n),var(n));
        else
          divF=divF+diff(F(n),var(n));
        end
      end
    else
      error('Number of rows not matching number of variables!');
    end
  case 'matrix'
    if ncol==nvar
      divF=sym(zeros(nrow,1));
      for i=1:nrow
        for n=1:nvar
          if isfunvar
            divF(i)=divF(i)+functionalDerivative(F(i,n),var(n));
          else
            divF(i)=divF(i)+diff(F(i,n),var(n));
          end
        end
      end
    else
      error('Number of columns not matching number of variables!');
    end
end

end