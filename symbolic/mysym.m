function [F] = mysym(name,size,fvar)
% Crate symbolic scalar/vector/matrix function

% Get number of rows and columns
nrow=size(1);
ncol=size(2);

% Get number of function variables
nfvar=length(fvar);

% Don't allow more components than number of variables
if nrow>nfvar
  error('Number of rows greater than number of variables not allowed!');
end
if ncol>nfvar
  error('Number of columns greater than number of variables not allowed!');
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

% Function variables with brackets
fvarfun=['(',char(fvar(1))];
for n=2:nfvar
  fvarfun=[fvarfun,',',char(fvar(n))];
end
fvarfun=[fvarfun,')'];

% Define symbolic scalar/vector/matrix function
switch whatis
  case 'scalar'
    F=str2sym(sprintf('%s%s',name,fvarfun));
  case 'vector'
    F=sym(zeros(nrow,ncol));
    for i=1:max(nrow,ncol)
      F(i)=str2sym(sprintf('%s_%s%s',name,fvar(i),fvarfun));
    end
  case 'matrix'
    F=sym(zeros(nrow,ncol));
    for i=1:nrow
      for j=1:ncol
        F(i,j)=str2sym(sprintf('%s_%s%s%s',name,fvar(i),fvar(j),fvarfun));
      end
    end
end

end