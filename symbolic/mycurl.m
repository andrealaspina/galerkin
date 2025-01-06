function [curlF] = mycurl(F,var)
% Calculate symbolic curl of scalar/vector function

% Get number of rows and columns
[nrow,ncol]=size(F);

% Get number of function variables
nfvar=length(symvar(F));

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
    error('Matrix functions not allowed!');
  end
end

% Check if functional variables
funvar=setdiff(var,symvar(var));
isfunvar=not(isempty(funvar));

% Force vertical vectors
if strcmp(whatis,'vector') && ncol>1
  nrow=ncol;
  F=transpose(F);
end

% Allow only 1D, 2D and 3D functions
if nrow>3
  error('Function size not allowed!');
end

% Allow only 2 and 3 (functional) variables
if max(length(funvar),nvar)<2 || max(length(funvar),nvar)>3
  error('Variables size not allowed!');
end

% Calculate curl
switch whatis
  case 'scalar'
    if nvar==2 || nvar==3
      if isfunvar
        curlF=[+functionalDerivative(F,var(2))
               -functionalDerivative(F,var(1))];
      else
        curlF=[+diff(F,var(2))
               -diff(F,var(1))];
      end
    end
  case 'vector'
    if nrow==2
      if isfunvar
        curlF=-functionalDerivative(F(1),var(2))+functionalDerivative(F(2),var(1));
      else
        curlF=-diff(F(1),var(2))+diff(F(2),var(1));
      end
    elseif nrow==3
      if isfunvar
        curlF=[-functionalDerivative(F(2),var(3))+functionalDerivative(F(3),var(2))
               +functionalDerivative(F(1),var(3))-functionalDerivative(F(3),var(1))
               -functionalDerivative(F(1),var(2))+functionalDerivative(F(2),var(1))];
      else
        curlF=[-diff(F(2),var(3))+diff(F(3),var(2))
               +diff(F(1),var(3))-diff(F(3),var(1))
               -diff(F(1),var(2))+diff(F(2),var(1))];
      end
    end
end

end