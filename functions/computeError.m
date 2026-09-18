function [Error]=computeError(Time,Mesh,Parameters,RefElement,Sizes,Solution,AnalyticSolution,...
  Norm,Condition,Domain)
% Compute the error of a solution with respect to an analytic one in the given norm

if strcmp(Norm,'Number')
  Error=computeNumberError(Time,Solution,AnalyticSolution,Domain);
else
  Error=0;
  parfor iElem=1:Sizes.NumElements
    Error=Error+computeElementError(iElem,Mesh,Time,Parameters,RefElement,Sizes,...
      Solution,AnalyticSolution,Norm,Condition,Domain);
  end
  Error=sqrt(Error); % Remove this to obtain the area/volume of the geometry
end

end

function [Error]=computeNumberError(Time,Solution,AnalyticSolution,Domain)
% Error of a scalar quantity (relative if no analytic component is zero, absolute otherwise)

uh=Solution;
uA=AnalyticSolution(getIndependentVariable(Time,Domain));
if not(any(uA<1e-12))
  Error=abs((uh-uA)./uA);
else
  Error=abs(uh-uA);
end

end

function [ErrorElem]=computeElementError(iElem,Mesh,Time,Parameters,RefElement,Sizes,...
  Solution,AnalyticSolution,Norm,Condition,Domain)
% Elemental contribution to the squared error

% Get general parameters
nsd=Sizes.NumSpaceDim;
k=Parameters.Degree;
Ce=Mesh.Elements(:,iElem)';
Xe=Mesh.Nodes(:,Ce)';
uA=AnalyticSolution;
t_or_w_or_n=getIndependentVariable(Time,Domain);

% Get solution
if k==0
  uhe=Solution(iElem,:);
else
  uhe=Solution(Ce,:);
end

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions(1,RefElement,RefElement,Xe,nsd);

% Compute variables at Gauss points
if k==0
  Xeg=mean(Xe);
else
  Xeg=Ne*Xe;
end
uheg=Ne*uhe;
uAeg=uA(Xeg(:,1),Xeg(:,2),Xeg(:,3),t_or_w_or_n);

% Compute L2 contribution
ErrorElemL2=sum(((uheg-uAeg).^2).*weg);

% Compute elemental contribution
switch Norm
  case 'L2'
    ErrorElem=ErrorElemL2;
  case 'Hdiv'
    if size(uhe,2)~=nsd
      error('Impossible to compute error in the Hdiv norm for this variable');
    end

    % Get analytic solution at nodes (potentially discontinuous across elements)
    uAe=Ne\uAeg;

    % Compute Hdiv contribution
    if nsd==2
      ErrorElemHdiv=sum((((Nex*uhe(:,1)+Ney*uhe(:,2))-...
                          (Nex*uAe(:,1)+Ney*uAe(:,2))).^2).*weg);
    elseif nsd==3
      ErrorElemHdiv=sum((((Nex*uhe(:,1)+Ney*uhe(:,2)+Nez*uhe(:,3))-...
                          (Nex*uAe(:,1)+Ney*uAe(:,2)+Nez*uAe(:,3))).^2).*weg);
    end
    ErrorElem=[ErrorElemL2,ErrorElemHdiv];
  case 'Hcurl'
    if size(uhe,2)~=nsd
      error('Impossible to compute error in the Hcurl norm for this variable');
    end

    % Get analytic solution at nodes (potentially discontinuous across elements)
    uAe=Ne\uAeg;

    % Compute Hcurl contribution
    if nsd==2
      ErrorElemHcurl=sum((((Nex*uhe(:,2)-Ney*uhe(:,1))-...
                           (Nex*uAe(:,2)-Ney*uAe(:,1))).^2).*weg);
    elseif nsd==3
      ErrorElemHcurl=sum(([((Ney*uhe(:,3)-Nez*uhe(:,2))-...
                            (Ney*uAe(:,3)-Nez*uAe(:,2))).^2,...
                           ((Nez*uhe(:,1)-Nex*uhe(:,3))-...
                            (Nez*uAe(:,1)-Nex*uAe(:,3))).^2,...
                           ((Nex*uhe(:,2)-Ney*uhe(:,1))-...
                            (Nex*uAe(:,2)-Ney*uAe(:,1))).^2]).*weg);
    end
    ErrorElem=[ErrorElemL2,ErrorElemHcurl];
end
ErrorElem=abs(ErrorElem);

% Discard the contribution if the condition is not met at the element center
if not(feval(str2func(Condition),mean(Xe(:,1)),mean(Xe(:,2)),mean(Xe(:,3))))
  ErrorElem=0;
end

end

function [t_or_w_or_n]=getIndependentVariable(Time,Domain)
% Time, angular frequency or mode number, depending on the domain

if strcmp(Domain,'Time')
  t_or_w_or_n=Time.Time;
elseif strcmp(Domain,'Frequency')
  t_or_w_or_n=2*pi*Time.Frequency;
elseif strcmp(Domain,'Mode')
  t_or_w_or_n=Time.Mode(end);
end

end
