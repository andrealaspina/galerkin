function [Error]=computeError(...
         Time,Mesh,RefElement,Sizes,Solution,AnalyticSolution,Norm,Condition,Domain)
         % Compute error

% Loop in the elements
if strcmp(Norm,'Number')
  Error=computeNumberError(Time,Solution,AnalyticSolution,Domain);
else
  Error=0;
  parfor iElem=1:Sizes.NumElements
    ErrorElem=computeElementError(iElem,Mesh,Time,RefElement,Sizes,...
      Solution,AnalyticSolution,Norm,Condition,Domain);
    Error=Error+ErrorElem;
  end
  Error=sqrt(Error); % Remove this to obtain area/volume of the geometry
end

end

function [Error]=computeNumberError(...
  Time,Solution,AnalyticSolution,Domain)

% Get general parameters
if strcmp(Domain,'Time')
  t_or_w=Time.Time;
elseif strcmp(Domain,'Frequency')
  t_or_w=2*pi*Time.Frequency;
end

% Get solution
uh=Solution;
uA=AnalyticSolution(t_or_w);

% Compute (relative) error if all components of the analytical solution are not zero
if not(any(uA<1e-12)) % Relative error if all components of the analytical solution are not zero
  Error=abs((uh-uA)./uA);
else                  % Aboslute error if some components of the analytical solution are zero
  Error=abs(uh-uA);
end

end

function [ErrorElem]=computeElementError(...
  iElem,Mesh,Time,RefElement,Sizes,Solution,AnalyticSolution,Norm,Condition,Domain)

% Get general parameters
nsd=Sizes.NumSpaceDim;
Ce=Mesh.Elements(:,iElem)';
Xe=Mesh.Nodes(:,Ce)';
uA=AnalyticSolution;
if strcmp(Domain,'Time')
  t_or_w=Time.Time;
elseif strcmp(Domain,'Frequency')
  t_or_w=2*pi*Time.Frequency;
end

% Get solution
uhe=Solution(Ce,:);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions(1,RefElement,RefElement,Xe,nsd);

% Compute variables at Gauss points
Xeg=Ne*Xe;
uheg=Ne*uhe;
uAeg=uA(Xeg(:,1),Xeg(:,2),Xeg(:,3),t_or_w);

% Compute L2 contribution
ErrorElemL2=sum(((uheg-uAeg).^2).*weg);

% Compute elemental contribution
switch Norm
  case 'L2'
    ErrorElem=ErrorElemL2;
  case 'Hdiv'
    % Check variable size
    if size(uhe,2)~=nsd
      error('Impossible to compute error in the Hdiv norm for this variable');
    end
    
    % Get analytic solution at nodes (potential discontinous solution across elements)
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
    % Check variable size
    if size(uhe,2)~=nsd
      error('Impossible to compute error in the Hcurl norm for this variable');
    end
    
    % Get analytic solution at nodes (potential discontinous solution across elements)
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

% Delete error contribution if the condition is not met
if not(feval(str2func(Condition),mean(Xe(:,1)),mean(Xe(:,2)),mean(Xe(:,3))))
  ErrorElem=0;
end

end