function Error=computeError(Time,Mesh,Parameters,RefElement,Sizes,Solution,AnalyticSolution,...
  Norm,Condition,Domain)
% Compute the error of a discrete solution against its analytic counterpart

if strcmp(Norm,'Number')
  Error=computeNumberError(Time,Solution,AnalyticSolution,Domain);
else
  Error=0;
  parfor iElem=1:Sizes.NumElements
    Error=Error+computeElementError(iElem,Mesh,Time,Parameters,RefElement,Sizes,...
      Solution,AnalyticSolution,Norm,Condition,Domain);
  end
  Error=sqrt(Error); % Drop the square root to obtain the area/volume of the geometry
end

end

function Error=computeNumberError(Time,Solution,AnalyticSolution,Domain)
% Compute the error of a scalar quantity of interest

% Time, angular frequency or mode number, depending on the domain
if strcmp(Domain,'Time')
  t=Time.Time;
elseif strcmp(Domain,'Frequency')
  t=2*pi*Time.Frequency;
elseif strcmp(Domain,'Mode')
  t=Time.Mode(end);
end

% Relative error, or absolute error if some component of the analytic solution vanishes
uA=AnalyticSolution(t);
if not(any(uA<1e-12))
  Error=abs((Solution-uA)./uA);
else
  Error=abs(Solution-uA);
end

end

function Error=computeElementError(iElem,Mesh,Time,Parameters,RefElement,Sizes,Solution,...
  AnalyticSolution,Norm,Condition,Domain)
% Compute the contribution of one element to the squared error

% Element data
nsd=Sizes.NumSpaceDim;
k=Parameters.Degree;
Ce=Mesh.Elements(:,iElem)';
Xe=Mesh.Nodes(:,Ce)';
uA=AnalyticSolution;

% Time, angular frequency or mode number, depending on the domain
if strcmp(Domain,'Time')
  t=Time.Time;
elseif strcmp(Domain,'Frequency')
  t=2*pi*Time.Frequency;
elseif strcmp(Domain,'Mode')
  t=Time.Mode(end);
end

% Discrete solution, nodal for k>0 and cellwise constant for k=0
if k==0
  uhe=Solution(iElem,:);
else
  uhe=Solution(Ce,:);
end

% Shape functions and weights at the Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions(1,RefElement,RefElement,Xe,nsd);

% Discrete and analytic solutions at the Gauss points
if k==0
  Xeg=mean(Xe);
else
  Xeg=Ne*Xe;
end
uheg=Ne*uhe;
uAeg=uA(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);

% Contribution in the L2 norm
Error=sum(((uheg-uAeg).^2).*weg);

% Contribution of the semi-norm, on the analytic solution interpolated at the element nodes
if not(strcmp(Norm,'L2'))
  if size(uhe,2)~=nsd
    error('Impossible to compute the error in the %s norm for this variable',Norm);
  end
  uAe=Ne\uAeg;
  D={Nex,Ney,Nez};
  switch Norm
    case 'Hdiv'
      Divh=0; DivA=0;
      for iSD=1:nsd
        Divh=Divh+D{iSD}*uhe(:,iSD);
        DivA=DivA+D{iSD}*uAe(:,iSD);
      end
      Error=[Error,sum(((Divh-DivA).^2).*weg)];
    case 'Hcurl'
      if nsd==2
        Curlh=Nex*uhe(:,2)-Ney*uhe(:,1);
        CurlA=Nex*uAe(:,2)-Ney*uAe(:,1);
      elseif nsd==3
        Curlh=[Ney*uhe(:,3)-Nez*uhe(:,2),Nez*uhe(:,1)-Nex*uhe(:,3),Nex*uhe(:,2)-Ney*uhe(:,1)];
        CurlA=[Ney*uAe(:,3)-Nez*uAe(:,2),Nez*uAe(:,1)-Nex*uAe(:,3),Nex*uAe(:,2)-Ney*uAe(:,1)];
      end
      Error=[Error,sum(((Curlh-CurlA).^2).*weg)];
    otherwise
      error('Norm not yet implemented');
  end
end
Error=abs(Error);

% Drop the contribution of the elements that do not meet the condition
if not(feval(str2func(Condition),mean(Xe(:,1)),mean(Xe(:,2)),mean(Xe(:,3))))
  Error=0;
end

end
