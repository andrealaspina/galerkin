classdef ElasticityModal_CG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) NumSpaceDim;

    % Discretization type
    DiscretizationType='CG';

    % Time derivative order
    TimeDerOrder=NaN;
    
    % Domain
    Domain='Mode';
    
  end
  
  methods
    
    %% Initialize unknowns
    function [Block]=initializeUnknowns(~,iD,Block,~,~,Sizes)
      Block(iD,iD).SolutionGlobal=zeros(Sizes(iD).NumGlobalNodes,Sizes(iD).NumGlobalComp,1);
    end

    %% Evaluate solution at fixed DOFs
    function [Block]=evaluateSolutionFixedDofs(~,iD,Block,Parameters,Mesh,Time,Sizes)
      Block(iD,iD).SolutionGlobal(Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp),:)=...
        Parameters(iD).Displacement(...
        Mesh(iD).Nodes(1,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',...
        Mesh(iD).Nodes(2,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',...
        Mesh(iD).Nodes(3,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',min(Time.Mode));
    end
    
    %% Build block
    function [Block,Elements]=buildBlock(~,iD,Block,Elements,~,Parameters,~,~,~,RefElement,Sizes)
      NodesElem=Elements(iD).Nodes;
      LhsCoef=zeros(Sizes(iD).NumElementLhsCoef,Sizes(iD).NumElements);
      RhsCoef=zeros(Sizes(iD).NumElementLhsCoef,Sizes(iD).NumElements);
      parfor iElem=1:Sizes(iD).NumElements
        [LhsGlobalElem,RhsGlobalElem]=...
          buildBlockElement(NodesElem{iElem},...
          Parameters,RefElement,Sizes);
        LhsCoef(:,iElem)=reshape(LhsGlobalElem',[],1);
        RhsCoef(:,iElem)=reshape(RhsGlobalElem',[],1);
      end
      Block(iD,iD).LhsGlobal=sparse(Block(iD,iD).LhsRowIndices,...
                                    Block(iD,iD).LhsColIndices,LhsCoef(:));
      Block(iD,iD).RhsGlobal=sparse(Block(iD,iD).LhsRowIndices,...
                                    Block(iD,iD).LhsColIndices,RhsCoef(:));
    end
    
    %% Store results
    function [Results]=storeResults(~,iD,~,Results,Block,~,~,~,Time,~,~)
      Results(iD).Mode=Time.Mode;
      Results(iD).Displacement=Block(iD,iD).SolutionGlobal;
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal]=buildBlockElement(...
  Nodes,Parameters,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
rho=Parameters.Density;
Xe=Nodes';
Xem=sum(Xe(1:NumElementFaces,:),1)/NumElementFaces;

% Initialize lhs
Kuu=zeros(nsd*NumElementNodes,nsd*NumElementNodes);

% Initialize rhs
Muu=zeros(nsd*NumElementNodes,nsd*NumElementNodes);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,~,~]=mapShapeFunctions(1,RefElement,RefElement,Xe,nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute basic matrices
NweT=(weg.*Ne)';
NwexT=(weg.*Nex)';
NweyT=(weg.*Ney)';
if nsd==3
  NwezT=(weg.*Nez)';
end
Me=NweT*Ne;
Me=(Me+Me')/2;

% Compute stress and its linearization
if nsd==2
  [~,dsdFeg]=...
    computeStress('no','yes','no',Parameters,Xem,...
     repmat([1,0,0,1],numel(weg),1),Sizes);
  dsxxdFxxeg=dsdFeg(:,1);  dsxxdFxyeg=dsdFeg(:,2);  dsxxdFyxeg=dsdFeg(:,3);  dsxxdFyyeg=dsdFeg(:,4);
                           dsxydFxyeg=dsdFeg(:,6);  dsxydFyxeg=dsdFeg(:,7);  dsxydFyyeg=dsdFeg(:,8);
                                                    dsyxdFyxeg=dsdFeg(:,11); dsyxdFyyeg=dsdFeg(:,12);
                                                                             dsyydFyyeg=dsdFeg(:,16);

elseif nsd==3
  [~,dsdFeg]=...
    computeStress('no','yes','no',Parameters,Xem,...
    repmat([1,0,0,0,1,0,0,0,1],numel(weg),1),Sizes);
  dsxxdFxxeg=dsdFeg(:,1);  dsxxdFxyeg=dsdFeg(:,2);  dsxxdFxzeg=dsdFeg(:,3);  dsxxdFyxeg=dsdFeg(:,4);  dsxxdFyyeg=dsdFeg(:,5);  dsxxdFyzeg=dsdFeg(:,6);  dsxxdFzxeg=dsdFeg(:,7);  dsxxdFzyeg=dsdFeg(:,8);  dsxxdFzzeg=dsdFeg(:,9);
                           dsxydFxyeg=dsdFeg(:,11); dsxydFxzeg=dsdFeg(:,12); dsxydFyxeg=dsdFeg(:,13); dsxydFyyeg=dsdFeg(:,14); dsxydFyzeg=dsdFeg(:,15); dsxydFzxeg=dsdFeg(:,16); dsxydFzyeg=dsdFeg(:,17); dsxydFzzeg=dsdFeg(:,18);
                                                    dsxzdFxzeg=dsdFeg(:,21); dsxzdFyxeg=dsdFeg(:,22); dsxzdFyyeg=dsdFeg(:,23); dsxzdFyzeg=dsdFeg(:,24); dsxzdFzxeg=dsdFeg(:,25); dsxzdFzyeg=dsdFeg(:,26); dsxzdFzzeg=dsdFeg(:,27);
                                                                             dsyxdFyxeg=dsdFeg(:,31); dsyxdFyyeg=dsdFeg(:,32); dsyxdFyzeg=dsdFeg(:,33); dsyxdFzxeg=dsdFeg(:,34); dsyxdFzyeg=dsdFeg(:,35); dsyxdFzzeg=dsdFeg(:,36);
                                                                                                      dsyydFyyeg=dsdFeg(:,41); dsyydFyzeg=dsdFeg(:,42); dsyydFzxeg=dsdFeg(:,43); dsyydFzyeg=dsdFeg(:,44); dsyydFzzeg=dsdFeg(:,45);
                                                                                                                               dsyzdFyzeg=dsdFeg(:,51); dsyzdFzxeg=dsdFeg(:,52); dsyzdFzyeg=dsdFeg(:,53); dsyzdFzzeg=dsdFeg(:,54);
                                                                                                                                                        dszxdFzxeg=dsdFeg(:,61); dszxdFzyeg=dsdFeg(:,62); dszxdFzzeg=dsdFeg(:,63);
                                                                                                                                                                                 dszydFzyeg=dsdFeg(:,71); dszydFzzeg=dsdFeg(:,72);
                                                                                                                                                                                                          dszzdFzzeg=dsdFeg(:,81);
end

% Compute lhs
if nsd==2
  Kuu(ne1,ne1)=+NwexT*(dsxxdFxxeg.*Nex+dsxxdFxyeg.*Ney)...
               +NweyT*(dsxxdFxyeg.*Nex+dsxydFxyeg.*Ney);
  Kuu(ne1,ne2)=+NwexT*(dsxxdFyxeg.*Nex+dsxxdFyyeg.*Ney)...
               +NweyT*(dsxydFyxeg.*Nex+dsxydFyyeg.*Ney);
  Kuu(ne2,ne1)=+NwexT*(dsxxdFyxeg.*Nex+dsxydFyxeg.*Ney)...
               +NweyT*(dsxxdFyyeg.*Nex+dsxydFyyeg.*Ney);
  Kuu(ne2,ne2)=+NwexT*(dsyxdFyxeg.*Nex+dsyxdFyyeg.*Ney)...
               +NweyT*(dsyxdFyyeg.*Nex+dsyydFyyeg.*Ney);
elseif nsd==3
  Kuu(ne1,ne1)=+NwexT*(dsxxdFxxeg.*Nex+dsxxdFxyeg.*Ney+dsxxdFxzeg.*Nez)...
               +NweyT*(dsxxdFxyeg.*Nex+dsxydFxyeg.*Ney+dsxydFxzeg.*Nez)...
               +NwezT*(dsxxdFxzeg.*Nex+dsxydFxzeg.*Ney+dsxzdFxzeg.*Nez);
  Kuu(ne1,ne2)=+NwexT*(dsxxdFyxeg.*Nex+dsxxdFyyeg.*Ney+dsxxdFyzeg.*Nez)...
               +NweyT*(dsxydFyxeg.*Nex+dsxydFyyeg.*Ney+dsxydFyzeg.*Nez)...
               +NwezT*(dsxzdFyxeg.*Nex+dsxzdFyyeg.*Ney+dsxzdFyzeg.*Nez);
  Kuu(ne1,ne3)=+NwexT*(dsxxdFzxeg.*Nex+dsxxdFzyeg.*Ney+dsxxdFzzeg.*Nez)...
               +NweyT*(dsxydFzxeg.*Nex+dsxydFzyeg.*Ney+dsxydFzzeg.*Nez)...
               +NwezT*(dsxzdFzxeg.*Nex+dsxzdFzyeg.*Ney+dsxzdFzzeg.*Nez);
  Kuu(ne2,ne1)=+NwexT*(dsxxdFyxeg.*Nex+dsxydFyxeg.*Ney+dsxzdFyxeg.*Nez)...
               +NweyT*(dsxxdFyyeg.*Nex+dsxydFyyeg.*Ney+dsxzdFyyeg.*Nez)...
               +NwezT*(dsxxdFyzeg.*Nex+dsxydFyzeg.*Ney+dsxzdFyzeg.*Nez);
  Kuu(ne2,ne2)=+NwexT*(dsyxdFyxeg.*Nex+dsyxdFyyeg.*Ney+dsyxdFyzeg.*Nez)...
               +NweyT*(dsyxdFyyeg.*Nex+dsyydFyyeg.*Ney+dsyydFyzeg.*Nez)...
               +NwezT*(dsyxdFyzeg.*Nex+dsyydFyzeg.*Ney+dsyzdFyzeg.*Nez);
  Kuu(ne2,ne3)=+NwexT*(dsyxdFzxeg.*Nex+dsyxdFzyeg.*Ney+dsyxdFzzeg.*Nez)...
               +NweyT*(dsyydFzxeg.*Nex+dsyydFzyeg.*Ney+dsyydFzzeg.*Nez)...
               +NwezT*(dsyzdFzxeg.*Nex+dsyzdFzyeg.*Ney+dsyzdFzzeg.*Nez);
  Kuu(ne3,ne1)=+NwexT*(dsxxdFzxeg.*Nex+dsxydFzxeg.*Ney+dsxzdFzxeg.*Nez)...
               +NweyT*(dsxxdFzyeg.*Nex+dsxydFzyeg.*Ney+dsxzdFzyeg.*Nez)...
               +NwezT*(dsxxdFzzeg.*Nex+dsxydFzzeg.*Ney+dsxzdFzzeg.*Nez);
  Kuu(ne3,ne2)=+NwexT*(dsyxdFzxeg.*Nex+dsyydFzxeg.*Ney+dsyzdFzxeg.*Nez)...
               +NweyT*(dsyxdFzyeg.*Nex+dsyydFzyeg.*Ney+dsyzdFzyeg.*Nez)...
               +NwezT*(dsyxdFzzeg.*Nex+dsyydFzzeg.*Ney+dsyzdFzzeg.*Nez);
  Kuu(ne3,ne3)=+NwexT*(dszxdFzxeg.*Nex+dszxdFzyeg.*Ney+dszxdFzzeg.*Nez)...
               +NweyT*(dszxdFzyeg.*Nex+dszydFzyeg.*Ney+dszydFzzeg.*Nez)...
               +NwezT*(dszxdFzzeg.*Nex+dszydFzzeg.*Ney+dszzdFzzeg.*Nez);
end

Muu(ne1,ne1)=rho*Me;
Muu(ne2,ne2)=rho*Me;
if nsd==3
  Muu(ne3,ne3)=rho*Me;
end

% Compute elemental contributions to lhs and rhs

% Lhs for global problem
LhsGlobal=(Kuu+Kuu.')/2;

% Rhs for global problem
RhsGlobal=(Muu+Muu.')/2;

end