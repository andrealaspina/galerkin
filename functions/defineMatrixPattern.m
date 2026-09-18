function [Block]=defineMatrixPattern(Simulation,Parameters,Mesh,Faces,Sizes)
% Define the sparsity pattern (row and column indices) of the block matrices and vectors

Block=struct();

% Pattern of the single discretizations
for iD=1:Simulation.NumDiscretizations
  NumNodes=Sizes(iD).NumNodes;
  NumElements=Sizes(iD).NumElements;
  NumElementNodes=Sizes(iD).NumElementNodes;
  NumElementFaces=Sizes(iD).NumElementFaces;
  NumFaceNodes=Sizes(iD).NumFaceNodes;
  NumGlobalComp=Sizes(iD).NumGlobalComp;
  if strcmp(Parameters(iD).DiscretizationType,'CG')
    NumElementDofs=NumGlobalComp*NumElementNodes;
    Block(iD,iD).RhsRowIndices=reshape((repmat(Mesh(iD).Elements',1,NumGlobalComp)...
      +reshape(repmat((0:NumGlobalComp-1)*NumNodes,NumElementNodes,1),1,[]))',[],1);
  else
    NumElementDofs=NumGlobalComp*NumElementFaces*NumFaceNodes;
    Block(iD,iD).RhsRowIndices=reshape(((repmat(reshape(Faces(iD,iD).Connectivity',[],1),...
      1,NumGlobalComp*NumFaceNodes)-1)*NumGlobalComp*NumFaceNodes...
      +(1:NumGlobalComp*NumFaceNodes))',[],1);
  end
  Block(iD,iD).LhsRowIndices=repelem(Block(iD,iD).RhsRowIndices,NumElementDofs);
  Block(iD,iD).LhsColIndices=reshape(repmat(reshape(Block(iD,iD).RhsRowIndices,...
    NumElementDofs,NumElements),NumElementDofs,1),[],1);
end

% Pattern of the coupled discretizations
for iD=1:Simulation.NumDiscretizations
  for iDC=setdiff(1:Simulation.NumDiscretizations,iD)
    if matchField(Faces(iD,iDC),'Interface')
      Interface=Faces(iD,iDC).Interface;
    else
      Interface=zeros(0,5);
    end
    NumFacesInterface=size(Interface,1);
    NumRows=getNumElementDofs(iD,Parameters,Sizes);
    NumCols=getNumElementDofs(iDC,Parameters,Sizes);
    RowIndices=zeros(NumRows*NumCols,NumFacesInterface);
    ColIndices=zeros(NumRows*NumCols,NumFacesInterface);
    parfor iFaceInterface=1:NumFacesInterface
      DofsRow=getElementDofs(Interface(iFaceInterface,1),iD,Parameters,Mesh,Faces,Sizes);
      DofsCol=getElementDofs(Interface(iFaceInterface,3),iDC,Parameters,Mesh,Faces,Sizes);
      RowIndices(:,iFaceInterface)=repelem(DofsRow,NumCols);
      ColIndices(:,iFaceInterface)=repmat(DofsCol,NumRows,1);
    end

    % Add a dummy entry in the last row and column to span all the DOFs of both discretizations
    Block(iD,iDC).LhsRowIndices=[RowIndices(:);Sizes(iD).NumGlobalNodes*Sizes(iD).NumGlobalComp];
    Block(iD,iDC).LhsColIndices=[ColIndices(:);Sizes(iDC).NumGlobalNodes*Sizes(iDC).NumGlobalComp];
  end
end

end

function [NumElementDofs]=getNumElementDofs(iD,Parameters,Sizes)
% Number of DOFs of an element of discretization iD

if strcmp(Parameters(iD).DiscretizationType,'CG')
  NumElementDofs=Sizes(iD).NumGlobalComp*Sizes(iD).NumElementNodes;
elseif strcmp(Parameters(iD).DiscretizationType,'HDG')
  NumElementDofs=Sizes(iD).NumGlobalComp*Sizes(iD).NumElementFaces*Sizes(iD).NumFaceNodes;
end

end

function [Dofs]=getElementDofs(iElem,iD,Parameters,Mesh,Faces,Sizes)
% Global DOFs of element iElem of discretization iD (component-major for CG, face-major for HDG)

if strcmp(Parameters(iD).DiscretizationType,'CG')
  Nodes=Mesh(iD).Elements(:,iElem);
  Dofs=reshape(Nodes+(0:Sizes(iD).NumGlobalComp-1)*Sizes(iD).NumNodes,[],1);
elseif strcmp(Parameters(iD).DiscretizationType,'HDG')
  NumFaceDofs=Sizes(iD).NumGlobalComp*Sizes(iD).NumFaceNodes;
  ElementFaces=Faces(iD,iD).Connectivity(iElem,:)';
  Dofs=reshape(((ElementFaces-1)*NumFaceDofs+(1:NumFaceDofs))',[],1);
end

end
