function [Block]=...
         defineMatrixPattern(Simulation,Parameters,Mesh,Faces,Sizes)
         % Define matrix pattern

% Initialize structure
Block=struct();

% Define matrix pattern of the single discretizations
for iD=1:Simulation.NumDiscretizations
  % Get info
  NumNodes=Sizes(iD).NumNodes;
  NumElements=Sizes(iD).NumElements;
  NumElementNodes=Sizes(iD).NumElementNodes;
  NumElementFaces=Sizes(iD).NumElementFaces;
  NumFaceNodes=Sizes(iD).NumFaceNodes;
  NumGlobalComp=Sizes(iD).NumGlobalComp;
  if strcmp(Parameters(iD).DiscretizationType,'CG')
    Block(iD,iD).RhsRowIndices=reshape((repmat(Mesh(iD).Elements',1,NumGlobalComp)...
      +reshape(repmat((0:NumGlobalComp-1)*NumNodes,NumElementNodes,1),1,[]))',[],1);
    Block(iD,iD).LhsRowIndices=repelem(Block(iD,iD).RhsRowIndices,NumGlobalComp*NumElementNodes);
    Block(iD,iD).LhsColIndices=reshape(repmat(reshape(Block(iD,iD).RhsRowIndices,...
      NumGlobalComp*NumElementNodes,NumElements),NumGlobalComp*NumElementNodes,1),[],1);
  else
    Block(iD,iD).RhsRowIndices=reshape(((repmat(reshape(Faces(iD,iD).Connectivity',[],1),...
      1,NumGlobalComp*NumFaceNodes)-1)*NumGlobalComp*NumFaceNodes...
      +(1:NumGlobalComp*NumFaceNodes))',[],1);
    Block(iD,iD).LhsRowIndices=repelem(Block(iD,iD).RhsRowIndices,...
      NumGlobalComp*NumElementFaces*NumFaceNodes);
    Block(iD,iD).LhsColIndices=reshape(repmat(reshape(Block(iD,iD).RhsRowIndices,[],NumElements),...
      NumGlobalComp*NumElementFaces*NumFaceNodes,1),[],1);
  end
end

% Define matrix pattern of the coupled discretizations
for iD=1:Simulation.NumDiscretizations
  iDOthers=setdiff(1:Simulation.NumDiscretizations,iD);
  for iDAux=1:length(iDOthers)
    iDC=iDOthers(iDAux);
    % Get info
    if matchField(Faces(iD,iDC),'Interface')
      NumFacesInterface=size(Faces(iD,iDC).Interface,1);
    else
      NumFacesInterface=0;
    end
    NumNodes=Sizes(iD).NumNodes;
    NumElementNodes=Sizes(iD).NumElementNodes;
    NumElementFaces=Sizes(iD).NumElementFaces;
    NumFaceNodes=Sizes(iD).NumFaceNodes;
    NumGlobalComp=Sizes(iD).NumGlobalComp;
    NumNodesC=Sizes(iDC).NumNodes;
    NumElementNodesC=Sizes(iDC).NumElementNodes;
    NumElementFacesC=Sizes(iDC).NumElementFaces;
    NumFaceNodesC=Sizes(iDC).NumFaceNodes;
    NumGlobalCompC=Sizes(iDC).NumGlobalComp;
    if strcmp(Parameters(iD).DiscretizationType,'CG')
      nRow=NumGlobalComp*NumElementNodes;
    elseif strcmp(Parameters(iD).DiscretizationType,'HDG')
      nRow=NumGlobalComp*NumElementFaces*NumFaceNodes;
    end
    if strcmp(Parameters(iDC).DiscretizationType,'CG')
      nCol=NumGlobalCompC*NumElementNodesC;
    elseif strcmp(Parameters(iDC).DiscretizationType,'HDG')
      nCol=NumGlobalCompC*NumElementFacesC*NumFaceNodesC;
    end
    nTot=nRow*nCol;
    ind_i=zeros(nTot,NumFacesInterface);
    ind_j=zeros(nTot,NumFacesInterface);
    parfor iFaceInterface=1:NumFacesInterface
      ParametersPar=Parameters; MeshPar=Mesh; FacesPar=Faces;
      indRowC=[];
      indColC=[];
      if strcmp(ParametersPar(iD).DiscretizationType,'CG') && ...
         strcmp(ParametersPar(iDC).DiscretizationType,'CG')
        C=MeshPar(iD).Elements';
        iElem=FacesPar(iD,iDC).Interface(iFaceInterface,1);
        Ce=C(iElem,:);
        indRowC=zeros(NumGlobalComp*NumElementNodes,1);
        aux=(1:NumElementNodes);
        for iComp=1:NumGlobalComp
          indRowC((iComp-1)*NumElementNodes+aux)=...
               Ce+(iComp-1)*NumNodes;
        end
        CC=MeshPar(iDC).Elements';
        iElemC=FacesPar(iD,iDC).Interface(iFaceInterface,3);
        CCe=CC(iElemC,:);
        indColC=zeros(NumGlobalCompC*NumElementNodesC,1);
        auxC=(1:NumElementNodesC);
        for iCompC=1:NumGlobalCompC
          indColC((iCompC-1)*NumElementNodesC+auxC)=...
              CCe+(iCompC-1)*NumNodesC;
        end
      elseif strcmp(ParametersPar(iD).DiscretizationType,'CG') && ...
         strcmp(ParametersPar(iDC).DiscretizationType,'HDG')
        C=MeshPar(iD).Elements';
        FC=FacesPar(iDC,iDC).Connectivity;
        iElem=FacesPar(iD,iDC).Interface(iFaceInterface,1);
        Ce=C(iElem,:);
        indRowC=zeros(NumGlobalComp*NumElementNodes,1);
        aux=(1:NumElementNodes);
        for iComp=1:NumGlobalComp
          indRowC((iComp-1)*NumElementNodes+aux)=...
              Ce+(iComp-1)*NumNodes;
        end
        iElemC=FacesPar(iD,iDC).Interface(iFaceInterface,3);
        FCe=FC(iElemC,:);
        indColC=zeros(NumGlobalCompC*NumElementFacesC*NumFaceNodesC,1);
        auxC=(1:NumGlobalCompC*NumFaceNodesC);
        for iElemFaceC=1:NumElementFaces
          indColC((iElemFaceC-1)*NumGlobalCompC*NumFaceNodesC+auxC)=...
             (FCe(iElemFaceC)-1)*NumGlobalCompC*NumFaceNodesC+auxC;
        end
      elseif strcmp(ParametersPar(iD).DiscretizationType,'HDG') && ...
             strcmp(ParametersPar(iDC).DiscretizationType,'CG')
        F=FacesPar(iD,iD).Connectivity;
        CC=MeshPar(iDC).Elements';
        iElem=FacesPar(iD,iDC).Interface(iFaceInterface,1);
        Fe=F(iElem,:);
        indRowC=zeros(NumGlobalComp*NumElementFaces*NumFaceNodes,1);
        aux=(1:NumGlobalComp*NumFaceNodes);
        for iElemFace=1:NumElementFaces
          indRowC((iElemFace-1)*NumGlobalComp*NumFaceNodes+aux)=...
              (Fe(iElemFace)-1)*NumGlobalComp*NumFaceNodes+aux;
        end
        iElemC=FacesPar(iD,iDC).Interface(iFaceInterface,3);
        CCe=CC(iElemC,:);
        indColC=zeros(NumGlobalCompC*NumElementNodesC,1);
        auxC=(1:NumElementNodesC);
        for iCompC=1:NumGlobalCompC
          indColC((iCompC-1)*NumElementNodesC+auxC)=...
            CCe+(iCompC-1)*NumNodesC;
        end
      end
      ind_ie=zeros(nTot,1);
      ind_je=zeros(nTot,1);
      iN=0;
      for iRow=1:nRow
        for iCol=1:nCol
          iN=iN+1;
          ind_ie(iN)=indRowC(iRow);
          ind_je(iN)=indColC(iCol);
        end
      end
      ind_i(:,iFaceInterface)=reshape(ind_ie',[],1);
      ind_j(:,iFaceInterface)=reshape(ind_je',[],1);
    end
    Block(iD,iDC).LhsRowIndices=reshape(ind_i,[],1);
    Block(iD,iDC).LhsColIndices=reshape(ind_j,[],1);
    
    % Add dummy indices to make sure to span the entire DOFs of the coupled discretizations    
    if strcmp(Parameters(iD).DiscretizationType,'CG') && ...
       strcmp(Parameters(iDC).DiscretizationType,'CG')
      Block(iD,iDC).LhsRowIndices=[Block(iD,iDC).LhsRowIndices;...
                                   NumGlobalComp*NumNodes];
      Block(iD,iDC).LhsColIndices=[Block(iD,iDC).LhsColIndices;...
                                   NumGlobalCompC*NumNodesC];
    elseif strcmp(Parameters(iD).DiscretizationType,'CG') && ...
           strcmp(Parameters(iDC).DiscretizationType,'HDG')
      NumFacesC=Sizes(iDC).NumFaces;
      Block(iD,iDC).LhsRowIndices=[Block(iD,iDC).LhsRowIndices;...
                                   NumGlobalComp*NumNodes];
      Block(iD,iDC).LhsColIndices=[Block(iD,iDC).LhsColIndices;...
                                   NumGlobalCompC*NumFacesC*NumFaceNodesC];
    elseif strcmp(Parameters(iD).DiscretizationType,'HDG') && ...
           strcmp(Parameters(iDC).DiscretizationType,'CG')
      NumFaces=Sizes(iD).NumFaces;
      Block(iD,iDC).LhsRowIndices=[Block(iD,iDC).LhsRowIndices;...
                                   NumGlobalComp*NumFaces*NumFaceNodes];
      Block(iD,iDC).LhsColIndices=[Block(iD,iDC).LhsColIndices;...
                                   NumGlobalCompC*NumNodesC];
    elseif strcmp(Parameters(iD).DiscretizationType,'HDG') && ...
           strcmp(Parameters(iDC).DiscretizationType,'HDG')
      NumFaces=Sizes(iD).NumFaces;
      NumFacesC=Sizes(iDC).NumFaces;
      Block(iD,iDC).LhsRowIndices=[Block(iD,iDC).LhsRowIndices;...
                                   NumGlobalComp*NumFaces*NumFaceNodes];
      Block(iD,iDC).LhsColIndices=[Block(iD,iDC).LhsColIndices;...
                                   NumGlobalCompC*NumFacesC*NumFaceNodesC];
    end
  end
end

end