function Block=defineMatrixPattern(Simulation,Parameters,Mesh,Faces,Sizes)
% Define the row and column indices of the sparse blocks of the monolithic system
%
% The indices of the coupling blocks are padded with a dummy entry to make sure that the assembled
% blocks span the whole set of degrees of freedom of the two coupled discretizations.

Block=struct();

% Pattern of the diagonal blocks
for iD=1:Simulation.NumDiscretizations
  NumGlobalComp=Sizes(iD).NumGlobalComp;
  if strcmp(Parameters(iD).DiscretizationType,'CG')
    Block(iD,iD).RhsRowIndices=reshape((repmat(Mesh(iD).Elements',1,NumGlobalComp)...
      +reshape(repmat((0:NumGlobalComp-1)*Sizes(iD).NumNodes,...
                      Sizes(iD).NumElementNodes,1),1,[]))',[],1);
  else
    Block(iD,iD).RhsRowIndices=reshape(((repmat(reshape(Faces(iD,iD).Connectivity',[],1),...
      1,NumGlobalComp*Sizes(iD).NumFaceNodes)-1)*NumGlobalComp*Sizes(iD).NumFaceNodes...
      +(1:NumGlobalComp*Sizes(iD).NumFaceNodes))',[],1);
  end
  Block(iD,iD).LhsRowIndices=repelem(Block(iD,iD).RhsRowIndices,Sizes(iD).NumElementRhsCoef);
  Block(iD,iD).LhsColIndices=reshape(repmat(...
    reshape(Block(iD,iD).RhsRowIndices,[],Sizes(iD).NumElements),...
    Sizes(iD).NumElementRhsCoef,1),[],1);
end

% Pattern of the coupling blocks, assembled face by face along the interface
for iD=1:Simulation.NumDiscretizations
  for iDC=setdiff(1:Simulation.NumDiscretizations,iD)

    % Interface faces and connectivities, extracted for an effective use of parfor
    Interface=zeros(0,3);
    if matchField(Faces(iD,iDC),'Interface')
      Interface=Faces(iD,iDC).Interface;
    end
    NumFaces=size(Interface,1);
    CG =strcmp(Parameters(iD ).DiscretizationType,'CG');
    CGC=strcmp(Parameters(iDC).DiscretizationType,'CG');
    if CG;  Connec =Mesh(iD ).Elements'; else; Connec =Faces(iD ,iD ).Connectivity; end
    if CGC; ConnecC=Mesh(iDC).Elements'; else; ConnecC=Faces(iDC,iDC).Connectivity; end

    % Sizes of the elemental coupling block and strides of its rows and columns
    nRow=Sizes(iD ).NumElementRhsCoef;
    nCol=Sizes(iDC).NumElementRhsCoef;
    Stride =Sizes(iD ).NumGlobalComp*Sizes(iD ).NumFaceNodes;
    StrideC=Sizes(iDC).NumGlobalComp*Sizes(iDC).NumFaceNodes;
    Comp =(0:Sizes(iD ).NumGlobalComp-1)*Sizes(iD ).NumNodes;
    CompC=(0:Sizes(iDC).NumGlobalComp-1)*Sizes(iDC).NumNodes;

    % Degrees of freedom coupled by each interface face
    RowIndices=zeros(nRow*nCol,NumFaces);
    ColIndices=zeros(nRow*nCol,NumFaces);
    parfor iFace=1:NumFaces
      if CG
        indRow=reshape(Connec(Interface(iFace,1),:)'+Comp,[],1);
      else
        indRow=reshape((Connec(Interface(iFace,1),:)-1)*Stride+(1:Stride)',[],1);
      end
      if CGC
        indCol=reshape(ConnecC(Interface(iFace,3),:)'+CompC,[],1);
      else
        indCol=reshape((ConnecC(Interface(iFace,3),:)-1)*StrideC+(1:StrideC)',[],1);
      end
      RowIndices(:,iFace)=repelem(indRow,nCol);
      ColIndices(:,iFace)=repmat(indCol,nRow,1);
    end
    Block(iD,iDC).LhsRowIndices=[RowIndices(:);...
      Sizes(iD ).NumGlobalNodes*Sizes(iD ).NumGlobalComp];
    Block(iD,iDC).LhsColIndices=[ColIndices(:);...
      Sizes(iDC).NumGlobalNodes*Sizes(iDC).NumGlobalComp];
  end
end

end
