%% Pre-processing

% Store input
Input=cell2struct(cellfun(@eval,who('-regexp','.*'),'UniformOutput',false),who('-regexp','.*'));

% Get file name
FileName=dbstack(1).name;

% Start diary
if matchField(Options,'SaveResults','yes')
  if not(matchField(Options,'SaveResultsFolder'))
    Options.SaveResultsFolder='';
  end
  if not(matchField(Options,'SaveResultsMatVersion'))
    Options.SaveResultsMatVersion='';
  end
  diary off;  diary(sprintf('%s/output/%s%s.txt',fileparts(which(mfilename)),...
                                                               Options.SaveResultsFolder,FileName));
end

% Avoid parpool to start automatically
ParAux=parallel.Settings;
ParAux.Pool.AutoCreate=false;

% Set digits
if matchField(Simulation,'Digits') && Simulation.Digits~=16
  mp.Digits(Simulation.Digits);
else
  Simulation.Digits=16;
  if contains(path,'advanpix')
    rmpath(extractBefore(path,strfind(path,'advanpix')+8));
  end
end
[Parameters(:).MP]=deal(MP);

% Get number of discretizatios
Simulation.NumDiscretizations=numel(Parameters);

% Create formulation object
Formulation=cell(1,Simulation.NumDiscretizations);
for iD=1:Simulation.NumDiscretizations
  Formulation{iD}=eval(Parameters(iD).Formulation);
end

% Get time/frequency domain
Simulation.Domain=Formulation{1}.Domain;

% Get number of simulations
switch Simulation.Type
  case 'SingleSimulation'
    Simulation.NumSimulations=1;
  case 'ConvergenceSpace'
    if exist('MeshFile','var')
      NumMeshesAux=size(MeshFile,2);
      EmptyMeshAux=strcmp(MeshFile,'None');
    elseif exist('Distmesh','var')
      NumMeshesAux=size(Distmesh.ElementSize,2);
      EmptyMeshAux=isnan(Distmesh.ElementSize);
    else
      NumMeshesAux=size(Mesh.MinElementSize,2);
      EmptyMeshAux=isnan(Mesh.MinElementSize);
    end
    Simulation.NumSimulations=numel(Parameters(1).Degree)*NumMeshesAux-sum(EmptyMeshAux(:));
  case 'ConvergenceTime'
    Simulation.NumSimulations=numel(Time.BDFOrder)*size(Time.TimeStepSize,2)-...
                                                                   sum(isnan(Time.TimeStepSize(:)));
  case 'ParametricStudy'
    Simulation.NumSimulations=numel(Parameters(1).(Simulation.ParameterStudy));
  case 'ScalingStrong'
    Simulation.NumSimulations=numel(Simulation.NumProcessors);
  case 'ScalingWeak'
    Simulation.NumSimulations=(numel(Simulation.NumProcessors)==size(MeshFile,2))*...
                               numel(Simulation.NumProcessors);
end

% Initialize results
Results=struct();

% Run each simulation
iS1=1;
iS2=0;
for iS=1:Simulation.NumSimulations
  
  % Set parameters for each simulation
  switch Simulation.Type
    case 'SingleSimulation'
      iS1=1;
      iS2=1;
    case 'ConvergenceSpace'
      iS2=iS2+1;
      if iS2>NumMeshesAux || (size(EmptyMeshAux,1)>1 && EmptyMeshAux(iS1,iS2))
        iS1=iS1+1;
        iS2=1;
      end
      for iD=1:Simulation.NumDiscretizations
        Parameters(iD).Degree=Input.Parameters(iD).Degree(iS1);
        if matchField(Parameters,'MeshRefinement')
          Parameters(iD).MeshRefinement=Input.Parameters(iD).MeshRefinement(iS2);
        end
      end
      for iD=1:Simulation.NumDiscretizations
        BoundariesNamesAux=fieldnames(Input.Boundaries(iD));
        if size(Input.Boundaries(iD).(BoundariesNamesAux{1}),1)>1
          for iBN=1:numel(BoundariesNamesAux)
            if not(isempty(Input.Boundaries(iD).(BoundariesNamesAux{iBN})))
              Boundaries(iD).(BoundariesNamesAux{iBN})=...
                                         Input.Boundaries(iD).(BoundariesNamesAux{iBN})(iS2,:); %#ok
            end
          end
        end
      end
      if exist('Distmesh','var')
        Distmesh.ElementSize=Input.Distmesh.ElementSize(min(end,iS1),min(end,iS2));
      elseif not(exist('MeshFile','var'))
        Mesh.MinElementSize=Input.Mesh.MinElementSize(min(end,iS1),min(end,iS2));
        Mesh.MaxElementSize=Input.Mesh.MaxElementSize(min(end,iS1),min(end,iS2));
        Mesh.MeshGradation=Input.Mesh.MeshGradation;
      end
      if strcmp(Time.TimeDependent,'yes')
        Time.TimeStepSize=Input.Time.TimeStepSize(min(end,iS1),min(end,iS2));
      end
    case 'ConvergenceTime'
      iS2=iS2+1;
      if iS2>size(Input.Time.TimeStepSize,2) || (size(isnan(Input.Time.TimeStepSize),1)>1 && ...
                                                      isnan(Input.Time.TimeStepSize(iS1,iS2)))
        iS1=iS1+1;
        iS2=1;
      end
      Time.BDFOrder=Input.Time.BDFOrder(iS1);
      Time.TimeStepSize=Input.Time.TimeStepSize(min(end,iS1),min(end,iS2));
      if not(exist('MeshFile','var'))
        Mesh.MinElementSize=Input.Mesh.MinElementSize(min(end,iS1),min(end,iS2));
        Mesh.MaxElementSize=Input.Mesh.MaxElementSize(min(end,iS1),min(end,iS2));
        Mesh.MeshGradation=Input.Mesh.MeshGradation;
      end
    case 'ParametricStudy'
      iS2=iS2+1;
      for iD=1:Simulation.NumDiscretizations
        Parameters(iD).(Simulation.ParameterStudy)=...
                                              Input.Parameters(iD).(Simulation.ParameterStudy)(iS2);
      end
    case 'ScalingStrong'
      iS2=iS2+1;
      Simulation.NumProcessors=Input.Simulation.NumProcessors(iS2);
    case 'ScalingWeak'
      iS2=iS2+1;
      Simulation.NumProcessors=Input.Simulation.NumProcessors(iS2);
      for iD=1:Simulation.NumDiscretizations
        BoundariesNamesAux=fieldnames(Input.Boundaries(iD));
        if size(Input.Boundaries(iD).(BoundariesNamesAux{1}),1)>1
          for iBN=1:numel(BoundariesNamesAux)
            if not(isempty(Input.Boundaries(iD).(BoundariesNamesAux{iBN})))
              Boundaries(iD).(BoundariesNamesAux{iBN})=...
                                              Input.Boundaries(iD).(BoundariesNamesAux{iBN})(iS2,:);
            end
          end
        end
      end
  end
  
  % Start parpool if needed
  PoolAux=gcp; if isempty(PoolAux); PoolAux=[]; PoolAux.NumWorkers=1; end
  if matchField(Simulation,'NumProcessors') && Simulation.NumProcessors~=PoolAux.NumWorkers
    if PoolAux.NumWorkers>1
      fprintf('\n'); delete(PoolAux);
    end
    if Simulation.NumProcessors>1
      fprintf('\n'); parpool(parcluster('local'),Simulation.NumProcessors);
    end
  end
  
  fprintf('\n--------------------------------------------------------------------------------')
  
  % Print file name
  fprintf('\nInput file: %s\n',FileName);
  
  % Skip pre-processing
  if matchField(Simulation,'SkipPreProcessing','yes')
    fprintf('\nPre-processing  file: output/%s%s.mat',Options.SaveResultsFolder,...
                                                                  Simulation.PreProcessingFileName);
    load(sprintf('%s/output/%s%s.mat',fileparts(which(mfilename)),Options.SaveResultsFolder,...
                        Simulation.PreProcessingFileName),'-regexp','^(?!(FileName|Simulation)$).');
  else
    
    % Timer
    TimerPreProcessingAux=tic;
    fprintf('\nPre-processing  started...');
    
    % Create model
    for iD=1:Simulation.NumDiscretizations
      Model(iD)=createpde(); %#ok
    end

    % Get geometry and linear mesh
    if exist('MeshFile','var')
      load(['geometry/',MeshFile{min(end,iS1),min(end,iS2)},'.mat']);
      Mesh(1:Simulation.NumDiscretizations)=Mesh(1:Simulation.NumDiscretizations);
      for iD=1:Simulation.NumDiscretizations
        Mesh(iD).DegreeK=Mesh(iD);
        Mesh(iD).Nodes=Mesh(iD).Nodes(:,1:max(max(Mesh(iD).Elements(1:size(Mesh(iD).Nodes,1)+1,:))));
        Mesh(iD).Elements=Mesh(iD).Elements(1:size(Mesh(iD).Nodes,1)+1,:);
        [Geometry(iD),Mesh(iD).Degree1]=geometryFromMesh(Model(iD),double(Mesh(iD).Nodes),...
                                                                          Mesh(iD).Elements); %#ok
        [~,iOrder]=ismember(Mesh(iD).Degree1.Elements',...
                            Mesh(iD).DegreeK.Elements(1:(size(Mesh(iD).Nodes,1)+1),:)','rows');
        Mesh(iD).DegreeK.Elements=Mesh(iD).DegreeK.Elements(:,iOrder);
      end
    elseif exist('Mesh','var') && matchField(Mesh,'Nodes') && matchField(Mesh,'Elements')
      Mesh(1:Simulation.NumDiscretizations)=Mesh(1:Simulation.NumDiscretizations);
      for iD=1:Simulation.NumDiscretizations
        Mesh(iD).DegreeK=Mesh(iD);
        Mesh(iD).Nodes=Mesh(iD).Nodes(:,1:max(max(Mesh(iD).Elements(1:size(Mesh(iD).Nodes,1)+1,:))));
        Mesh(iD).Elements=Mesh(iD).Elements(1:size(Mesh(iD).Nodes,1)+1,:);
        [Geometry(iD),Mesh(iD).Degree1]=geometryFromMesh(Model(iD),double(Mesh(iD).Nodes),...
                                                                          Mesh(iD).Elements);
        [~,iOrder]=ismember(Mesh(iD).Degree1.Elements',...
                            Mesh(iD).DegreeK.Elements(1:(size(Mesh(iD).Nodes,1)+1),:)','rows');
        Mesh(iD).DegreeK.Elements=Mesh(iD).DegreeK.Elements(:,iOrder);
      end
    elseif exist('Distmesh','var')
      Mesh=struct();
      for iD=1:Simulation.NumDiscretizations
        [Mesh(iD).Nodes,Mesh(iD).Elements]=distmesh(Distmesh(iD).DistanceFunction,...
                                       Distmesh(iD).ElementSizeFunction,Distmesh(iD).ElementSize,...
                                       Distmesh(iD).BoundingBox,Distmesh(iD).FixedNodePositions);
        [Geometry(iD),Mesh(iD).Degree1]=geometryFromMesh(Model(iD),Mesh(iD).Nodes',...
                                                                   Mesh(iD).Elements');
      end
    elseif exist('Geometry','var') && isa(Geometry,'numeric')
      Mesh(1:Simulation.NumDiscretizations)=Mesh(1:Simulation.NumDiscretizations);
      for iD=1:Simulation.NumDiscretizations
        Model(iD).Geometry=geometryFromEdges(Model(iD),Geometry(:,:,iD));
        Mesh(iD).Degree1=generateMesh(Model(iD),'GeometricOrder','Linear',...
                                  'Hmin' ,Mesh(iD).MinElementSize,'Hmax',Mesh(iD).MaxElementSize,...
                                  'Hgrad',Mesh(iD).MeshGradation);
      end
      clear Geometry;
      for iD=1:Simulation.NumDiscretizations
        Geometry(iD)=Model(iD).Geometry; %#ok
      end
    elseif  exist('Geometry','var') && (isa(Geometry,'pde.DiscreteGeometry') || ...
                                        isa(Geometry,'pde.AnalyticGeometry'))
      Mesh(1:Simulation.NumDiscretizations)=Mesh(1:Simulation.NumDiscretizations);
      for iD=1:Simulation.NumDiscretizations
        Model(iD).Geometry=Geometry(iD);
        Mesh(iD).Degree1=generateMesh(Model(iD),'GeometricOrder','Linear',...
                                  'Hmin' ,Mesh(iD).MinElementSize,'Hmax',Mesh(iD).MaxElementSize,...
                                  'Hgrad',Mesh(iD).MeshGradation);
      end
    end

    % Space dimension
    Sizes=struct();
    for iD=1:Simulation.NumDiscretizations
      Sizes(iD).NumSpaceDim=1+(Geometry(iD).NumFaces>0)+(Geometry(iD).NumCells>0);
    end
    
    % Assign default uniform nodes distribution
    if not(matchField(Parameters,'NodesDistribution'))
      [Parameters(:).NodesDistribution]=deal('Uniform');
    end

    % Get discretization type
    for iD=1:Simulation.NumDiscretizations
      Parameters(iD).DiscretizationType=Formulation{iD}.DiscretizationType;
    end
    
    % Create reference element
    for iD1=1:Simulation.NumDiscretizations
      for iD2=1:Simulation.NumDiscretizations
        RefElement(iD1,iD2)=createReferenceElement(Sizes(iD1).NumSpaceDim,...
                       Parameters(iD1).Degree,max(Parameters(iD1).Degree,Parameters(iD2).Degree),...
                       Parameters(iD1).NodesDistribution,Simulation.Digits); %#ok
      end
    end
    for iD=1:Simulation.NumDiscretizations
      RefElement(iD,iD).Degree1=createReferenceElement(Sizes(iD).NumSpaceDim,1,1,...
                                               Parameters(iD1).NodesDistribution,Simulation.Digits);
      RefElement(iD,iD).MeshDegree1=createReferenceElement(Sizes(iD).NumSpaceDim,...
                        1,Parameters(iD).Degree,Parameters(iD).NodesDistribution,Simulation.Digits);
      RefElement(iD,iD).MeshDegreeK=createReferenceElement(Sizes(iD).NumSpaceDim,....
                                                Parameters(iD).Degree,Parameters(iD).Degree,...
                                                Parameters(iD).NodesDistribution,Simulation.Digits);
    end
    for iD=1:Simulation.NumDiscretizations
      if strcmp(Parameters(iD).DiscretizationType,'HDG') && ...
         strcmp(Parameters(iD).PostProcessingHDG,'yes')
        RefElement(iD,iD).PostLow=createReferenceElement(Sizes(iD).NumSpaceDim,...
                                                Parameters(iD).Degree  ,Parameters(iD).Degree+1,...
                                                Parameters(iD).NodesDistribution,Simulation.Digits);
        RefElement(iD,iD).Post=createReferenceElement(Sizes(iD).NumSpaceDim,...
                                                Parameters(iD).Degree+1,Parameters(iD).Degree+1,...
                                                Parameters(iD).NodesDistribution,Simulation.Digits);
      end
    end
    
    % Number of components
    for iD=1:Simulation.NumDiscretizations
      Sizes(iD).NumGlobalComp=Formulation{iD}.NumGlobalComp(Sizes(iD).NumSpaceDim);
      if strcmp(Parameters(iD).DiscretizationType,'HDG')
        Sizes(iD).NumVoigtComp=Formulation{iD}.NumVoigtComp(Sizes(iD).NumSpaceDim); 
        Sizes(iD).NumLocalComp=Formulation{iD}.NumLocalComp(Sizes(iD).NumSpaceDim);
        if strcmp(Parameters(iD).PostProcessingHDG,'yes')
          Sizes(iD).NumPostComp=Formulation{iD}.NumPostComp(Sizes(iD).NumSpaceDim);
        end
      end
    end
    
    % Initialize regions
    if not(exist('Regions','var'))
      Regions=[];
    end
    
    % Get faces and BCs
    [Faces,Elements,BCs]=getFacesAndBCs(Simulation,Geometry,Mesh,Boundaries,Regions,...
                                                                       Parameters,RefElement,Sizes);
    
    % Degree-independent sizes
    for iD=1:Simulation.NumDiscretizations
      Sizes(iD).NumElements=size(Mesh(iD).Degree1.Elements,2);
      Sizes(iD).NumElementFaces=size(RefElement(iD,iD).Degree1.FaceNodesElem,1);
    end
    
    % Generate high order mesh
    Mesh(1:Simulation.NumDiscretizations)=Mesh(1:Simulation.NumDiscretizations);
    for iD=1:Simulation.NumDiscretizations
      if Parameters(iD).Degree==1
        Mesh(iD).Nodes=Mesh(iD).Degree1.Nodes;
        Mesh(iD).Elements=Mesh(iD).Degree1.Elements;
      else
        if not(matchField(Mesh,'DegreeK')) || ...
           size(Mesh(iD).DegreeK.Elements,1)==size(Mesh(iD).Degree1.Elements,1)
          Mesh(iD).DegreeK=generateMeshHighOrder(Mesh(iD).Degree1,RefElement(iD,iD).MeshDegree1,...
                                                                  RefElement(iD,iD).MeshDegreeK);
        end
        Mesh(iD).Nodes=Mesh(iD).DegreeK.Nodes;
        Mesh(iD).Elements=Mesh(iD).DegreeK.Elements;
      end
      Mesh(iD).MaxElementSize=Mesh(iD).Degree1.MaxElementSize;
      Mesh(iD).MinElementSize=Mesh(iD).Degree1.MinElementSize;
      Mesh(iD).MeshGradation= Mesh(iD).Degree1.MeshGradation;
    end
    
    % Generate mesh for HDG post-processing
    for iD=1:Simulation.NumDiscretizations
      if strcmp(Parameters(iD).DiscretizationType,'HDG') && ...
         strcmp(Parameters(iD).PostProcessingHDG,'yes')
        MeshPostAux=generateMeshHighOrder(Mesh(iD),RefElement(iD,iD).PostLow,...
                                                   RefElement(iD,iD).Post);
        Mesh(iD).Post.Nodes=MeshPostAux.Nodes;
        Mesh(iD).Post.Elements=MeshPostAux.Elements;
      end
    end
    
    % Sizes
    for iD=1:Simulation.NumDiscretizations
      Sizes(iD).NumNodes=size(Mesh(iD).Nodes,2);
      Sizes(iD).NumElementNodes=size(RefElement(iD,iD).NodesCoordElem,1);
      Sizes(iD).NumFaceNodes=size(RefElement(iD,iD).NodesCoordFace,1);
      if Sizes(iD).NumSpaceDim==2
        Sizes(iD).NumEdgeNodes=Sizes(iD).NumFaceNodes;
      elseif Sizes(iD).NumSpaceDim==3
        Sizes(iD).NumEdgeNodes=size(RefElement(iD,iD).ShapeFunctionsEdge,2);
      end
      if strcmp(Parameters(iD).DiscretizationType,'HDG') && ...
         strcmp(Parameters(iD).PostProcessingHDG,'yes')
        Sizes(iD).NumElementNodesPost=size(RefElement(iD,iD).Post.NodesCoordElem,1);
        Sizes(iD).NumFaceNodesPost=size(RefElement(iD,iD).Post.NodesCoordFace,1);
      end
    end
    
    % Modify mesh for HDG
    for iD=1:Simulation.NumDiscretizations
      if strcmp(Parameters(iD).DiscretizationType,'HDG')
        Mesh(iD).Nodes=Mesh(iD).Nodes(:,Mesh(iD).Elements(1:end));
        Mesh(iD).Elements=reshape(1:Sizes(iD).NumElementNodes*Sizes(iD).NumElements,...
                                   [Sizes(iD).NumElementNodes,Sizes(iD).NumElements]);
        if strcmp(Parameters(iD).PostProcessingHDG,'yes')
          Mesh(iD).Post.Nodes=Mesh(iD).Post.Nodes(:,Mesh(iD).Post.Elements(1:end));
          Mesh(iD).Post.Elements=reshape(1:Sizes(iD).NumElementNodesPost*Sizes(iD).NumElements,...
                                          [Sizes(iD).NumElementNodesPost,Sizes(iD).NumElements]);
        end
      end
    end
    
    % Add z-component for 2D
    for iD=1:Simulation.NumDiscretizations
      if Sizes(iD).NumSpaceDim==2
        Mesh(iD).Nodes=[Mesh(iD).Nodes;zeros(1,size(Mesh(iD).Nodes,2))];
        if strcmp(Parameters(iD).DiscretizationType,'HDG') && ...
           strcmp(Parameters(iD).PostProcessingHDG,'yes')
          Mesh(iD).Post.Nodes=[Mesh(iD).Post.Nodes;zeros(1,size(Mesh(iD).Post.Nodes,2))];
        end
      end
    end
    
    % Store initial mesh
    for iD=1:Simulation.NumDiscretizations
      Mesh(iD).NodesInitial=Mesh(iD).Nodes;
    end
    
    % Store element data for effective use of parfor
    for iD=1:Simulation.NumDiscretizations
      Elements(iD).Nodes=mat2cell(Mesh(iD).Nodes(:,Mesh(iD).Elements(:)),...
        3,ones(Sizes(iD).NumElements,1)*Sizes(iD).NumElementNodes)';
      Elements(iD).Element=mat2cell(Mesh(iD).Elements,...
        Sizes(iD).NumElementNodes,ones(Sizes(iD).NumElements,1))';
      Elements(iD).GlobalLocal=Faces(iD,iD).GlobalLocal;
      Elements(iD).NodesInitial=mat2cell(Mesh(iD).NodesInitial(:,Mesh(iD).Elements(:)),...
        3,ones(Sizes(iD).NumElements,1)*Sizes(iD).NumElementNodes)';
    end
    
    % Number of nodes
    for iD=1:Simulation.NumDiscretizations
      if strcmp(Parameters(iD).DiscretizationType,'CG')
        Sizes(iD).NumGlobalNodes=Sizes(iD).NumNodes;
      elseif strcmp(Parameters(iD).DiscretizationType,'HDG')
        Sizes(iD).NumFaces=max(Faces(iD,iD).Connectivity(:));
        Sizes(iD).NumGlobalNodes=Sizes(iD).NumFaces*Sizes(iD).NumFaceNodes;
        Sizes(iD).NumLocalNodes= Sizes(iD).NumElements*Sizes(iD).NumElementNodes;
      end
    end
    
    % Sizes of elemental contributions
    for iD=1:Simulation.NumDiscretizations
      if strcmp(Parameters(iD).DiscretizationType,'CG')
        Sizes(iD).NumElementRhsCoef=Sizes(iD).NumGlobalComp*Sizes(iD).NumElementNodes;
      elseif strcmp(Parameters(iD).DiscretizationType,'HDG')
        Sizes(iD).NumElementRhsCoef=Sizes(iD).NumGlobalComp*...
                                    Sizes(iD).NumElementFaces*Sizes(iD).NumFaceNodes;
      end
      Sizes(iD).NumElementLhsCoef=Sizes(iD).NumElementRhsCoef^2;
    end
    
    % Interface related sizes
    if matchField(Boundaries,'Interface')
      for iD1=1:Simulation.NumDiscretizations
        for iD2=1:Simulation.NumDiscretizations
          Sizes(iD1).NumFacesInterface(iD2)=size(Faces(iD1,iD2).Interface,1);
          Sizes(iD1).NumElementLhsCoupCoef(iD2)=Sizes(iD1).NumElementRhsCoef*...
                                                Sizes(iD2).NumElementRhsCoef;
        end
      end
    end
    
    % Define matrix pattern
    Block=defineMatrixPattern(Simulation,Parameters,Mesh,Faces,Sizes);
    
    % Get fixed DOFs
    Block=getFixedDofs(Block,Simulation,Mesh,Boundaries,Parameters,RefElement,Sizes);
    System.DOFsFixed=[];
    for iD=1:Simulation.NumDiscretizations
      System.DOFsFixed=[System.DOFsFixed,Block(iD,iD).DOFsFixed+...
                        sum(arrayfun(@(S) S.NumGlobalNodes*S.NumGlobalComp,Sizes(1:iD-1)))];
    end
    
    % Get periodic DOFs
    Block=getPeriodicDofs(Block,Simulation,Mesh,Boundaries,Parameters,RefElement,Sizes);
    System.DOFsMaster=[]; System.DOFsSlave=[];
    for iD=1:Simulation.NumDiscretizations
      System.DOFsMaster=[System.DOFsMaster,Block(iD,iD).DOFsMaster+...
                         sum(arrayfun(@(S) S.NumGlobalNodes*S.NumGlobalComp,Sizes(1:iD-1)))];
      System.DOFsSlave=[System.DOFsSlave,Block(iD,iD).DOFsSlave+...
                        sum(arrayfun(@(S) S.NumGlobalNodes*S.NumGlobalComp,Sizes(1:iD-1)))];
    end
        
    % Default settings
    if strcmp(Time.TimeDependent,'yes')
      for iD=1:Simulation.NumDiscretizations
        Parameters(iD).TimeDerOrder=Formulation{iD}.TimeDerOrder;
      end
      if not(matchField(Parameters,'PredictorDegree'))
        [Parameters(:).PredictorDegree]=deal(0);
      end
    end
    if strcmp(System.Nonlinear,'no')
      System.Tolerance=1e+12;
      System.MaxIterations=1;
    end
    if strcmp(Time.TimeDependent,'yes')
      Time.NumTimeSteps=double(round((Time.FinalTime-Time.InitialTime)/Time.TimeStepSize));
    else
      Time.InitialTime=0;
      Time.FinalTime=0;
      Time.TimeStepSize=0;
      Time.NumTimeSteps=1;
    end
    Time.Time=Time.InitialTime;
    
    % BDF parameters
    if strcmp(Time.TimeDependent,'yes')
      Time.BDF1stDer=['[     1,      -1,      0,       0,     0,       0,        0;',...
                      '    3/2,      -2,    1/2,       0,     0,       0,        0;',...
                      '   11/6,      -3,    3/2,    -1/3,     0,       0,        0;',...
                      '  25/12,      -4,      3,    -4/3,   1/4,       0,        0;',...
                      ' 137/60,      -5,      5,   -10/3,   5/4,    -1/5,        0;',...
                      '  49/20,      -6,   15/2,   -20/3,  15/4,    -6/5,      1/6]'];
      Time.BDF2ndDer=['[     1,      -2,      1,       0,     0,       0,        0,     0;',...
                      '      2,      -5,      4,      -1,     0,       0,        0,     0;',...
                      '  35/12,   -26/3,   19/2,   -14/3, 11/12,       0,        0,     0;',...
                      '   15/4,   -77/6,  107/6,     -13, 61/12,    -5/6,        0,     0;',...
                      ' 203/45,   -87/5,  117/4,  -254/9,  33/2,   -27/5,  137/180,     0;',...
                      ' 469/90, -223/10, 879/20, -949/18,    41, -201/10, 1019/180, -7/10]'];
      if Simulation.Digits==16
        Time.BDF1stDer=eval(Time.BDF1stDer);
        Time.BDF2ndDer=eval(Time.BDF2ndDer);
      else
        Time.BDF1stDer=mp(Time.BDF1stDer);
        Time.BDF2ndDer=mp(Time.BDF2ndDer);
      end
    end
    
    % Colors
    Colors={'b','r','k','g','m','c','y','b','r','k','g','m','c','y'};
    
    % Styles
    Styles={'-','--',':','-.','-','--',':','-.','-','--',':','-.'};
    
    % Markers
    Markers={'o','s','d','p','^','x','*','h','o','s','d','p'};
    
    % Get axis limits
    MinMaxAux=minmax(double([Mesh.Nodes]));
    XMinMaxAux=MinMaxAux(1,:)+[-1,+1]*(MinMaxAux(1,1)==MinMaxAux(1,2));
    YMinMaxAux=MinMaxAux(2,:)+[-1,+1]*(MinMaxAux(2,1)==MinMaxAux(2,2));
    ZMinMaxAux=MinMaxAux(3,:)+[-1,+1]*(MinMaxAux(3,1)==MinMaxAux(3,2));
    
    % Plot geometry
    if strcmp(Options.PlotGeometry,'yes')
      figure('Color','w');
      for iD=1:Simulation.NumDiscretizations
        subplot(1+(numel(Parameters)-1)*(abs(diff(XMinMaxAux))> abs(diff(YMinMaxAux))),...
                1+(numel(Parameters)-1)*(abs(diff(XMinMaxAux))<=abs(diff(YMinMaxAux))),iD)
        if Sizes(iD).NumSpaceDim==2
          Plots=pdegplot(Geometry(iD),'EdgeLabels','on');
          set(Plots(1),'Color',Colors{iD})
        elseif Sizes(iD).NumSpaceDim==3
          Plots=pdegplot(Geometry(iD),'FaceLabels','on','FaceAlpha',0.25);
          set(Plots(1),'FaceColor',Colors{iD})
          ChildrenAux=get(gca,'Children'); delete(ChildrenAux([4,5]));
        end
        grid on; box on; axis equal;
        title(['Geometry',sprintf(' (%d)',iD)*(Simulation.NumDiscretizations>1)]);
        xlim(XMinMaxAux); xlabel('x'); ylim(YMinMaxAux); ylabel('y'); zlim(ZMinMaxAux); zlabel('z');
        pause(eps);
      end
    end
    
    % Plot mesh
    if strcmp(Options.PlotMesh,'yes')
      figure('Color','w');
      for iD=1:Simulation.NumDiscretizations
        subplot(...
          1+(Simulation.NumDiscretizations-1)*(abs(diff(XMinMaxAux))> abs(diff(YMinMaxAux))),...
          1+(Simulation.NumDiscretizations-1)*(abs(diff(XMinMaxAux))<=abs(diff(YMinMaxAux))),iD)
        if Sizes(iD).NumSpaceDim==2
          pdeplot(Mesh(iD).Degree1.Nodes,Mesh(iD).Degree1.Elements,'EdgeColor',Colors{iD});
        elseif Sizes(iD).NumSpaceDim==3
          pdeplot3D(Mesh(iD).Degree1.Nodes,Mesh(iD).Degree1.Elements,...
                                              'FaceColor','w','EdgeColor',Colors{iD},'FaceAlpha',1);
          ChildrenAux=get(gca,'Children'); delete(ChildrenAux([2,3]));
        end
        hold on
        BNAux=fieldnames(BCs(iD));
        Plots=matlab.graphics.chart.primitive.Line.empty();
        for iBN=1:numel(BNAux)
          if Sizes(iD).NumSpaceDim==2
            Plot.(BNAux{iBN})=plot(Mesh(iD).Degree1.Nodes(1,BCs(iD).(BNAux{iBN})),...
                                   Mesh(iD).Degree1.Nodes(2,BCs(iD).(BNAux{iBN})),...
                                   'ok','MarkerFaceColor',Colors{iBN},'DisplayName',BNAux{iBN});
          elseif Sizes(iD).NumSpaceDim==3
            Plot.(BNAux{iBN})=plot3(Mesh(iD).Degree1.Nodes(1,BCs(iD).(BNAux{iBN})),...
                                    Mesh(iD).Degree1.Nodes(2,BCs(iD).(BNAux{iBN})),...
                                    Mesh(iD).Degree1.Nodes(3,BCs(iD).(BNAux{iBN})),...
                                    'ok','MarkerFaceColor',Colors{iBN},'DisplayName',BNAux{iBN});
          end
          if not(isempty(Plot.(BNAux{iBN})))
            Plots(iBN)=Plot.(BNAux{iBN});
          end
        end
        hold off; box on; axis equal;
        title(['Mesh',sprintf(' (%d)',iD)*(Simulation.NumDiscretizations>1)]);
        legend(Plots(arrayfun(@(i) isa(Plots(i),'matlab.graphics.chart.primitive.Line'),...
          1:numel(Plots))),'Location','SouthEast')
        set(gca,'Visible','on')
        xlim(XMinMaxAux); xlabel('x'); ylim(YMinMaxAux); ylabel('y'); zlim(ZMinMaxAux); zlabel('z');
        pause(eps);
      end
    end
    
    % Plot mesh distortion
    if matchField(Options,'PlotMeshDistortion','yes')
      NumColAux=10;
      ColAux=jet(NumColAux);
      figure('Color','w');
      [Mesh(:).Distortion]=deal([]);
      for iD=1:Simulation.NumDiscretizations
        Mesh(iD).Distortion=1-meshQuality(Mesh(iD).Degree1)';
        for iMD=1:NumColAux
          iElemAux=find(Mesh(iD).Distortion>(1/NumColAux*(iMD-1)) & ...
                        Mesh(iD).Distortion<(1/NumColAux*(iMD-0)));
          if not(isempty(iElemAux))
            if Sizes(iD).NumSpaceDim==2
              pdeplot(Mesh(iD).Degree1.Nodes,Mesh(iD).Degree1.Elements(:,iElemAux),...
                        'EdgeColor',ColAux(iMD,:));
            elseif Sizes(iD).NumSpaceDim==3
              pdeplot3D(Mesh(iD).Degree1.Nodes,Mesh(iD).Degree1.Elements(:,iElemAux),...
                        'FaceColor',ColAux(iMD,:),'EdgeColor','k','FaceAlpha',1);
              ChildrenAux=get(gca,'Children'); delete(ChildrenAux([2,3]));
            end
          end
          hold on
        end
      end
      hold off; box on; axis equal; colormap(jet(NumColAux)); caxis([0,1]);
      colorbar('YLim',[0,1],'YTick',linspace(0,1,NumColAux+1))
      title('Mesh distortion');
      set(gca,'Visible','on')
      xlim(XMinMaxAux); xlabel('x'); ylim(YMinMaxAux); ylabel('y'); zlim(ZMinMaxAux); zlabel('z');
      pause(eps);
    end
    
    % Timer
    Timer.PreProcessing=toc(TimerPreProcessingAux);
    Memory.PreProcessing=whos; Memory.PreProcessing=sum(vertcat(Memory.PreProcessing.bytes));
    fprintf('\nPre-processing  completed in %.1f sec (%.1f MB)',Timer.PreProcessing,...
                                                                Memory.PreProcessing/1e6);
    
    % Stop here for pre-processing only
    if matchField(Simulation,'ComputeOnlyPreProcessing','yes')
      fprintf('\n'); break;
    end
    
  end % End of skip pre-processing
  
  %% Processing
  
  % Timer
  TimerProcessingAux=tic;
  fprintf('\nProcessing      started...');
  Timer.Evaluation=[];
  Timer.Solution=[];
  Timer.Local=[];
  
  % Compute quantity of interest at the start of the simulation
  if matchField(Options,'ComputeQuantityStart')
    eval(Options.ComputeQuantityStart);
  end
  
  % Initialization of unknowns
  for iD=1:Simulation.NumDiscretizations
    Block=Formulation{iD}.initializeUnknowns(iD,Block,Parameters,Time,Sizes);
  end

  % Compute initial conditions
  for iD=1:Simulation.NumDiscretizations
    Block=Formulation{iD}.computeInitialConditions(iD,Block,Parameters,Mesh,Faces,Time,...
      RefElement,Sizes);
  end
  
  % Initial and final time step
  if matchField(Simulation,'RestartInitialFinalSteps')
    Time.InitialTimeStep=Simulation.RestartInitialFinalSteps(1);
    Time.FinalTimeStep=Simulation.RestartInitialFinalSteps(2);
  else
    Time.InitialTimeStep=1;
    Time.FinalTimeStep=Time.NumTimeSteps;
  end

  % Send constant variables to the workers to avoid repeated data transfer (only for the time loop)
  RefElement=parallel.pool.Constant(RefElement);
  
  if strcmp(Time.TimeDependent,'yes')
    fprintf('\nTime integration ...............................................................');
  end
  
  % Time loop
  for iT=Time.InitialTimeStep:Time.FinalTimeStep
    Time.TimeStep=double(iT);
    
    % Load restart file
    if matchField(Simulation,'RestartInitialFinalSteps') && ...
       Time.TimeStep==Time.InitialTimeStep && Time.InitialTimeStep>1
      fprintf('\nRestart file:    output/%s%s.mat',Options.SaveResultsFolder,...
                                                                        Simulation.RestartFileName);
      load(sprintf('%s/output/%s%s.mat',fileparts(which(mfilename)),Options.SaveResultsFolder,...
                 Simulation.RestartFileName),'-regexp','^(?!(FileName|Simulation|Time|Options)$).');
      RefElement=parallel.pool.Constant(RefElement);
    end
    
    % Update time
    if strcmp(Time.TimeDependent,'yes')
      Time.Time=Time.InitialTime+Time.TimeStep*Time.TimeStepSize;
    end
    
    % Print time step info
    if strcmp(Time.TimeDependent,'yes')
      fprintf('\nTime = %.2e',Time.Time);
      fprintf(' Step = %.*d/%d',ceil(log10(Time.NumTimeSteps)+eps),Time.TimeStep,Time.NumTimeSteps);
    end
    
    % Select appropriate BDF coefficients
    if strcmp(Time.TimeDependent,'yes')
      if Time.BDFOrder>2
        Time.BDFOrderEff=Time.BDFOrder;
      else
        Time.BDFOrderEff=min(Time.BDFOrder,Time.TimeStep);
      end
      Time.BDF1stDerEff=Time.BDF1stDer(Time.BDFOrderEff,1:Time.BDFOrderEff+1)';
      Time.BDF2ndDerEff=Time.BDF2ndDer(Time.BDFOrderEff,1:Time.BDFOrderEff+2)';
    else
      Time.BDFOrderEff=1;
    end
    
    % Prediction of unknowns
    if strcmp(Time.TimeDependent,'yes')
      for iD=1:Simulation.NumDiscretizations
        PAux=min(Parameters(iD).PredictorDegree,size(Block(iD,iD).SolutionOld,3)-1);
        Caux=round((-1)^PAux*((-PAux:1)'.^(0:PAux+1)'\[zeros(PAux+1,1);factorial(PAux+1)]));
        SolutionPredAux=sum(Block(iD,iD).SolutionOld(:,:,1:PAux+1).*reshape(Caux(2:end),1,1,PAux+1),3);
        if strcmp(Parameters(iD).DiscretizationType,'CG')
          Block(iD,iD).SolutionGlobal=SolutionPredAux;
        elseif strcmp(Parameters(iD).DiscretizationType,'HDG')
          Block(iD,iD).SolutionLocal=SolutionPredAux;
        end
      end
    end
    
    % Initialization of system solutions
    for iD=1:Simulation.NumDiscretizations
      Block(iD,iD).SysSolutionGlobal=MP*zeros(Sizes(iD).NumGlobalNodes*Sizes(iD).NumGlobalComp,1);
      if strcmp(Parameters(iD).DiscretizationType,'HDG')
        Block(iD,iD).SysSolutionLocal=MP*zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumLocalComp);
      end
    end
    
    % Evaluate solution at fixed DOFs
    for iD=1:Simulation.NumDiscretizations
      Block=Formulation{iD}.evaluateSolutionFixedDofs(iD,Block,Parameters,Mesh,Time,Sizes);
    end
    
    % Compute quantity of interest at the current time step
    if matchField(Options,'ComputeQuantityTimeStep')
      eval(Options.ComputeQuantityTimeStep);
    end
    
    % Initialization of residual
    System.ResidualNorm=Inf;
    
    % Initialization of iteration count
    System.Iteration=0;
    
    % Newton iterations
    while System.ResidualNorm > System.Tolerance && ...
          System.Iteration    < System.MaxIterations
      
      % Update iteration count
      System.Iteration=System.Iteration+1;
      
      % Timer
      TimerEvaluationAux=tic;
      
      % Store element data for effective use of parfor
      for iD=1:Simulation.NumDiscretizations
        if strcmp(Parameters(iD).DiscretizationType,'CG')
          Elements(iD).SolutionGlobal=mat2cell(Block(iD,iD).SolutionGlobal(...
            Mesh(iD).Elements(:),:),ones(Sizes(iD).NumElements,1)*...
            Sizes(iD).NumElementNodes,Sizes(iD).NumGlobalComp);
          if strcmp(Time.TimeDependent,'yes')
            Elements(iD).SolutionOld=mat2cell(Block(iD,iD).SolutionOld(...
              Mesh(iD).Elements(:),:,:),ones(Sizes(iD).NumElements,1)*...
              Sizes(iD).NumElementNodes,Sizes(iD).NumGlobalComp,size(Block(iD,iD).SolutionOld,3));
          else
            Elements(iD).SolutionOld=cell(Sizes(iD).NumElements,1);
          end
        elseif strcmp(Parameters(iD).DiscretizationType,'HDG')          
          Elements(iD).SolutionGlobal=mat2cell(Block(iD,iD).SolutionGlobal(...
            cell2mat(Faces(iD,iD).GlobalLocal),:),ones(Sizes(iD).NumElements,1)*...
            Sizes(iD).NumElementFaces*Sizes(iD).NumFaceNodes*Sizes(iD).NumGlobalComp,1);
          Elements(iD).SolutionLocal=mat2cell(Block(iD,iD).SolutionLocal(...
            Mesh(iD).Elements(:),:),ones(Sizes(iD).NumElements,1)*...
            Sizes(iD).NumElementNodes,Sizes(iD).NumLocalComp);
          if strcmp(Time.TimeDependent,'yes')
            Elements(iD).SolutionOld=mat2cell(Block(iD,iD).SolutionOld(...
              Mesh(iD).Elements(:),:,:),ones(Sizes(iD).NumElements,1)*...
              Sizes(iD).NumElementNodes,Sizes(iD).NumLocalComp,size(Block(iD,iD).SolutionOld,3));
          else
            Elements(iD).SolutionOld=cell(Sizes(iD).NumElements,1);
          end
        end
      end
      
      % Compute quantity of interest at the current iteration
      if matchField(Options,'ComputeQuantityIteration')
        eval(Options.ComputeQuantityIteration);
      end
      
      % Build block matrices and residuals
      for iD1=1:Simulation.NumDiscretizations
        for iD2=1:Simulation.NumDiscretizations
          if iD1==iD2
            [Block,Elements]=Formulation{iD1}.buildBlock(iD1,Block,Elements,Simulation,...
              Parameters,Mesh,Faces,Time,RefElement,Sizes);
          else
            Block=Formulation{iD1}.doCoupling(iD1,iD2,Block,Elements,Simulation,...
              Parameters,Mesh,Faces,Time,RefElement,Sizes);
          end
        end
      end
      
      % Timer
      Timer.Evaluation=[Timer.Evaluation,toc(TimerEvaluationAux)];
      
      % Build system
      System.Lhs=struct();
      for iD=1:Simulation.NumDiscretizations
        System.Lhs(iD).LhsGlobal=horzcat(Block(iD,:).LhsGlobal);
      end
      System.Lhs=vertcat(System.Lhs.LhsGlobal); %Block=rmfield(Block,'LhsGlobal');
      System.Rhs=vertcat(Block.RhsGlobal);

      % Strongly impose fixed BCs
      System.Rhs(System.DOFsFixed,1)=0;
      System.Lhs(System.DOFsFixed,:)=0;
      System.Lhs(System.DOFsFixed,System.DOFsFixed)=eye(numel(System.DOFsFixed));

      % Strongly impose periodic BCs
      System.Rhs(System.DOFsMaster,:)=System.Rhs(System.DOFsMaster,:)+System.Rhs(System.DOFsSlave,:);
      System.Lhs(System.DOFsMaster,:)=System.Lhs(System.DOFsMaster,:)+System.Lhs(System.DOFsSlave,:);
      System.Rhs(System.DOFsSlave,1)=0;
      System.Lhs(System.DOFsSlave,:)=0;
      System.Lhs(System.DOFsSlave,System.DOFsSlave) =+eye(numel(System.DOFsSlave));
      System.Lhs(System.DOFsSlave,System.DOFsMaster)=-eye(numel(System.DOFsSlave));

      % Symmetrize matrix
      if matchField(System,'SymmetrizeMatrix','yes') && ...
        (not(matchField(System,'SymmetrizeMatrixOnlyOnce')) || ...
        (matchField(System,'SymmetrizeMatrixOnlyOnce','yes') && ...
         Time.TimeStep<=Time.BDFOrderEff && System.Iteration==1))
          System.Lhs=(System.Lhs+System.Lhs.')/2;
      end
      
      % GLOBAL PROBLEM -----------------------------------------------------------------------------
      % Equilibrate
      if matchField(Solver,'Equilibrate','yes')
        [System.EquilP,System.EquilR,System.EquilC]=equilibrate(System.Lhs);
        System.Lhs=System.EquilR*System.EquilP*System.Lhs*System.EquilC;
        System.Rhs=System.EquilR*System.EquilP*System.Rhs;
      end
      
      % Preconditioning
      if matchField(Solver,'Preconditioner')
        if (Time.TimeStep==1    || strcmp(Solver.PreconditionerAllTimeSteps,'yes')) && ...
           (System.Iteration==1 || strcmp(Solver.PreconditionerAllIterations,'yes'))
          switch Solver.Preconditioner
            case 'none'
              System.PrecondL=[]; System.PrecondU=[];
            case 'ichol'
              System.PrecondL=ichol(System.Lhs);
            case 'ilu'
              [System.PrecondL,System.PrecondU]=ilu(System.Lhs);
          end
        end
        fprintf('\n\n');
      end
      
      % Distribute
      if matchField(System,'Distribute','yes')
        System.Lhs=distributed(System.Lhs);
        System.Rhs=distributed(System.Rhs);
      end
      
      % Timer
      TimerSolutionAux=tic;
      
      % Decompose
      if matchField(System,'Decompose','yes')
        if Time.TimeStep<Time.BDFOrder
          System.LhsDecomposed=System.Lhs;
        elseif Time.TimeStep==Time.BDFOrder
          System.LhsDecomposed=decomposition(System.Lhs);
        end
        System.Lhs=System.LhsDecomposed;
      end
      
      % Store old solution for relaxation
      if matchField(System,'RelaxationParameter')
        if System.Iteration>1
          System.SolutionOld=System.Solution;
        else
          System.SolutionOld=sparse(size(System.Rhs,1),1);
        end
      end
      
      % Solve
      switch Solver.Type
        case 'backslash'
          System.Solution=System.Lhs\System.Rhs;
        case 'pcg'
          System.Solution=pcg(System.Lhs,System.Rhs,...
              Solver.Tolerance,Solver.MaxIterations,System.PrecondL,System.PrecondL');
        case 'minres'
          System.Solution=minres(System.Lhs,System.Rhs,...
              Solver.Tolerance,Solver.MaxIterations,System.PrecondL,System.PrecondL');
        case 'gmres'
           System.Solution=gmres(System.Lhs,System.Rhs,...
             Solver.Restart,Solver.Tolerance,min(Solver.MaxIterations,size(System.Rhs,1)),...
             System.PrecondL,System.PrecondU);
      end
      %System=rmfield(System,'Lhs');
      
      % Timer
      Timer.Solution=[Timer.Solution,toc(TimerSolutionAux)];
      
      % Relax solution
      if matchField(System,'RelaxationParameter')
        System.Solution=System.RelaxationParameter*System.Solution...
                   +(1-System.RelaxationParameter)*System.SolutionOld;
      end
      
      % Gather
      if matchField(System,'Distribute','yes')
        System.Solution=gather(System.Solution);
      end
      
      % Un-equilibrate
      if matchField(Solver,'Equilibrate','yes')
        System.Solution=System.EquilC*System.Solution;
      end
      
      % Extract
      Block(1,1).SysSolutionGlobal=System.Solution((1:size(Block(1,1).RhsGlobal,1)),1);
      for iD=2:Simulation.NumDiscretizations
        Block(iD,iD).SysSolutionGlobal=System.Solution(....
          sum(arrayfun(@(S) S.NumGlobalNodes*S.NumGlobalComp,Sizes(1:iD-1)))...
          +(1:size(Block(iD,iD).RhsGlobal,1)),1);
      end
      %System=rmfield(System,'Solution');
      %Block=rmfield(Block,'RhsGlobal');
      
      % Reshape
      for iD=1:Simulation.NumDiscretizations
        if strcmp(Parameters(iD).DiscretizationType,'CG')
          Block(iD,iD).SysSolutionGlobal=reshape(Block(iD,iD).SysSolutionGlobal,...
                                                [Sizes(iD).NumGlobalNodes,Sizes(iD).NumGlobalComp]);
        end
      end
      
      % Update
      for iD=1:Simulation.NumDiscretizations
        if matchField(System,'IncrementalForm','no')
          Block(iD,iD).SolutionGlobal=Block(iD,iD).SysSolutionGlobal;
        else
          Block(iD,iD).SolutionGlobal=Block(iD,iD).SolutionGlobal+Block(iD,iD).SysSolutionGlobal;
        end
      end
      % --------------------------------------------------------------------------------------------
      
      % LOCAL PROBLEMS -----------------------------------------------------------------------------
      % Timer
      TimerLocalAux=tic;
      
      for iD=1:Simulation.NumDiscretizations
        if strcmp(Parameters(iD).DiscretizationType,'HDG')
          
          % Store element data for effective use of parfor
          Elements(iD).SysSolutionGlobal=mat2cell(Block(iD,iD).SysSolutionGlobal(...
            cell2mat(Faces(iD,iD).GlobalLocal),:),ones(Sizes(iD).NumElements,1)*...
            Sizes(iD).NumElementFaces*Sizes(iD).NumFaceNodes*Sizes(iD).NumGlobalComp,1);
          
          % Solve
          for iElem=1:Sizes(iD).NumElements
            Block(iD,iD).SysSolutionLocal((iElem-1)*Sizes(iD).NumElementNodes+...
                                          (1:Sizes(iD).NumElementNodes),:)=...
                              reshape(Elements(iD).VecLocal{iElem}-...
                              Elements(iD).MatLocal{iElem}*Elements(iD).SysSolutionGlobal{iElem},...
                              Sizes(iD).NumElementNodes,Sizes(iD).NumLocalComp);
          end
          
          % Update
          if matchField(System,'IncrementalForm','no')
            Block(iD,iD).SolutionLocal=Block(iD,iD).SysSolutionLocal;
          else
            Block(iD,iD).SolutionLocal=Block(iD,iD).SolutionLocal+Block(iD,iD).SysSolutionLocal;
          end
        end
      end
      
      % Timer
      Timer.Local=[Timer.Local,toc(TimerLocalAux)];
      % --------------------------------------------------------------------------------------------
      
      % Compute norm of the residual
      System.ResidualNorm=norm(System.Rhs);
      
      % Print Newton iteration info
      if strcmp(System.Nonlinear,'yes')
        fprintf('\n%2sIter = %d','',System.Iteration);
        fprintf('%*s|Res| = %.1e',7-numel(num2str(System.Iteration)),'',System.ResidualNorm);
        fprintf('%6s[te,ts%s] = [%.1f,%.1f%s]','',...
          repmat(',tl',1,any(strcmp({Parameters.DiscretizationType},'HDG'))),...
          Timer.Evaluation(end),Timer.Solution(end),repmat(sprintf(',%.1f',Timer.Local(end)),1,...
          any(strcmp({Parameters.DiscretizationType},'HDG'))));
      end
      
      % Compute quantity of interest at the end of the current iteration
      if matchField(Options,'ComputeQuantityIterationEnd')
        eval(Options.ComputeQuantityIterationEnd);
      end
      
    end % End of the Newton iterations
    
    % Shift old solutions
    if strcmp(Time.TimeDependent,'yes') && Time.TimeStep<Time.NumTimeSteps
      for iD=1:Simulation.NumDiscretizations
        Block(iD,iD).SolutionOld(:,:,2:min(Time.BDFOrderEff+1,Time.BDFOrder)...
                                                               +(Parameters(iD).TimeDerOrder-1))=...
        Block(iD,iD).SolutionOld(:,:,1:min(Time.BDFOrderEff+1,Time.BDFOrder)...
                                                               -(2-Parameters(iD).TimeDerOrder));
        if strcmp(Parameters(iD).DiscretizationType,'CG')
          Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionGlobal;
        elseif strcmp(Parameters(iD).DiscretizationType,'HDG')
          Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionLocal;
        end
      end
    end
    
    % Store results of current time step
    if matchField(Options,'StoreTimeSteps')
      iST=find(Time.TimeStep==eval(Options.StoreTimeSteps));
      if not(isempty(iST))
        for iD=1:Simulation.NumDiscretizations
          Results=Formulation{iD}.storeResults(iD,iST,Results,Block,Simulation,Parameters,Mesh,...
            Time,RefElement,Sizes);
        end
      end
    end
    
    % Compute quantity of interest
    if matchField(Options,'ComputeQuantity')
      eval(Options.ComputeQuantity);
    end
    
    % Plot solution of current time step
    if matchField(Options,'PlotSolution') && ...
       matchField(Options,'PlotSolutionTimeSteps') && ...
       any(Time.TimeStep==eval(Options.PlotSolutionTimeSteps))
      iPlotT=find(Time.TimeStep==eval(Options.StoreTimeSteps));
      plotSolution(Simulation,Mesh,Faces,Results,Parameters,RefElement,Sizes,Options,iPlotT);
    end
    
    % Export to Paraview of current time step
    if matchField(Options,'Export2Paraview') && ...
       matchField(Options,'StoreTimeSteps') && ...
       matchField(Options,'Export2ParaviewTimeSteps') && ...
       any(Time.TimeStep==eval(Options.Export2ParaviewTimeSteps))
      for iD=1:Simulation.NumDiscretizations
        if Simulation.NumDiscretizations==1
          FileNameAux=sprintf('%s_step_%.*d',FileName,...
            ceil(log10(Time.NumTimeSteps)+eps),Time.TimeStep);
        else
          FileNameAux=sprintf('%s_discr_%d_step_%.*d',FileName,iD,...
            ceil(log10(Time.NumTimeSteps)+eps),Time.TimeStep);
        end
        export2Paraview(Mesh(iD),Results(iD),Parameters(iD),Time,Sizes(iD),Options,FileNameAux);
      end
    end
  
    % Save results of current time step
    if matchField(Options,'SaveResults','yes') && ...
       matchField(Options,'StoreTimeSteps') && ...
       matchField(Options,'SaveResultsTimeSteps') && ...
       any(Time.TimeStep==eval(Options.SaveResultsTimeSteps))
      save(sprintf('%s/output/%s%s_step_%.*d',fileparts(which(mfilename)),...
            Options.SaveResultsFolder,FileName,ceil(log10(Time.NumTimeSteps)+eps),Time.TimeStep),...
           'Results',Options.SaveResultsMatVersion);
    end
    
  end % End of the time loop
  
  if strcmp(Time.TimeDependent,'yes')
    fprintf('\n................................................................................');
  end

  % Retrieve the original variables from the workers
  RefElement=RefElement.Value;
  
  % Store element data for effective use of parfor
  for iD=1:Simulation.NumDiscretizations
    if strcmp(Parameters(iD).DiscretizationType,'CG')
      Elements(iD).SolutionGlobal=mat2cell(Block(iD,iD).SolutionGlobal(...
        Mesh(iD).Elements(:),:),ones(Sizes(iD).NumElements,1)*...
        Sizes(iD).NumElementNodes,Sizes(iD).NumGlobalComp);
      if strcmp(Time.TimeDependent,'yes')
        Elements(iD).SolutionOld=mat2cell(Block(iD,iD).SolutionOld(...
          Mesh(iD).Elements(:),:,:),ones(Sizes(iD).NumElements,1)*...
          Sizes(iD).NumElementNodes,Sizes(iD).NumGlobalComp,size(Block(iD,iD).SolutionOld,3));
      else
        Elements(iD).SolutionOld=cell(Sizes(iD).NumElements,1);
      end
    elseif strcmp(Parameters(iD).DiscretizationType,'HDG')
      Elements(iD).SolutionGlobal=mat2cell(Block(iD,iD).SolutionGlobal(...
        cell2mat(Faces(iD,iD).GlobalLocal),:),ones(Sizes(iD).NumElements,1)*...
        Sizes(iD).NumElementFaces*Sizes(iD).NumFaceNodes*Sizes(iD).NumGlobalComp,1);
      Elements(iD).SolutionLocal=mat2cell(Block(iD,iD).SolutionLocal(...
        Mesh(iD).Elements(:),:),ones(Sizes(iD).NumElements,1)*...
        Sizes(iD).NumElementNodes,Sizes(iD).NumLocalComp);
      if strcmp(Time.TimeDependent,'yes')
        Elements(iD).SolutionOld=mat2cell(Block(iD,iD).SolutionOld(...
          Mesh(iD).Elements(:),:,:),ones(Sizes(iD).NumElements,1)*...
          Sizes(iD).NumElementNodes,Sizes(iD).NumLocalComp,size(Block(iD,iD).SolutionOld,3));
      else
        Elements(iD).SolutionOld=cell(Sizes(iD).NumElements,1);
      end
    end
  end
  
  % HDG POST-PROCESSING ----------------------------------------------------------------------------
  for iD=1:Simulation.NumDiscretizations
    if strcmp(Parameters(iD).DiscretizationType,'HDG') && ...
       strcmp(Parameters(iD).PostProcessingHDG,'yes')
      % Build lhs and rhs
      Elements(iD).LhsPost=[]; Elements(iD).RhsPost=[];
      Elements(iD)=Formulation{iD}.doPostProcess(Elements(iD),Simulation,Parameters(iD),...
                                                     Faces(iD,iD),Time,RefElement(iD,iD),Sizes(iD));
      
      % Solve
      SolutionPost=MP*zeros(Sizes(iD).NumElementNodesPost*Sizes(iD).NumPostComp,...
                                                          Sizes(iD).NumElements);
      LhsPostAux=Elements(iD).LhsPost; RhsPostAux=Elements(iD).RhsPost;      
      parfor iElem=1:Sizes(iD).NumElements
        SolutionPostAux=reshape(LhsPostAux{iElem}\RhsPostAux{iElem},...
                                Sizes(iD).NumElementNodesPost,Sizes(iD).NumPostComp)'; %#ok
        SolutionPost(:,iElem)=SolutionPostAux(:);
      end
      Block(iD,iD).SolutionPost=SolutionPost;
      clear SolutionPost;
      
      % Reshape
      Block(iD,iD).SolutionPost=reshape(Block(iD,iD).SolutionPost,Sizes(iD).NumPostComp,[])';
    end
  end
  % ------------------------------------------------------------------------------------------------
  
  % Store final results
  if matchField(Options,'StoreTimeSteps')
    iST=find(Time.TimeStep==eval(Options.StoreTimeSteps));
  else
    iST=1;
  end
  for iD=1:Simulation.NumDiscretizations
    Results=Formulation{iD}.storeResults(iD,iST,Results,Block,Simulation,Parameters,Mesh,Time,...
      RefElement,Sizes);
  end
  
  % Compute Fast Fourier Transform
  if matchField(Options,'ComputeFFT')
    for iD=1:Simulation.NumDiscretizations
      for iF=1:size(Options.ComputeFFT,1)
        if not(isempty(Results(iD).([Options.ComputeFFT{iF,1}])))
          [Results(iD).Frequency,Results(iD).([Options.ComputeFFT{iF,1},'FFT'])]=...
                                computeFFT(Results(iD).Time,Results(iD).(Options.ComputeFFT{iF,1}));
        end
      end
    end
  end
  
  % Timer
  Timer.Processing=toc(TimerProcessingAux);
  Memory.Processing=whos; Memory.Processing=sum(vertcat(Memory.Processing.bytes));
  fprintf('\nProcessing      completed in %.1f sec (%.1f MB)',Timer.Processing,...
                                                              Memory.Processing/1e6);
  
  %% Post-processing
  
  % Timer
  TimerPostProcessingAux=tic;
  fprintf('\nPost-processing started...');
  
  % Compute quantity of interest for error computation
  if matchField(Options,'ComputeQuantityError')
    eval(Options.ComputeQuantityError);
  end
  
  % Plot final solution
  if matchField(Options,'PlotSolution')
    if matchField(Results,'Time')
      iPlotT=numel(Results(1).Time);
    else
      iPlotT=1;
    end
    plotSolution(Simulation,Mesh,Faces,Results,Parameters,RefElement,Sizes,Options,iPlotT);
  end
  
  % Compute error
  if matchField(Options,'ComputeError')
    if size(Options.ComputeError,2)==1
      Options.ComputeError(:,2)={'L2'};
    end
    if size(Options.ComputeError,2)==2
      Options.ComputeError(:,3)={'@(x,y,z) true'};
    end
    for iD=1:Simulation.NumDiscretizations
      for iE=1:size(Options.ComputeError,1)
        if not(isempty(Results(iD).([Options.ComputeError{iE,1}])))
          if not(strcmp(Options.ComputeError{iE,1}(end-3:end),'Post'))
            Results(iD).(sprintf('%sError%s',Options.ComputeError{iE,1:2}))(iS1,iS2)=...
              sqrt(sum(computeError(Time,Mesh(iD),RefElement(iD,iD),Sizes(iD),...
              Results(iD).(Options.ComputeError{iE,1})(:,:,end),...
              Parameters(iD).(Options.ComputeError{iE,1}),...
              Options.ComputeError{iE,2},Options.ComputeError{iE,3},Simulation.Domain).^2));
          else
            Results(iD).(sprintf('%sError%s',Options.ComputeError{iE,1:2}))(iS1,iS2)=...
              sqrt(sum(computeError(Time,Mesh(iD).Post,RefElement(iD,iD).Post,Sizes(iD),...
              Results(iD).(Options.ComputeError{iE,1})(:,:,end),...
              Parameters(iD).(Options.ComputeError{iE,1}(1:end-4)),...
              Options.ComputeError{iE,2},Options.ComputeError{iE,3},Simulation.Domain).^2));
          end
        end
      end
    end
  end
  
  % Compute quantity of interest for after error computation
  if matchField(Options,'ComputeQuantityPostError')
    eval(Options.ComputeQuantityPostError);
  end
  
  % Compute convergence rate for space convergence
  if strcmp(Simulation.Type,'ConvergenceSpace')
    for iD=1:Simulation.NumDiscretizations
      Results(iD).Degree(iS1,1)=Parameters.Degree;
      if isa(Options.CharacteristicElementSize,'char')
        Results(iD).ElementSize(iS1,iS2)=...
                                       Mesh(iD).([Options.CharacteristicElementSize,'ElementSize']);
      elseif isa(Options.CharacteristicElementSize,'numeric')
        Results(iD).ElementSize(iS1,iS2)=Options.CharacteristicElementSize(iS1,iS2);
      end
      if strcmp(Time.TimeDependent,'yes')
        Results(iD).TimeStepSize(iS1,1)=Time.TimeStepSize;
      end
      for iE=1:size(Options.ComputeError,1)
        if not(isempty(Results(iD).([Options.ComputeError{iE,1}])))
          if iS2==1
            Results(iD).(sprintf('%sError%sRate',Options.ComputeError{iE,1:2}))(iS1,iS2)=NaN;
          else
            Results(iD).(sprintf('%sError%sRate',Options.ComputeError{iE,1:2}))(iS1,iS2)=...
              log(Results(iD).(sprintf('%sError%s',Options.ComputeError{iE,1:2}))(iS1,iS2)/...
                  Results(iD).(sprintf('%sError%s',Options.ComputeError{iE,1:2}))(iS1,iS2-1))/...
              log(Results(iD).ElementSize(iS1,iS2)/...
                  Results(iD).ElementSize(iS1,iS2-1));
          end
        end
      end
    end
    
    % Plot space convergence
    if strcmp(Options.PlotConvergence,'yes')
      if iS==1
        figure('Color','w')
        NumPlotsAux=size(Options.ComputeError,1);
        if NumPlotsAux>1
          SizeAux=get(0,'ScreenSize');
          set(gcf,'Position',[(3-ceil(NumPlotsAux/4))/3*SizeAux(3),0,...
                                                       ceil(NumPlotsAux/4)/3*SizeAux(3),SizeAux(4)])
        end
      elseif iS2>1
        for iE=1:NumPlotsAux
          for iD=1:Simulation.NumDiscretizations
            subplot(ceil(NumPlotsAux/ceil(NumPlotsAux/4)),ceil(NumPlotsAux/4),iE)
            if not(isempty(Results(iD).([Options.ComputeError{iE,1}])))
              plot(log2(Results(iD).ElementSize(iS1,iS2-1:iS2)),...
                  log10(Results(iD).(sprintf('%sError%s',...
                        Options.ComputeError{iE,1:2}))(iS1,iS2-1:iS2)),...
                  sprintf('%s',Colors{iS1},Styles{iD},Markers{iS1}));
              text(log2(Results(iD).ElementSize(iS1,iS2)),...
                  log10(Results(iD).(sprintf('%sError%s',...
                        Options.ComputeError{iE,1:2}))(iS1,iS2)),...
                  sprintf('%.1f',Results(iD).(sprintf('%sError%sRate',...
                        Options.ComputeError{iE,1:2}))(iS1,iS2)),...
                 'HorizontalAlignment','Center','VerticalAlignment','Top',...
                 'Color',sprintf('%s',Colors{iS1}))
            end
            hold on
          end
          grid on; box on;
          title(Options.ComputeError{iE,1});
          xlabel('log_{2}(Element size)');
          ylabel(['log_{10}(Error-',Options.ComputeError{iE,2},')']);
          pause(eps);
        end
      end
    end
  end
  
  % Compute convergence rate for time convergence
  if strcmp(Simulation.Type,'ConvergenceTime')
    for iD=1:Simulation.NumDiscretizations
      Results(iD).BDFOrder(iS1,1)=Time.BDFOrder;
      Results(iD).TimeStepSize(iS1,iS2)=Time.TimeStepSize;
      if isa(Options.CharacteristicElementSize,'char')
        Results(iD).ElementSize(iS1,1)=Mesh(iD).([Options.CharacteristicElementSize,'ElementSize']);
      elseif isa(Options.CharacteristicElementSize,'numeric')
        Results(iD).ElementSize(iS1,1)=Options.CharacteristicElementSize(iS1);
      end
      for iE=1:size(Options.ComputeError,1)
        if not(isempty(Results(iD).([Options.ComputeError{iE,1}])))
          if iS2==1
            Results(iD).(sprintf('%sError%sRate',Options.ComputeError{iE,1:2}))(iS1,iS2)=NaN;
          else
            Results(iD).(sprintf('%sError%sRate',Options.ComputeError{iE,1:2}))(iS1,iS2)=...
              log(Results(iD).(sprintf('%sError%s',Options.ComputeError{iE,1:2}))(iS1,iS2)/...
                  Results(iD).(sprintf('%sError%s',Options.ComputeError{iE,1:2}))(iS1,iS2-1))/...
              log(Results(iD).TimeStepSize(iS1,iS2)/...
                  Results(iD).TimeStepSize(iS1,iS2-1));
          end
        end
      end
    end
    
    % Plot time convergence
    if strcmp(Options.PlotConvergence,'yes')
      if iS==1
        figure('Color','w')
        NumPlotsAux=size(Options.ComputeError,1);
        if size(Options.ComputeError,1)>1
          SizeAux=get(0,'ScreenSize');
          set(gcf,'Position',[(3-ceil(NumPlotsAux/4))/3*SizeAux(3),0,...
                                                       ceil(NumPlotsAux/4)/3*SizeAux(3),SizeAux(4)])
        end
      elseif iS2>1
        for iE=1:size(Options.ComputeError,1)
          for iD=1:Simulation.NumDiscretizations
            subplot(ceil(NumPlotsAux/ceil(NumPlotsAux/4)),ceil(NumPlotsAux/4),iE)
            if not(isempty(Results(iD).([Options.ComputeError{iE,1}])))
              plot(log2(Results(iD).TimeStepSize(iS1,iS2-1:iS2)),...
                  log10(Results(iD).(sprintf('%sError%s',...
                        Options.ComputeError{iE,1:2}))(iS1,iS2-1:iS2)),...
                  sprintf('%s',Colors{iS1},Styles{iD},Markers{iS1}));
              text(log2(Results(iD).TimeStepSize(iS1,iS2)),...
                  log10(Results(iD).(sprintf('%sError%s',...
                        Options.ComputeError{iE,1:2}))(iS1,iS2)),...
                  sprintf('%.1f',Results(iD).(sprintf('%sError%sRate',...
                        Options.ComputeError{iE,1:2}))(iS1,iS2)),...
                 'HorizontalAlignment','Center','VerticalAlignment','Top',...
                 'Color',sprintf('%s',Colors{iS1}))
            end
            hold on
          end
          grid on; box on;
          title(Options.ComputeError{iE,1});
          xlabel('log_{2}(Time step size)');
          ylabel(['log_{10}(Error-',Options.ComputeError{iE,2},')']);
          pause(eps);
        end
      end
    end
  end
  
  % Save parameter of parametric study
  if strcmp(Simulation.Type,'ParametricStudy')
    for iD=1:Simulation.NumDiscretizations
      Results(iD).(Simulation.ParameterStudy)(iS2)=Parameters(iD).(Simulation.ParameterStudy);
    end
  end
  
  % Compute speedup for strong scaling
  if strcmp(Simulation.Type,'ScalingStrong')
    for iD=1:Simulation.NumDiscretizations
      Results(iD).NumProcessors(iS2)=Simulation.NumProcessors;
      Results(iD).Timer(iS2)=sum(Timer.(Options.CharacteristicTimer));
      Results(iD).Speedup(iS2)=Results(iD).Timer(1)/Results(iD).Timer(iS2);
    end
    
    % Plot speedup
    if strcmp(Options.PlotSpeedup,'yes')
      if iS==1
        figure('Color','w')
      elseif iS2>1
        plot(Results(1).NumProcessors(iS2-1:iS2),Results(1).NumProcessors(iS2-1:iS2),'k--')
        hold on
        plot(Results(1).NumProcessors(iS2-1:iS2),Results(1).Speedup(1,iS2-1:iS2),'b-o');
        text(Results(1).NumProcessors(iS2),Results(1).Speedup(iS2),...
             sprintf('%.1f',Results(1).Speedup(iS2)),...
             'HorizontalAlignment','Center','VerticalAlignment','Top','Color','b')
        grid on; box on;
        title('Strong scaling');
        xlabel('Number of processors');
        ylabel(['Speedup (Timer.',Options.CharacteristicTimer,')']);
        pause(eps);
      end
    end
  end
  
  % Compute speedup for weak scaling
  if strcmp(Simulation.Type,'ScalingWeak')
    for iD=1:Simulation.NumDiscretizations
      Results(iD).NumProcessors(iS2)=Simulation.NumProcessors;
      Results(iD).NumElements(iS2)=sum([Sizes.NumElements]);
      Results(iD).Timer(iS2)=sum(Timer.(Options.CharacteristicTimer));
      Results(iD).Efficiency(iS2)=Results(iD).Timer(1)/Results(iD).Timer(iS2);
    end
    
    % Plot efficiency
    if strcmp(Options.PlotEfficiency,'yes')
      if iS==1
        figure('Color','w')
      elseif iS2>1
        plot(Results(1).NumProcessors(iS2-1:iS2),[1,1],'k--')
        hold on
        plot(Results(1).NumProcessors(iS2-1:iS2),Results(1).Efficiency(1,iS2-1:iS2),'b-o');
        text(Results(1).NumProcessors(iS2),Results(1).Efficiency(iS2),...
             sprintf('%.1f',Results(1).Efficiency(iS2)),...
             'HorizontalAlignment','Center','VerticalAlignment','Top','Color','b')
        grid on; box on;
        title('Weak scaling');
        xlabel('Number of processors');
        ylabel(['Efficiency (Timer.',Options.CharacteristicTimer,')']);
        ylim([0,1.2])
        pause(eps);
      end
    end
  end
  
  % Compute quantity of interest at the end of the simulation
  if matchField(Options,'ComputeQuantityEnd')
    eval(Options.ComputeQuantityEnd);
  end
  
  % Export to Paraview
  if matchField(Options,'Export2Paraview')
    for iD=1:Simulation.NumDiscretizations
      if Simulation.NumDiscretizations==1
        FileNameAux=FileName;
      else
        FileNameAux=[FileName,'_discr_',num2str(iD)];
      end
      export2Paraview(Mesh(iD),Results(iD),Parameters(iD),Time,Sizes(iD),Options,FileNameAux);
      if strcmp(Parameters(iD).DiscretizationType,'HDG') && ...
         strcmp(Parameters(iD).PostProcessingHDG,'yes')
        FileNamePostAux=[FileNameAux,'_post'];
        export2Paraview(Mesh(iD).Post,Results(iD),Parameters(iD),Time,Sizes(iD),Options,...
          FileNamePostAux);
      end
    end
  end
  
  % Timer
  Timer.PostProcessing=toc(TimerPostProcessingAux);
  Memory.PostProcessing=whos; Memory.PostProcessing=sum(vertcat(Memory.PostProcessing.bytes));
  fprintf('\nPost-processing completed in %.1f sec (%.1f MB)',Timer.PostProcessing,...
                                                              Memory.PostProcessing/1e6);
  
  % Print errors and convergence rates
  if matchField(Options,'ComputeError')
    fprintf('\n\nVariable%21sError%6sNorm','','');
    if strcmp(Simulation.Type,'ConvergenceSpace') || strcmp(Simulation.Type,'ConvergenceTime')
      fprintf('%3sRate','');
    end
    for iD=1:Simulation.NumDiscretizations
      for iE=1:size(Options.ComputeError,1)
        if not(isempty(Results(iD).([Options.ComputeError{iE,1}])))
          fprintf('\n%s',Options.ComputeError{iE,1});
          fprintf('%s',(Simulation.NumDiscretizations>1)*sprintf('_%d',iD));
          fprintf('%*s',27-numel(Options.ComputeError{iE,1}),'');
          fprintf('%.2e',Results(iD).(sprintf('%sError%s',Options.ComputeError{iE,1:2}))(iS1,iS2));
          fprintf('%3s%s','',Options.ComputeError{iE,2});
          if strcmp(Simulation.Type,'ConvergenceSpace') || strcmp(Simulation.Type,'ConvergenceTime')
            fprintf('%*s%.2f',7-numel(Options.ComputeError{iE,2}),'',...
                      Results(iD).(sprintf('%sError%sRate',Options.ComputeError{iE,1:2}))(iS1,iS2));
          end
        end
      end
    end
  end
  
  % Print simulation info
  if strcmp(Simulation.Type,'ConvergenceSpace')
    fprintf('\n\nSimulation with %s completed',strjoin(arrayfun(@(iD)...
           sprintf('k = %d and h = %.2e',Parameters(iD).Degree,Results(iD).ElementSize(iS1,iS2)),...
           1:Simulation.NumDiscretizations,'UniformOutput',false),'\n       and with '));
  elseif strcmp(Simulation.Type,'ConvergenceTime')
    fprintf('\n\nSimulation with BDFo = %d and dt = %.2e completed',...
                                                    Time.BDFOrder,Results(1).TimeStepSize(iS1,iS2));
  elseif strcmp(Simulation.Type,'ParametricStudy')
    fprintf('\n\nSimulation with %s completed',strjoin(arrayfun(@(iD)...
                   sprintf('%s = %.2e',Simulation.ParameterStudy,...
                   Results(iD).(Simulation.ParameterStudy)(iS2)),1:Simulation.NumDiscretizations,...
                   'UniformOutput',false),'\n       and with '));
  elseif strcmp(Simulation.Type,'ScalingStrong')
    fprintf('\n\nSimulation with np = %d completed',Results(1).NumProcessors(iS2));
  elseif strcmp(Simulation.Type,'ScalingWeak')
    fprintf('\n\nSimulation with np = %d and nel = %d completed',...
                                         Results(1).NumProcessors(iS2),Results(1).NumElements(iS2));
  end
  
  fprintf('\n--------------------------------------------------------------------------------\n');
  
  % Clear variables
  if iS<Simulation.NumSimulations
    clear RefElement Faces Block Mesh;
  end
  
end % End of the single simulation

% Perform test
if matchField(Options,'Test')
  Simulation.TestPassed=eval(Options.Test);
  if Simulation.TestPassed; fprintf('\nTEST PASSED\n');
  else;                     fprintf('\nTEST FAILED\n'); end
end

% Hold off potential figures
if not(isempty(findobj('type','figure')))
  hold off
end

% Clear useless variables
if not(matchField(Simulation,'ComputeOnlyPreProcessing','yes'))
  ClearVariables=[who('i*'); who('Plot*'); who('-regexp','Aux$'); who('-regexp','Par$')];
  clear(ClearVariables{:});
end

% Save results
if matchField(Options,'SaveResults','yes')
  save(sprintf('%s/output/%s%s.mat',fileparts(which(mfilename)),Options.SaveResultsFolder,...
                  FileName),'-regexp','^(?!(CPUTime|Passed|Test)$).',Options.SaveResultsMatVersion);
  fprintf('\nSaved output file in: output/%s%s.mat',Options.SaveResultsFolder,FileName);
  fprintf('\nSaved  diary file in: output/%s%s.txt\n',Options.SaveResultsFolder,FileName);
  diary(sprintf('%s/output/%s%s.txt',fileparts(which(mfilename)),Options.SaveResultsFolder,...
                                                                              FileName)); diary off;
end