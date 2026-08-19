%% Pre-processing

% Store the input deck before any of its fields is overwritten
Input=cell2struct(cellfun(@eval,who,'UniformOutput',false),who);

% Get the name of the input file
FileName=dbstack(1).name;

% Start the diary
if matchField(Options,'SaveResults','yes')
  if not(matchField(Options,'SaveResultsFolder'))
    Options.SaveResultsFolder='';
  end
  if not(matchField(Options,'SaveResultsMatVersion'))
    Options.SaveResultsMatVersion='';
  end
  diary off;
  diary(sprintf('%s/output/%s%s.txt',...
    fileparts(which(mfilename)),Options.SaveResultsFolder,FileName));
end

% Prevent the parallel pool from starting automatically
ParallelSettings=parallel.Settings;
ParallelSettings.Pool.AutoCreate=false;

% Set the working precision (multiprecision arithmetic beyond 16 digits)
if matchField(Simulation,'Digits') && Simulation.Digits~=16
  mp.Digits(Simulation.Digits);
else
  Simulation.Digits=16;
  if contains(path,'advanpix')
    rmpath(extractBefore(path,strfind(path,'advanpix')+8));
  end
end
[Parameters(:).MP]=deal(MP);

% Create the formulation object of each discretization
nD=numel(Parameters);
Simulation.NumDiscretizations=nD;
Formulation=cell(1,nD);
for iD=1:nD
  Formulation{iD}=eval(Parameters(iD).Formulation);
end
Simulation.Domain=Formulation{1}.Domain;

% Get the number of simulations
switch Simulation.Type
  case 'SingleSimulation'
    Simulation.NumSimulations=1;
  case 'ConvergenceSpace'
    if exist('Mesh','var') && matchField(Mesh,'File')
      NumMeshes=size(Mesh.File,2);
      EmptyMesh=strcmp(Mesh.File,'None');
    elseif exist('Distmesh','var')
      NumMeshes=size(Distmesh.ElementSize,2);
      EmptyMesh=isnan(Distmesh.ElementSize);
    else
      NumMeshes=size(Mesh.MinElementSize,2);
      EmptyMesh=isnan(Mesh.MinElementSize);
    end
    Simulation.NumSimulations=numel(Parameters(1).Degree)*NumMeshes-sum(EmptyMesh(:));
  case 'ConvergenceTime'
    Simulation.NumSimulations=numel(Time.BDFOrder)*size(Time.TimeStepSize,2)...
      -sum(isnan(Time.TimeStepSize(:)));
  case 'ParametricStudy'
    Simulation.NumSimulations=numel(Parameters(1).(Simulation.ParameterStudy));
  case 'ScalingStrong'
    Simulation.NumSimulations=numel(Simulation.NumProcessors);
  case 'ScalingWeak'
    Simulation.NumSimulations=(numel(Simulation.NumProcessors)==size(Mesh.File,2))...
      *numel(Simulation.NumProcessors);
end

% Initialize the results
Results=struct();

% Run each simulation
iS1=1;
iS2=0;
for iS=1:Simulation.NumSimulations

  % Set the parameters of the current simulation
  switch Simulation.Type
    case 'SingleSimulation'
      iS1=1;
      iS2=1;
    case 'ConvergenceSpace'
      iS2=iS2+1;
      if iS2>NumMeshes || (size(EmptyMesh,1)>1 && EmptyMesh(iS1,iS2))
        iS1=iS1+1;
        iS2=1;
      end
      for iD=1:nD
        Parameters(iD).Degree=Input.Parameters(iD).Degree(iS1);
        if matchField(Parameters,'MeshRefinement')
          Parameters(iD).MeshRefinement=Input.Parameters(iD).MeshRefinement(iS2);
        end
      end
      if exist('Distmesh','var')
        Distmesh.ElementSize=Input.Distmesh.ElementSize(min(end,iS1),min(end,iS2));
      elseif not(matchField(Input,'Mesh') && matchField(Input.Mesh,'File'))
        Mesh.MinElementSize=Input.Mesh.MinElementSize(min(end,iS1),min(end,iS2));
        Mesh.MaxElementSize=Input.Mesh.MaxElementSize(min(end,iS1),min(end,iS2));
        Mesh.MeshGradation=Input.Mesh.MeshGradation;
      end
      if strcmp(Time.TimeDependent,'yes')
        Time.TimeStepSize=Input.Time.TimeStepSize(min(end,iS1),min(end,iS2));
      end
    case 'ConvergenceTime'
      iS2=iS2+1;
      if iS2>size(Input.Time.TimeStepSize,2) || (size(Input.Time.TimeStepSize,1)>1 && ...
                                                 isnan(Input.Time.TimeStepSize(iS1,iS2)))
        iS1=iS1+1;
        iS2=1;
      end
      Time.BDFOrder=Input.Time.BDFOrder(iS1);
      Time.TimeStepSize=Input.Time.TimeStepSize(min(end,iS1),min(end,iS2));
      if not(matchField(Input,'Mesh') && matchField(Input.Mesh,'File'))
        Mesh.MinElementSize=Input.Mesh.MinElementSize(min(end,iS1),min(end,iS2));
        Mesh.MaxElementSize=Input.Mesh.MaxElementSize(min(end,iS1),min(end,iS2));
        Mesh.MeshGradation=Input.Mesh.MeshGradation;
      end
    case 'ParametricStudy'
      iS2=iS2+1;
      for iD=1:nD
        Parameters(iD).(Simulation.ParameterStudy)=...
          Input.Parameters(iD).(Simulation.ParameterStudy)(iS2);
      end
    case {'ScalingStrong','ScalingWeak'}
      iS2=iS2+1;
      Simulation.NumProcessors=Input.Simulation.NumProcessors(iS2);
  end

  % Select the boundaries of the current simulation
  if any(strcmp(Simulation.Type,{'ConvergenceSpace','ScalingWeak'}))
    for iD=1:nD
      BoundaryNames=fieldnames(Input.Boundaries(iD));
      if size(Input.Boundaries(iD).(BoundaryNames{1}),1)>1
        for iBN=1:numel(BoundaryNames)
          if not(isempty(Input.Boundaries(iD).(BoundaryNames{iBN})))
            Boundaries(iD).(BoundaryNames{iBN})=...
              Input.Boundaries(iD).(BoundaryNames{iBN})(iS2,:); %#ok
          end
        end
      end
    end
  end

  % Resize the parallel pool if needed
  Pool=gcp; if isempty(Pool); Pool=[]; Pool.NumWorkers=1; end
  if matchField(Simulation,'NumProcessors') && Simulation.NumProcessors~=Pool.NumWorkers
    if Pool.NumWorkers>1
      fprintf('\n'); delete(Pool);
    end
    if Simulation.NumProcessors>1
      fprintf('\n'); parpool(parcluster('local'),Simulation.NumProcessors);
    end
  end

  fprintf('\n--------------------------------------------------------------------------------')
  fprintf('\nInput file: %s\n',FileName);

  % Load a stored pre-processing or perform it
  if matchField(Simulation,'SkipPreProcessing','yes')
    fprintf('\nPre-processing  file: output/%s%s.mat',...
      Options.SaveResultsFolder,Simulation.PreProcessingFileName);
    load(sprintf('%s/output/%s%s.mat',fileparts(which(mfilename)),...
      Options.SaveResultsFolder,Simulation.PreProcessingFileName),...
      '-regexp','^(?!(FileName|Simulation)$).');
  else
    TimerPreProcessing=tic;
    fprintf('\nPre-processing  started...');

    % Create the model of each discretization
    for iD=1:nD
      Model(iD)=createpde(); %#ok
    end

    % Get the geometry and the linear mesh
    if matchField(Input,'Mesh') && matchField(Input.Mesh,'File')
      load(['geometry/',Input.Mesh.File{min(end,iS1),min(end,iS2)},'.mat']);
    end
    if exist('Mesh','var') && matchField(Mesh,'Nodes') && matchField(Mesh,'Elements')
      for iD=1:nD
        Mesh(iD).DegreeK=Mesh(iD);
        NumVertexNodes=size(Mesh(iD).Nodes,1)+1;
        Mesh(iD).Nodes=Mesh(iD).Nodes(:,1:max(max(Mesh(iD).Elements(1:NumVertexNodes,:))));
        Mesh(iD).Elements=Mesh(iD).Elements(1:NumVertexNodes,:);
        [Geometry(iD),Mesh(iD).Degree1]=geometryFromMesh(Model(iD),...
          double(Mesh(iD).Nodes),Mesh(iD).Elements); %#ok
        [~,Order]=ismember(Mesh(iD).Degree1.Elements',...
          Mesh(iD).DegreeK.Elements(1:NumVertexNodes,:)','rows');
        Mesh(iD).DegreeK.Elements=Mesh(iD).DegreeK.Elements(:,Order);
      end
    elseif exist('Distmesh','var')
      Mesh=struct();
      for iD=1:nD
        [Mesh(iD).Nodes,Mesh(iD).Elements]=distmesh(Distmesh(iD).DistanceFunction,...
          Distmesh(iD).ElementSizeFunction,Distmesh(iD).ElementSize,...
          Distmesh(iD).BoundingBox,Distmesh(iD).FixedNodePositions);
        [Geometry(iD),Mesh(iD).Degree1]=geometryFromMesh(Model(iD),...
          Mesh(iD).Nodes',Mesh(iD).Elements');
      end
    elseif exist('Geometry','var') && (isnumeric(Geometry) || ...
           isa(Geometry,'pde.DiscreteGeometry') || isa(Geometry,'pde.AnalyticGeometry'))
      for iD=1:nD
        if isnumeric(Geometry)
          Model(iD).Geometry=geometryFromEdges(Model(iD),Geometry(:,:,iD));
        else
          Model(iD).Geometry=Geometry(iD);
        end
        Mesh(iD).Degree1=generateMesh(Model(iD),'GeometricOrder','Linear',...
          'Hmin',Mesh(iD).MinElementSize,'Hmax',Mesh(iD).MaxElementSize,...
          'Hgrad',Mesh(iD).MeshGradation);
      end
      if isnumeric(Geometry)
        clear Geometry;
        for iD=1:nD
          Geometry(iD)=Model(iD).Geometry; %#ok
        end
      end
    end

    % Get the space dimension and the discretization type
    Sizes=struct();
    for iD=1:nD
      Sizes(iD).NumSpaceDim=1+(Geometry(iD).NumFaces>0)+(Geometry(iD).NumCells>0);
      Parameters(iD).DiscretizationType=Formulation{iD}.DiscretizationType;
    end

    % Assign the default uniform nodes distribution
    if not(matchField(Parameters,'NodesDistribution'))
      [Parameters(:).NodesDistribution]=deal('Uniform');
    end

    % Create the reference elements
    for iD1=1:nD
      for iD2=1:nD
        RefElement(iD1,iD2)=createReferenceElement(Sizes(iD1).NumSpaceDim,...
          Parameters(iD1).Degree,max(Parameters(iD1).Degree,Parameters(iD2).Degree),...
          Parameters(iD1).NodesDistribution,Simulation.Digits); %#ok
      end
    end
    for iD=1:nD
      RefElement(iD,iD).Degree1=createReferenceElement(Sizes(iD).NumSpaceDim,1,1,...
        Parameters(iD).NodesDistribution,Simulation.Digits);
      RefElement(iD,iD).MeshDegree1=createReferenceElement(Sizes(iD).NumSpaceDim,...
        1,Parameters(iD).Degree,Parameters(iD).NodesDistribution,Simulation.Digits);
      RefElement(iD,iD).MeshDegreeK=createReferenceElement(Sizes(iD).NumSpaceDim,...
        Parameters(iD).Degree,Parameters(iD).Degree,...
        Parameters(iD).NodesDistribution,Simulation.Digits);
      if strcmp(Parameters(iD).DiscretizationType,'HDG') && ...
         strcmp(Parameters(iD).PostProcessingHDG,'yes')
        RefElement(iD,iD).PostLow=createReferenceElement(Sizes(iD).NumSpaceDim,...
          Parameters(iD).Degree,Parameters(iD).Degree+1,...
          Parameters(iD).NodesDistribution,Simulation.Digits);
        RefElement(iD,iD).Post=createReferenceElement(Sizes(iD).NumSpaceDim,...
          Parameters(iD).Degree+1,Parameters(iD).Degree+1,...
          Parameters(iD).NodesDistribution,Simulation.Digits);
      end
    end

    % Get the number of components
    for iD=1:nD
      Sizes(iD).NumGlobalComp=Formulation{iD}.NumGlobalComp(Sizes(iD).NumSpaceDim);
      if strcmp(Parameters(iD).DiscretizationType,'HDG')
        Sizes(iD).NumVoigtComp=Formulation{iD}.NumVoigtComp(Sizes(iD).NumSpaceDim);
        Sizes(iD).NumLocalComp=Formulation{iD}.NumLocalComp(Sizes(iD).NumSpaceDim);
        if strcmp(Parameters(iD).PostProcessingHDG,'yes')
          Sizes(iD).NumPostComp=Formulation{iD}.NumPostComp(Sizes(iD).NumSpaceDim);
        end
      end
    end

    % Get the faces and the boundary conditions
    if not(exist('Regions','var'))
      Regions=[];
    end
    [Faces,Elements,BCs]=getFacesAndBCs(Simulation,Geometry,Mesh,Boundaries,Regions,...
      Parameters,RefElement,Sizes);

    % Get the degree-independent sizes
    for iD=1:nD
      Sizes(iD).NumElements=size(Mesh(iD).Degree1.Elements,2);
      Sizes(iD).NumElementFaces=size(RefElement(iD,iD).Degree1.FaceNodesElem,1);
    end

    % Build the computational mesh and the corresponding sizes
    for iD=1:nD

      % Generate the high-order mesh
      if Parameters(iD).Degree<=1
        Mesh(iD).Nodes=Mesh(iD).Degree1.Nodes;
        Mesh(iD).Elements=Mesh(iD).Degree1.Elements;
      else
        if not(matchField(Mesh,'DegreeK')) || ...
           size(Mesh(iD).DegreeK.Elements,1)==size(Mesh(iD).Degree1.Elements,1)
          Mesh(iD).DegreeK=generateMeshHighOrder(Mesh(iD).Degree1,...
            RefElement(iD,iD).MeshDegree1,RefElement(iD,iD).MeshDegreeK);
        end
        Mesh(iD).Nodes=Mesh(iD).DegreeK.Nodes;
        Mesh(iD).Elements=Mesh(iD).DegreeK.Elements;
      end
      Mesh(iD).MaxElementSize=Mesh(iD).Degree1.MaxElementSize;
      Mesh(iD).MinElementSize=Mesh(iD).Degree1.MinElementSize;
      Mesh(iD).MeshGradation =Mesh(iD).Degree1.MeshGradation;

      % Generate the mesh for the HDG post-processing
      HDG=strcmp(Parameters(iD).DiscretizationType,'HDG');
      HDGPost=HDG && strcmp(Parameters(iD).PostProcessingHDG,'yes');
      if HDGPost
        MeshPost=generateMeshHighOrder(Mesh(iD),RefElement(iD,iD).PostLow,...
                                                RefElement(iD,iD).Post);
        Mesh(iD).Post.Nodes=MeshPost.Nodes;
        Mesh(iD).Post.Elements=MeshPost.Elements;
      end

      % Get the mesh-dependent sizes
      Sizes(iD).NumNodes=size(Mesh(iD).Nodes,2);
      Sizes(iD).NumElementNodes=size(RefElement(iD,iD).NodesCoordElem,1);
      Sizes(iD).NumFaceNodes=size(RefElement(iD,iD).NodesCoordFace,1);
      if Sizes(iD).NumSpaceDim==2
        Sizes(iD).NumEdgeNodes=Sizes(iD).NumFaceNodes;
      elseif Sizes(iD).NumSpaceDim==3
        Sizes(iD).NumEdgeNodes=size(RefElement(iD,iD).ShapeFunctionsEdge,2);
      end
      if HDGPost
        Sizes(iD).NumElementNodesPost=size(RefElement(iD,iD).Post.NodesCoordElem,1);
        Sizes(iD).NumFaceNodesPost=size(RefElement(iD,iD).Post.NodesCoordFace,1);
      end

      % Number of nodes stored per element in the mesh connectivity
      NumMeshNodes=max(Sizes(iD).NumSpaceDim+1,Sizes(iD).NumElementNodes);

      % Duplicate the nodes shared by neighboring elements for HDG
      if HDG
        Mesh(iD).Nodes=Mesh(iD).Nodes(:,Mesh(iD).Elements(1:end));
        Mesh(iD).Elements=reshape(1:NumMeshNodes*Sizes(iD).NumElements,...
                                   [NumMeshNodes,Sizes(iD).NumElements]);
        if HDGPost
          Mesh(iD).Post.Nodes=Mesh(iD).Post.Nodes(:,Mesh(iD).Post.Elements(1:end));
          Mesh(iD).Post.Elements=reshape(...
            1:Sizes(iD).NumElementNodesPost*Sizes(iD).NumElements,...
            [Sizes(iD).NumElementNodesPost,Sizes(iD).NumElements]);
        end
      end

      % Add the z-component in 2D
      if Sizes(iD).NumSpaceDim==2
        Mesh(iD).Nodes=[Mesh(iD).Nodes;zeros(1,size(Mesh(iD).Nodes,2))];
        if HDGPost
          Mesh(iD).Post.Nodes=[Mesh(iD).Post.Nodes;zeros(1,size(Mesh(iD).Post.Nodes,2))];
        end
      end

      % Add the elements' centroids for FCFV
      if HDG && Parameters(iD).Degree==0
        Mesh(iD).Centroids=squeeze(mean(reshape(...
          Mesh(iD).Nodes(:,Mesh(iD).Elements(1:Sizes(iD).NumSpaceDim+1,:)),...
          3,Sizes(iD).NumSpaceDim+1,[]),2));
      end

      % Store the initial mesh
      Mesh(iD).NodesInitial=Mesh(iD).Nodes;

      % Store the element data for an effective use of parfor
      Elements(iD).Nodes=mat2cell(Mesh(iD).Nodes(:,Mesh(iD).Elements(:)),...
        3,ones(Sizes(iD).NumElements,1)*NumMeshNodes)';
      Elements(iD).Element=mat2cell(Mesh(iD).Elements,...
        NumMeshNodes,ones(Sizes(iD).NumElements,1))';
      Elements(iD).GlobalLocal=Faces(iD,iD).GlobalLocal;
      Elements(iD).NodesInitial=mat2cell(Mesh(iD).NodesInitial(:,Mesh(iD).Elements(:)),...
        3,ones(Sizes(iD).NumElements,1)*NumMeshNodes)';

      % Get the number of degrees of freedom
      if HDG
        Sizes(iD).NumFaces=max(Faces(iD,iD).Connectivity(:));
        Sizes(iD).NumGlobalNodes=Sizes(iD).NumFaces*Sizes(iD).NumFaceNodes;
        Sizes(iD).NumLocalNodes=Sizes(iD).NumElements*Sizes(iD).NumElementNodes;
        Sizes(iD).NumElementRhsCoef=Sizes(iD).NumGlobalComp...
          *Sizes(iD).NumElementFaces*Sizes(iD).NumFaceNodes;
      else
        Sizes(iD).NumGlobalNodes=Sizes(iD).NumNodes;
        Sizes(iD).NumElementRhsCoef=Sizes(iD).NumGlobalComp*Sizes(iD).NumElementNodes;
      end
      Sizes(iD).NumElementLhsCoef=Sizes(iD).NumElementRhsCoef^2;
    end

    % Get the sizes of the interface contributions
    if matchField(Boundaries,'Interface')
      for iD1=1:nD
        for iD2=1:nD
          Sizes(iD1).NumFacesInterface(iD2)=size(Faces(iD1,iD2).Interface,1);
          Sizes(iD1).NumElementLhsCoupCoef(iD2)=Sizes(iD1).NumElementRhsCoef...
                                               *Sizes(iD2).NumElementRhsCoef;
        end
      end
    end

    % Define the pattern of the monolithic matrix
    Block=defineMatrixPattern(Simulation,Parameters,Mesh,Faces,Sizes);

    % Get the fixed and the periodic degrees of freedom
    Block=getFixedDofs(Block,Simulation,Mesh,Boundaries,Parameters,RefElement,Sizes);
    Block=getPeriodicDofs(Block,Simulation,Mesh,Boundaries,Parameters,RefElement,Sizes);
    System.DOFsFixed=[]; System.DOFsMaster=[]; System.DOFsSlave=[];
    for iD=1:nD
      Offset=sum(arrayfun(@(S) S.NumGlobalNodes*S.NumGlobalComp,Sizes(1:iD-1)));
      System.DOFsFixed =[System.DOFsFixed ,Block(iD,iD).DOFsFixed +Offset];
      System.DOFsMaster=[System.DOFsMaster,Block(iD,iD).DOFsMaster+Offset];
      System.DOFsSlave =[System.DOFsSlave ,Block(iD,iD).DOFsSlave +Offset];
    end

    % Set the default settings of the nonlinear solver
    if strcmp(System.Nonlinear,'no')
      System.Tolerance=1e+12;
      System.MaxIterations=1;
    end

    % Set the default settings of the time integration
    if strcmp(Time.TimeDependent,'yes')
      for iD=1:nD
        Parameters(iD).TimeDerOrder=Formulation{iD}.TimeDerOrder;
      end
      if not(matchField(Parameters,'PredictorDegree'))
        [Parameters(:).PredictorDegree]=deal(0);
      end
      Time.NumTimeSteps=double(round((Time.FinalTime-Time.InitialTime)/Time.TimeStepSize));
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
    else
      Time.InitialTime=0;
      Time.FinalTime=0;
      Time.TimeStepSize=0;
      Time.NumTimeSteps=1;
    end
    Time.Time=Time.InitialTime;

    % Set the plotting styles and the axis limits
    Colors={'b','r','k','g','m','c','y','b','r','k','g','m','c','y'};
    Styles={'-','--',':','-.','-','--',':','-.','-','--',':','-.'};
    Markers={'o','s','d','p','^','x','*','h','o','s','d','p'};
    Limits=minmax(double([Mesh.Nodes]));
    Limits=Limits+[-1,+1].*(Limits(:,1)==Limits(:,2));
    Wide=abs(diff(Limits(1,:)))>abs(diff(Limits(2,:)));

    % Plot the geometry
    if strcmp(Options.PlotGeometry,'yes')
      figure('Color','w');
      for iD=1:nD
        subplot(1+(nD-1)*Wide,1+(nD-1)*not(Wide),iD)
        if Sizes(iD).NumSpaceDim==2
          Plots=pdegplot(Geometry(iD),'EdgeLabels','on');
          set(Plots(1),'Color',Colors{iD})
        elseif Sizes(iD).NumSpaceDim==3
          Plots=pdegplot(Geometry(iD),'FaceLabels','on','FaceAlpha',0.25);
          set(Plots(1),'FaceColor',Colors{iD})
          Children=get(gca,'Children'); delete(Children([4,5]));
        end
        grid on; box on; axis equal;
        title(['Geometry',repmat(sprintf(' (%d)',iD),1,nD>1)]);
        xlim(Limits(1,:)); xlabel('x');
        ylim(Limits(2,:)); ylabel('y');
        zlim(Limits(3,:)); zlabel('z');
        pause(eps);
      end
    end

    % Plot the mesh and the boundaries
    if strcmp(Options.PlotMesh,'yes')
      figure('Color','w');
      for iD=1:nD
        subplot(1+(nD-1)*Wide,1+(nD-1)*not(Wide),iD)
        if Sizes(iD).NumSpaceDim==2
          pdeplot(Mesh(iD).Degree1.Nodes,Mesh(iD).Degree1.Elements,'EdgeColor',Colors{iD});
        elseif Sizes(iD).NumSpaceDim==3
          pdeplot3D(Mesh(iD).Degree1.Nodes,Mesh(iD).Degree1.Elements,...
            'FaceColor','w','EdgeColor',Colors{iD},'FaceAlpha',1);
          Children=get(gca,'Children'); delete(Children([2,3]));
        end
        hold on
        BoundaryNames=fieldnames(BCs(iD));
        Plots=matlab.graphics.chart.primitive.Line.empty();
        for iBN=1:numel(BoundaryNames)
          NodesBC=Mesh(iD).Degree1.Nodes(:,BCs(iD).(BoundaryNames{iBN}));
          if Sizes(iD).NumSpaceDim==2
            PlotBC=plot(NodesBC(1,:),NodesBC(2,:),...
              'ok','MarkerFaceColor',Colors{iBN},'DisplayName',BoundaryNames{iBN});
          elseif Sizes(iD).NumSpaceDim==3
            PlotBC=plot3(NodesBC(1,:),NodesBC(2,:),NodesBC(3,:),...
              'ok','MarkerFaceColor',Colors{iBN},'DisplayName',BoundaryNames{iBN});
          end
          if not(isempty(PlotBC))
            Plots(iBN)=PlotBC;
          end
        end
        hold off; box on; axis equal;
        title(['Mesh',repmat(sprintf(' (%d)',iD),1,nD>1)]);
        legend(Plots(arrayfun(@(i) isa(Plots(i),'matlab.graphics.chart.primitive.Line'),...
          1:numel(Plots))),'Location','SouthEast')
        set(gca,'Visible','on')
        xlim(Limits(1,:)); xlabel('x');
        ylim(Limits(2,:)); ylabel('y');
        zlim(Limits(3,:)); zlabel('z');
        pause(eps);
      end
    end

    % Plot the mesh distortion
    if matchField(Options,'PlotMeshDistortion','yes')
      NumColors=10;
      ColorMap=jet(NumColors);
      figure('Color','w');
      [Mesh(:).Distortion]=deal([]);
      for iD=1:nD
        Mesh(iD).Distortion=1-meshQuality(Mesh(iD).Degree1)';
        for iCol=1:NumColors
          iElems=find(Mesh(iD).Distortion>(iCol-1)/NumColors & ...
                      Mesh(iD).Distortion< iCol   /NumColors);
          if not(isempty(iElems))
            if Sizes(iD).NumSpaceDim==2
              pdeplot(Mesh(iD).Degree1.Nodes,Mesh(iD).Degree1.Elements(:,iElems),...
                'EdgeColor',ColorMap(iCol,:));
            elseif Sizes(iD).NumSpaceDim==3
              pdeplot3D(Mesh(iD).Degree1.Nodes,Mesh(iD).Degree1.Elements(:,iElems),...
                'FaceColor',ColorMap(iCol,:),'EdgeColor','k','FaceAlpha',1);
              Children=get(gca,'Children'); delete(Children([2,3]));
            end
          end
          hold on
        end
      end
      hold off; box on; axis equal; colormap(ColorMap); caxis([0,1]);
      colorbar('YLim',[0,1],'YTick',linspace(0,1,NumColors+1))
      title('Mesh distortion');
      set(gca,'Visible','on')
      xlim(Limits(1,:)); xlabel('x');
      ylim(Limits(2,:)); ylabel('y');
      zlim(Limits(3,:)); zlabel('z');
      pause(eps);
    end

    % Report the pre-processing time and memory
    Timer.PreProcessing=toc(TimerPreProcessing);
    Memory.PreProcessing=whos; Memory.PreProcessing=sum(vertcat(Memory.PreProcessing.bytes));
    fprintf('\nPre-processing  completed in %.1f sec (%.1f MB)',...
      Timer.PreProcessing,Memory.PreProcessing/1e6);

    % Stop here if only the pre-processing is requested
    if matchField(Simulation,'ComputeOnlyPreProcessing','yes')
      fprintf('\n'); break;
    end
  end

  %% Processing

  TimerProcessing=tic;
  fprintf('\nProcessing      started...');
  Timer.Evaluation=[];
  Timer.Solution=[];
  Timer.Local=[];

  % Compute the quantity of interest at the start of the simulation
  if matchField(Options,'ComputeQuantityStart')
    eval(Options.ComputeQuantityStart);
  end

  % Initialize the unknowns
  for iD=1:nD
    Block=Formulation{iD}.initializeUnknowns(iD,Block,Parameters,Time,Sizes);
  end

  % Compute the initial conditions
  for iD=1:nD
    Block=Formulation{iD}.computeInitialConditions(iD,Block,Parameters,Mesh,Faces,Time,...
      RefElement,Sizes);
  end

  % Get the first and the last time step
  if matchField(Simulation,'RestartInitialFinalSteps')
    Time.InitialTimeStep=Simulation.RestartInitialFinalSteps(1);
    Time.FinalTimeStep=Simulation.RestartInitialFinalSteps(2);
  else
    Time.InitialTimeStep=1;
    Time.FinalTimeStep=Time.NumTimeSteps;
  end

  if strcmp(Time.TimeDependent,'yes')
    fprintf('\nTime integration ...............................................................');
  end

  % Advance in time
  for iT=Time.InitialTimeStep:Time.FinalTimeStep
    Time.TimeStep=double(iT);

    % Load the restart file
    if matchField(Simulation,'RestartInitialFinalSteps') && ...
       Time.TimeStep==Time.InitialTimeStep && Time.InitialTimeStep>1
      fprintf('\nRestart file:    output/%s%s.mat',...
        Options.SaveResultsFolder,Simulation.RestartFileName);
      load(sprintf('%s/output/%s%s.mat',fileparts(which(mfilename)),...
        Options.SaveResultsFolder,Simulation.RestartFileName),...
        '-regexp','^(?!(FileName|Simulation|Time|Options)$).');
    end

    % Update the time and select the BDF coefficients
    if strcmp(Time.TimeDependent,'yes')
      Time.Time=Time.InitialTime+Time.TimeStep*Time.TimeStepSize;
      fprintf('\nTime = %.2e',Time.Time);
      fprintf(' Step = %.*d/%d',ceil(log10(Time.NumTimeSteps)+eps),...
        Time.TimeStep,Time.NumTimeSteps);
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

    % Predict the unknowns by extrapolation of the previous solutions
    if strcmp(Time.TimeDependent,'yes')
      for iD=1:nD
        kP=min(Parameters(iD).PredictorDegree,size(Block(iD,iD).SolutionOld,3)-1);
        CP=round((-1)^kP*((-kP:1)'.^(0:kP+1)'\[zeros(kP+1,1);factorial(kP+1)]));
        SolutionPred=sum(Block(iD,iD).SolutionOld(:,:,1:kP+1).*reshape(CP(2:end),1,1,kP+1),3);
        if strcmp(Parameters(iD).DiscretizationType,'CG')
          Block(iD,iD).SolutionGlobal=SolutionPred;
        elseif strcmp(Parameters(iD).DiscretizationType,'HDG')
          Block(iD,iD).SolutionLocal=SolutionPred;
        end
      end
    end

    % Initialize the system solutions and evaluate the solution at the fixed DOFs
    for iD=1:nD
      Block(iD,iD).SysSolutionGlobal=MP*zeros(Sizes(iD).NumGlobalNodes*Sizes(iD).NumGlobalComp,1);
      if strcmp(Parameters(iD).DiscretizationType,'HDG')
        Block(iD,iD).SysSolutionLocal=MP*zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumLocalComp);
      end
    end
    for iD=1:nD
      Block=Formulation{iD}.evaluateSolutionFixedDofs(iD,Block,Parameters,Mesh,Time,Sizes);
    end

    % Compute the quantity of interest at the current time step
    if matchField(Options,'ComputeQuantityTimeStep')
      eval(Options.ComputeQuantityTimeStep);
    end

    % Iterate until convergence of the Newton-Raphson method
    System.ResidualNorm=Inf;
    System.Iteration=0;
    while System.ResidualNorm>System.Tolerance && System.Iteration<System.MaxIterations
      System.Iteration=System.Iteration+1;
      TimerEvaluation=tic;

      % Store the element data for an effective use of parfor
      for iD=1:nD
        ElementNodes=Mesh(iD).Elements(1:Sizes(iD).NumElements*Sizes(iD).NumElementNodes);
        NodeBlocks=ones(Sizes(iD).NumElements,1)*Sizes(iD).NumElementNodes;
        if strcmp(Parameters(iD).DiscretizationType,'CG')
          NodalSolution='SolutionGlobal';
          NumNodalComp=Sizes(iD).NumGlobalComp;
        else
          NodalSolution='SolutionLocal';
          NumNodalComp=Sizes(iD).NumLocalComp;
          Elements(iD).SolutionGlobal=mat2cell(...
            Block(iD,iD).SolutionGlobal(cell2mat(Faces(iD,iD).GlobalLocal),:,:),...
            ones(Sizes(iD).NumElements,1)*Sizes(iD).NumElementFaces*Sizes(iD).NumFaceNodes...
            *Sizes(iD).NumGlobalComp,size(Block(iD,iD).SolutionGlobal,3));
        end
        Elements(iD).(NodalSolution)=mat2cell(Block(iD,iD).(NodalSolution)(ElementNodes,:,:),...
          NodeBlocks,NumNodalComp,size(Block(iD,iD).(NodalSolution),3));
        if strcmp(Time.TimeDependent,'yes')
          Elements(iD).SolutionOld=mat2cell(Block(iD,iD).SolutionOld(ElementNodes,:,:),...
            NodeBlocks,NumNodalComp,size(Block(iD,iD).SolutionOld,3));
        else
          Elements(iD).SolutionOld=cell(Sizes(iD).NumElements,1);
        end
      end

      % Compute the quantity of interest at the current iteration
      if matchField(Options,'ComputeQuantityIteration')
        eval(Options.ComputeQuantityIteration);
      end

      % Build the block matrices and residuals
      for iD1=1:nD
        for iD2=1:nD
          if iD1==iD2
            [Block,Elements]=Formulation{iD1}.buildBlock(iD1,Block,Elements,Simulation,...
              Parameters,Mesh,Faces,Time,RefElement,Sizes);
          else
            Block=Formulation{iD1}.doCoupling(iD1,iD2,Block,Elements,Simulation,...
              Parameters,Mesh,Faces,Time,RefElement,Sizes);
          end
        end
      end
      Timer.Evaluation=[Timer.Evaluation,toc(TimerEvaluation)];

      % Assemble the monolithic system
      System.Lhs=struct();
      for iD=1:nD
        System.Lhs(iD).LhsGlobal=horzcat(Block(iD,:).LhsGlobal);
      end
      System.Lhs=vertcat(System.Lhs.LhsGlobal);
      System.Rhs=vertcat(Block.RhsGlobal);

      % Impose the fixed boundary conditions strongly
      System.Lhs(System.DOFsFixed,:)=0;
      System.Lhs(System.DOFsFixed,System.DOFsFixed)=eye(numel(System.DOFsFixed));
      System.Rhs(System.DOFsFixed,:)=0;
      if size(System.Rhs,2)==size(System.Lhs,2)
        System.Rhs(System.DOFsFixed,System.DOFsFixed)=eye(numel(System.DOFsFixed));
      end

      % Impose the periodic boundary conditions strongly
      System.Rhs(System.DOFsMaster,:)=System.Rhs(System.DOFsMaster,:)...
                                     +System.Rhs(System.DOFsSlave,:);
      System.Lhs(System.DOFsMaster,:)=System.Lhs(System.DOFsMaster,:)...
                                     +System.Lhs(System.DOFsSlave,:);
      System.Rhs(System.DOFsSlave,:)=0;
      System.Lhs(System.DOFsSlave,:)=0;
      System.Lhs(System.DOFsSlave,System.DOFsSlave) =+eye(numel(System.DOFsSlave));
      System.Lhs(System.DOFsSlave,System.DOFsMaster)=-eye(numel(System.DOFsSlave));

      % Symmetrize the matrix
      if matchField(System,'SymmetrizeMatrix','yes') && ...
         (not(matchField(System,'SymmetrizeMatrixOnlyOnce')) || ...
          (matchField(System,'SymmetrizeMatrixOnlyOnce','yes') && ...
           Time.TimeStep<=Time.BDFOrderEff && System.Iteration==1))
        System.Lhs=(System.Lhs+System.Lhs.')/2;
        if size(System.Rhs,2)==size(System.Lhs,2)
          System.Rhs=(System.Rhs+System.Rhs.')/2;
        end
      end

      % Equilibrate the system
      if matchField(Solver,'Equilibrate','yes')
        [System.EquilP,System.EquilR,System.EquilC]=equilibrate(System.Lhs);
        System.Lhs=System.EquilR*System.EquilP*System.Lhs*System.EquilC;
        System.Rhs=System.EquilR*System.EquilP*System.Rhs;
      end

      % Compute the preconditioner
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

      % Distribute the system among the workers
      if matchField(System,'Distribute','yes')
        System.Lhs=distributed(System.Lhs);
        System.Rhs=distributed(System.Rhs);
      end

      TimerSolution=tic;

      % Decompose the matrix once and reuse the decomposition
      if matchField(System,'Decompose','yes')
        if Time.TimeStep<Time.BDFOrder
          System.LhsDecomposed=System.Lhs;
        elseif Time.TimeStep==Time.BDFOrder
          System.LhsDecomposed=decomposition(System.Lhs);
        end
        System.Lhs=System.LhsDecomposed;
      end

      % Store the old solution for the relaxation
      if matchField(System,'RelaxationParameter')
        if System.Iteration>1
          System.SolutionOld=System.Solution;
        else
          System.SolutionOld=sparse(size(System.Rhs,1),1);
        end
      end

      % Solve the global problem
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
          System.Solution=gmres(System.Lhs,System.Rhs,Solver.Restart,Solver.Tolerance,...
            min(Solver.MaxIterations,size(System.Rhs,1)),System.PrecondL,System.PrecondU);
        case 'eigs'
          [Eigenvectors,Eigenvalues]=eigs(System.Lhs,System.Rhs,max(Time.Mode),'smallestabs');
          Eigenvectors=Eigenvectors(:,Time.Mode); [~,iMax]=max(abs(Eigenvectors));
          System.Solution(:,1,:)=Eigenvectors./...
            arrayfun(@(n) Eigenvectors(iMax(n),n),1:numel(Time.Mode));
          Results.NaturalFrequency(1,1,:)=sqrt(diag(Eigenvalues(Time.Mode,Time.Mode)))';
      end
      System.Solution=full(System.Solution);
      Timer.Solution=[Timer.Solution,toc(TimerSolution)];

      % Relax the solution
      if matchField(System,'RelaxationParameter')
        System.Solution=System.RelaxationParameter*System.Solution...
                    +(1-System.RelaxationParameter)*System.SolutionOld;
      end

      % Gather the solution from the workers
      if matchField(System,'Distribute','yes')
        System.Solution=gather(System.Solution);
      end

      % Undo the equilibration
      if matchField(Solver,'Equilibrate','yes')
        System.Solution=System.EquilC*System.Solution;
      end

      % Extract, reshape and update the global solution of each discretization
      for iD=1:nD
        Offset=sum(arrayfun(@(S) S.NumGlobalNodes*S.NumGlobalComp,Sizes(1:iD-1)));
        Block(iD,iD).SysSolutionGlobal=...
          System.Solution(Offset+(1:size(Block(iD,iD).RhsGlobal,1)),1,:);
        if strcmp(Parameters(iD).DiscretizationType,'CG')
          Block(iD,iD).SysSolutionGlobal=reshape(Block(iD,iD).SysSolutionGlobal,...
            Sizes(iD).NumGlobalNodes,Sizes(iD).NumGlobalComp,[]);
        end
        if matchField(System,'IncrementalForm','no')
          Block(iD,iD).SolutionGlobal=Block(iD,iD).SysSolutionGlobal;
        else
          Block(iD,iD).SolutionGlobal=Block(iD,iD).SolutionGlobal+Block(iD,iD).SysSolutionGlobal;
        end
      end

      % Solve the local problems of each HDG discretization
      TimerLocal=tic;
      for iD=1:nD
        if strcmp(Parameters(iD).DiscretizationType,'HDG')
          Elements(iD).SysSolutionGlobal=mat2cell(...
            Block(iD,iD).SysSolutionGlobal(cell2mat(Faces(iD,iD).GlobalLocal),:,:),...
            ones(Sizes(iD).NumElements,1)*Sizes(iD).NumElementFaces*Sizes(iD).NumFaceNodes...
            *Sizes(iD).NumGlobalComp,size(Block(iD,iD).SysSolutionGlobal,3));
          for iElem=1:Sizes(iD).NumElements
            Block(iD,iD).SysSolutionLocal(...
              (iElem-1)*Sizes(iD).NumElementNodes+(1:Sizes(iD).NumElementNodes),:)=...
              reshape(Elements(iD).VecLocal{iElem}...
                     -Elements(iD).MatLocal{iElem}*Elements(iD).SysSolutionGlobal{iElem},...
                      Sizes(iD).NumElementNodes,Sizes(iD).NumLocalComp);
          end
          if matchField(System,'IncrementalForm','no')
            Block(iD,iD).SolutionLocal=Block(iD,iD).SysSolutionLocal;
          else
            Block(iD,iD).SolutionLocal=Block(iD,iD).SolutionLocal+Block(iD,iD).SysSolutionLocal;
          end
        end
      end
      Timer.Local=[Timer.Local,toc(TimerLocal)];

      % Compute the norm of the residual and print the iteration info
      System.ResidualNorm=norm(System.Rhs,'fro');
      if strcmp(System.Nonlinear,'yes')
        HasLocal=any(strcmp({Parameters.DiscretizationType},'HDG'));
        fprintf('\n%2sIter = %d','',System.Iteration);
        fprintf('%*s|Res| = %.1e',7-numel(num2str(System.Iteration)),'',System.ResidualNorm);
        fprintf('%6s[te,ts%s] = [%.1f,%.1f%s]','',repmat(',tl',1,HasLocal),...
          Timer.Evaluation(end),Timer.Solution(end),...
          repmat(sprintf(',%.1f',Timer.Local(end)),1,HasLocal));
      end

      % Compute the quantity of interest at the end of the current iteration
      if matchField(Options,'ComputeQuantityIterationEnd')
        eval(Options.ComputeQuantityIterationEnd);
      end
    end

    % Shift the old solutions
    if strcmp(Time.TimeDependent,'yes') && Time.TimeStep<Time.NumTimeSteps
      Shift=min(Time.BDFOrderEff+1,Time.BDFOrder);
      for iD=1:nD
        Block(iD,iD).SolutionOld(:,:,2:Shift+(Parameters(iD).TimeDerOrder-1))=...
        Block(iD,iD).SolutionOld(:,:,1:Shift-(2-Parameters(iD).TimeDerOrder));
        if strcmp(Parameters(iD).DiscretizationType,'CG')
          Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionGlobal;
        elseif strcmp(Parameters(iD).DiscretizationType,'HDG')
          Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionLocal;
        end
      end
    end

    % Store the results of the current time step
    if matchField(Options,'StoreTimeSteps')
      iST=find(Time.TimeStep==eval(Options.StoreTimeSteps));
      if not(isempty(iST))
        for iD=1:nD
          Results=Formulation{iD}.storeResults(iD,iST,Results,Block,Simulation,Parameters,Mesh,...
            Time,RefElement,Sizes);
        end
      end
    end

    % Compute the quantity of interest
    if matchField(Options,'ComputeQuantity')
      eval(Options.ComputeQuantity);
    end

    % Plot the solution of the current time step
    if matchField(Options,'PlotSolution') && matchField(Options,'PlotSolutionTimeSteps') && ...
       any(Time.TimeStep==eval(Options.PlotSolutionTimeSteps))
      iPlotT=find(Time.TimeStep==eval(Options.StoreTimeSteps));
      plotSolution(Simulation,Mesh,Faces,Results,Parameters,RefElement,Sizes,Options,iPlotT);
    end

    % Export the solution of the current time step to Paraview
    if matchField(Options,'Export2Paraview') && matchField(Options,'StoreTimeSteps') && ...
       matchField(Options,'Export2ParaviewTimeSteps') && ...
       any(Time.TimeStep==eval(Options.Export2ParaviewTimeSteps))
      for iD=1:nD
        export2Paraview(Mesh(iD),Results(iD),Parameters(iD),Time,Sizes(iD),Options,...
          sprintf('%s%s_step_%.*d',FileName,repmat(sprintf('_discr_%d',iD),1,nD>1),...
                  ceil(log10(Time.NumTimeSteps)+eps),Time.TimeStep));
      end
    end

    % Save the results of the current time step
    if matchField(Options,'SaveResults','yes') && matchField(Options,'StoreTimeSteps') && ...
       matchField(Options,'SaveResultsTimeSteps') && ...
       any(Time.TimeStep==eval(Options.SaveResultsTimeSteps))
      save(sprintf('%s/output/%s%s_step_%.*d',fileparts(which(mfilename)),...
        Options.SaveResultsFolder,FileName,ceil(log10(Time.NumTimeSteps)+eps),Time.TimeStep),...
        'Results',Options.SaveResultsMatVersion);
    end
  end

  if strcmp(Time.TimeDependent,'yes')
    fprintf('\n................................................................................');
  end

  % Store the element data for an effective use of parfor
  for iD=1:nD
    ElementNodes=Mesh(iD).Elements(1:Sizes(iD).NumElements*Sizes(iD).NumElementNodes);
    NodeBlocks=ones(Sizes(iD).NumElements,1)*Sizes(iD).NumElementNodes;
    if strcmp(Parameters(iD).DiscretizationType,'CG')
      NodalSolution='SolutionGlobal';
      NumNodalComp=Sizes(iD).NumGlobalComp;
    else
      NodalSolution='SolutionLocal';
      NumNodalComp=Sizes(iD).NumLocalComp;
      Elements(iD).SolutionGlobal=mat2cell(...
        Block(iD,iD).SolutionGlobal(cell2mat(Faces(iD,iD).GlobalLocal),:,:),...
        ones(Sizes(iD).NumElements,1)*Sizes(iD).NumElementFaces*Sizes(iD).NumFaceNodes...
        *Sizes(iD).NumGlobalComp,size(Block(iD,iD).SolutionGlobal,3));
    end
    Elements(iD).(NodalSolution)=mat2cell(Block(iD,iD).(NodalSolution)(ElementNodes,:,:),...
      NodeBlocks,NumNodalComp,size(Block(iD,iD).(NodalSolution),3));
    if strcmp(Time.TimeDependent,'yes')
      Elements(iD).SolutionOld=mat2cell(Block(iD,iD).SolutionOld(ElementNodes,:,:),...
        NodeBlocks,NumNodalComp,size(Block(iD,iD).SolutionOld,3));
    else
      Elements(iD).SolutionOld=cell(Sizes(iD).NumElements,1);
    end
  end

  % Solve the HDG post-processing problem of each element
  for iD=1:nD
    if strcmp(Parameters(iD).DiscretizationType,'HDG') && ...
       strcmp(Parameters(iD).PostProcessingHDG,'yes')
      Elements(iD).LhsPost=[]; Elements(iD).RhsPost=[];
      Elements(iD)=Formulation{iD}.doPostProcess(Elements(iD),Simulation,Parameters(iD),...
        Faces(iD,iD),Time,RefElement(iD,iD),Sizes(iD));
      SolutionPost=MP*zeros(Sizes(iD).NumElementNodesPost*Sizes(iD).NumPostComp,...
                            Sizes(iD).NumElements);
      LhsPost=Elements(iD).LhsPost; RhsPost=Elements(iD).RhsPost;
      parfor iElem=1:Sizes(iD).NumElements
        SolutionPostElem=reshape(LhsPost{iElem}\RhsPost{iElem},...
          Sizes(iD).NumElementNodesPost,Sizes(iD).NumPostComp)'; %#ok
        SolutionPost(:,iElem)=SolutionPostElem(:);
      end
      Block(iD,iD).SolutionPost=SolutionPost;
      clear SolutionPost;
      Block(iD,iD).SolutionPost=reshape(Block(iD,iD).SolutionPost,Sizes(iD).NumPostComp,[])';
    end
  end

  % Store the final results
  if matchField(Options,'StoreTimeSteps')
    iST=find(Time.TimeStep==eval(Options.StoreTimeSteps));
  else
    iST=1;
  end
  for iD=1:nD
    Results=Formulation{iD}.storeResults(iD,iST,Results,Block,Simulation,Parameters,Mesh,Time,...
      RefElement,Sizes);
  end

  % Compute the Fast Fourier Transform
  if matchField(Options,'ComputeFFT')
    for iD=1:nD
      for iF=1:size(Options.ComputeFFT,1)
        if not(isempty(Results(iD).(Options.ComputeFFT{iF,1})))
          [Results(iD).Frequency,Results(iD).([Options.ComputeFFT{iF,1},'FFT'])]=...
            computeFFT(Results(iD).Time,Results(iD).(Options.ComputeFFT{iF,1}));
        end
      end
    end
  end

  % Report the processing time and memory
  Timer.Processing=toc(TimerProcessing);
  Memory.Processing=whos; Memory.Processing=sum(vertcat(Memory.Processing.bytes));
  fprintf('\nProcessing      completed in %.1f sec (%.1f MB)',...
    Timer.Processing,Memory.Processing/1e6);

  %% Post-processing

  TimerPostProcessing=tic;
  fprintf('\nPost-processing started...');

  % Compute the quantity of interest before the error computation
  if matchField(Options,'ComputeQuantityError')
    eval(Options.ComputeQuantityError);
  end

  % Plot the final solution
  if matchField(Options,'PlotSolution')
    if matchField(Results,'Time')
      iPlotT=numel(Results(1).Time);
    else
      iPlotT=size(System.Solution,3);
    end
    plotSolution(Simulation,Mesh,Faces,Results,Parameters,RefElement,Sizes,Options,iPlotT);
  end

  % Compute the error
  if matchField(Options,'ComputeError')
    if size(Options.ComputeError,2)==1
      Options.ComputeError(:,2)={'L2'};
    end
    if size(Options.ComputeError,2)==2
      Options.ComputeError(:,3)={'@(x,y,z) true'};
    end
    for iD=1:nD
      for iE=1:size(Options.ComputeError,1)
        Variable=Options.ComputeError{iE,1};
        if not(isempty(Results(iD).(Variable)))
          if not(strcmp(Variable(end-3:end),'Post'))
            MeshError=Mesh(iD); RefElementError=RefElement(iD,iD);
            VariableExact=Variable;
          else
            MeshError=Mesh(iD).Post; RefElementError=RefElement(iD,iD).Post;
            VariableExact=Variable(1:end-4);
          end
          Results(iD).(sprintf('%sError%s',Options.ComputeError{iE,1:2}))(iS1,iS2)=...
            sqrt(sum(computeError(Time,MeshError,Parameters(iD),RefElementError,Sizes(iD),...
              Results(iD).(Variable)(:,:,end),Parameters(iD).(VariableExact),...
              Options.ComputeError{iE,2},Options.ComputeError{iE,3},Simulation.Domain).^2));
        end
      end
    end
  end

  % Compute the quantity of interest after the error computation
  if matchField(Options,'ComputeQuantityPostError')
    eval(Options.ComputeQuantityPostError);
  end

  % Compute the convergence rate in space or in time
  if any(strcmp(Simulation.Type,{'ConvergenceSpace','ConvergenceTime'}))
    if strcmp(Simulation.Type,'ConvergenceSpace')
      Abscissa='ElementSize'; AbscissaName='Element size';
    else
      Abscissa='TimeStepSize'; AbscissaName='Time step size';
    end
    for iD=1:nD
      if strcmp(Simulation.Type,'ConvergenceSpace')
        Results(iD).Degree(iS1,1)=Parameters(iD).Degree;
        if isa(Options.CharacteristicElementSize,'char')
          Results(iD).ElementSize(iS1,iS2)=...
            Mesh(iD).([Options.CharacteristicElementSize,'ElementSize']);
        elseif isa(Options.CharacteristicElementSize,'numeric')
          Results(iD).ElementSize(iS1,iS2)=Options.CharacteristicElementSize(iS1,iS2);
        end
        if strcmp(Time.TimeDependent,'yes')
          Results(iD).TimeStepSize(iS1,1)=Time.TimeStepSize;
        end
      else
        Results(iD).BDFOrder(iS1,1)=Time.BDFOrder;
        Results(iD).TimeStepSize(iS1,iS2)=Time.TimeStepSize;
        if isa(Options.CharacteristicElementSize,'char')
          Results(iD).ElementSize(iS1,1)=...
            Mesh(iD).([Options.CharacteristicElementSize,'ElementSize']);
        elseif isa(Options.CharacteristicElementSize,'numeric')
          Results(iD).ElementSize(iS1,1)=Options.CharacteristicElementSize(iS1);
        end
      end
      for iE=1:size(Options.ComputeError,1)
        if not(isempty(Results(iD).(Options.ComputeError{iE,1})))
          Error=sprintf('%sError%s',Options.ComputeError{iE,1:2});
          if iS2==1
            Results(iD).([Error,'Rate'])(iS1,iS2)=NaN;
          else
            Results(iD).([Error,'Rate'])(iS1,iS2)=...
              log(Results(iD).(Error)(iS1,iS2)/Results(iD).(Error)(iS1,iS2-1))...
             /log(Results(iD).(Abscissa)(iS1,iS2)/Results(iD).(Abscissa)(iS1,iS2-1));
          end
        end
      end
    end

    % Plot the convergence
    if strcmp(Options.PlotConvergence,'yes')
      if iS==1
        figure('Color','w')
        NumPlots=size(Options.ComputeError,1);
        if NumPlots>1
          ScreenSize=get(0,'ScreenSize');
          set(gcf,'Position',[(3-ceil(NumPlots/4))/3*ScreenSize(3),0,...
                                 ceil(NumPlots/4) /3*ScreenSize(3),ScreenSize(4)])
        end
      elseif iS2>1
        for iE=1:NumPlots
          for iD=1:nD
            subplot(ceil(NumPlots/ceil(NumPlots/4)),ceil(NumPlots/4),iE)
            if not(isempty(Results(iD).(Options.ComputeError{iE,1})))
              Error=sprintf('%sError%s',Options.ComputeError{iE,1:2});
              plot( log2(Results(iD).(Abscissa)(iS1,iS2-1:iS2)),...
                   log10(Results(iD).(Error)(iS1,iS2-1:iS2)),...
                   [Colors{iS1},Styles{iD},Markers{iS1}]);
              text( log2(Results(iD).(Abscissa)(iS1,iS2)),...
                   log10(Results(iD).(Error)(iS1,iS2)),...
                   sprintf('%.1f',Results(iD).([Error,'Rate'])(iS1,iS2)),...
                   'HorizontalAlignment','Center','VerticalAlignment','Top','Color',Colors{iS1})
            end
            hold on
          end
          grid on; box on;
          title(Options.ComputeError{iE,1});
          xlabel(sprintf('log_{2}(%s)',AbscissaName));
          ylabel(['log_{10}(Error-',Options.ComputeError{iE,2},')']);
          pause(eps);
        end
      end
    end
  end

  % Store the parameter of the parametric study
  if strcmp(Simulation.Type,'ParametricStudy')
    for iD=1:nD
      Results(iD).(Simulation.ParameterStudy)(iS2)=Parameters(iD).(Simulation.ParameterStudy);
    end
  end

  % Compute the parallel performance
  if any(strcmp(Simulation.Type,{'ScalingStrong','ScalingWeak'}))
    Weak=strcmp(Simulation.Type,'ScalingWeak');
    if Weak
      Scaling='Weak'; Metric='Efficiency';
    else
      Scaling='Strong'; Metric='Speedup';
    end
    for iD=1:nD
      Results(iD).NumProcessors(iS2)=Simulation.NumProcessors;
      if Weak
        Results(iD).NumElements(iS2)=sum([Sizes.NumElements]);
      end
      Results(iD).Timer(iS2)=sum(Timer.(Options.CharacteristicTimer));
      Results(iD).(Metric)(iS2)=Results(iD).Timer(1)/Results(iD).Timer(iS2);
    end

    % Plot the parallel performance
    if strcmp(Options.(['Plot',Metric]),'yes')
      if iS==1
        figure('Color','w')
      elseif iS2>1
        if Weak
          plot(Results(1).NumProcessors(iS2-1:iS2),[1,1],'k--')
        else
          plot(Results(1).NumProcessors(iS2-1:iS2),Results(1).NumProcessors(iS2-1:iS2),'k--')
        end
        hold on
        plot(Results(1).NumProcessors(iS2-1:iS2),Results(1).(Metric)(1,iS2-1:iS2),'b-o');
        text(Results(1).NumProcessors(iS2),Results(1).(Metric)(iS2),...
             sprintf('%.1f',Results(1).(Metric)(iS2)),...
             'HorizontalAlignment','Center','VerticalAlignment','Top','Color','b')
        grid on; box on;
        title(sprintf('%s scaling',Scaling));
        xlabel('Number of processors');
        ylabel(sprintf('%s (Timer.%s)',Metric,Options.CharacteristicTimer));
        if Weak
          ylim([0,1.2])
        end
        pause(eps);
      end
    end
  end

  % Compute the quantity of interest at the end of the simulation
  if matchField(Options,'ComputeQuantityEnd')
    eval(Options.ComputeQuantityEnd);
  end

  % Export the final solution to Paraview
  if matchField(Options,'Export2Paraview')
    for iD=1:nD
      FileNameExport=[FileName,repmat(sprintf('_discr_%d',iD),1,nD>1)];
      export2Paraview(Mesh(iD),Results(iD),Parameters(iD),Time,Sizes(iD),Options,FileNameExport);
      if strcmp(Parameters(iD).DiscretizationType,'HDG') && ...
         strcmp(Parameters(iD).PostProcessingHDG,'yes')
        export2Paraview(Mesh(iD).Post,Results(iD),Parameters(iD),Time,Sizes(iD),Options,...
          [FileNameExport,'_post']);
      end
    end
  end

  % Report the post-processing time and memory
  Timer.PostProcessing=toc(TimerPostProcessing);
  Memory.PostProcessing=whos; Memory.PostProcessing=sum(vertcat(Memory.PostProcessing.bytes));
  fprintf('\nPost-processing completed in %.1f sec (%.1f MB)',...
    Timer.PostProcessing,Memory.PostProcessing/1e6);

  % Print the errors and the convergence rates
  if matchField(Options,'ComputeError')
    Converging=any(strcmp(Simulation.Type,{'ConvergenceSpace','ConvergenceTime'}));
    fprintf('\n\nVariable%21sError%6sNorm','','');
    fprintf('%s',repmat(sprintf('%3sRate',''),1,Converging));
    for iD=1:nD
      for iE=1:size(Options.ComputeError,1)
        if not(isempty(Results(iD).(Options.ComputeError{iE,1})))
          Error=sprintf('%sError%s',Options.ComputeError{iE,1:2});
          fprintf('\n%s',Options.ComputeError{iE,1});
          fprintf('%s',repmat(sprintf('_%d',iD),1,nD>1));
          fprintf('%*s',27-numel(Options.ComputeError{iE,1}),'');
          fprintf('%.2e',Results(iD).(Error)(iS1,iS2));
          fprintf('%3s%s','',Options.ComputeError{iE,2});
          if Converging
            fprintf('%*s%.2f',7-numel(Options.ComputeError{iE,2}),'',...
              Results(iD).([Error,'Rate'])(iS1,iS2));
          end
        end
      end
    end
  end

  % Print the info of the simulation
  switch Simulation.Type
    case 'ConvergenceSpace'
      fprintf('\n\nSimulation with %s completed',strjoin(arrayfun(@(iD) sprintf(...
        'k = %d and h = %.2e',Parameters(iD).Degree,Results(iD).ElementSize(iS1,iS2)),...
        1:nD,'UniformOutput',false),'\n       and with '));
    case 'ConvergenceTime'
      fprintf('\n\nSimulation with BDFo = %d and dt = %.2e completed',...
        Time.BDFOrder,Results(1).TimeStepSize(iS1,iS2));
    case 'ParametricStudy'
      fprintf('\n\nSimulation with %s completed',strjoin(arrayfun(@(iD) sprintf(...
        '%s = %.2e',Simulation.ParameterStudy,Results(iD).(Simulation.ParameterStudy)(iS2)),...
        1:nD,'UniformOutput',false),'\n       and with '));
    case 'ScalingStrong'
      fprintf('\n\nSimulation with np = %d completed',Results(1).NumProcessors(iS2));
    case 'ScalingWeak'
      fprintf('\n\nSimulation with np = %d and nel = %d completed',...
        Results(1).NumProcessors(iS2),Results(1).NumElements(iS2));
  end

  fprintf('\n--------------------------------------------------------------------------------\n');

  % Clear the quantities that the next simulation rebuilds
  if iS<Simulation.NumSimulations
    clear RefElement Faces Block Mesh;
  end
end

% Perform the test
if matchField(Options,'Test')
  Simulation.TestPassed=eval(Options.Test);
  if Simulation.TestPassed; fprintf('\nTEST PASSED\n');
  else;                     fprintf('\nTEST FAILED\n'); end
end

% Hold off the potential figures
if not(isempty(findobj('type','figure')))
  hold off
end

% Clear the temporary variables of the workflow
if not(matchField(Simulation,'ComputeOnlyPreProcessing','yes'))
  clear('-regexp','^i[A-Z0-9]','^Plot','^Timer[A-Z]',...
    '^(Abscissa|AbscissaName|BoundaryNames|CP|Children|ColorMap|Colors|Converging)$',...
    '^(Eigenvalues|Eigenvectors|ElementNodes|EmptyMesh|Error|FileNameExport|HDG|HDGPost)$',...
    '^(HasLocal|LhsPost|Limits|Markers|MeshError|MeshPost|Metric|NodalSolution|NodeBlocks)$',...
    '^(NodesBC|NumColors|NumMeshNodes|NumMeshes|NumNodalComp|NumPlots|NumVertexNodes)$',...
    '^(Offset|Order|ParallelSettings|Pool|RefElementError|RhsPost|Scaling|ScreenSize|Shift)$',...
    '^(SolutionPost|SolutionPostElem|SolutionPred|Styles|Variable|VariableExact|Weak|Wide)$',...
    '^(kP|nD)$');
end

% Save the results
if matchField(Options,'SaveResults','yes')
  save(sprintf('%s/output/%s%s.mat',fileparts(which(mfilename)),Options.SaveResultsFolder,...
    FileName),'-regexp','^(?!(CPUTime|Passed|Test)$).',Options.SaveResultsMatVersion);
  fprintf('\nSaved output file in: output/%s%s.mat',Options.SaveResultsFolder,FileName);
  fprintf('\nSaved  diary file in: output/%s%s.txt\n',Options.SaveResultsFolder,FileName);
  diary(sprintf('%s/output/%s%s.txt',fileparts(which(mfilename)),Options.SaveResultsFolder,...
    FileName)); diary off;
end
