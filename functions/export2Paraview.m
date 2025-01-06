function []=export2Paraview(...
         Formulation,Mesh,Results,Parameters,Sizes,Options,FileName)
         % Export solution to Paraview

% Get general parameters
nsd=Sizes.NumSpaceDim;
k=Parameters.Degree;
C=Mesh.Elements';
X=Mesh.Nodes';

% Extrapolate mesh info
NumElements=size(C,1);
NumElementNodes=size(C,2);
NumNodes=size(X,1);

% Understand if it is a postprocessed solution on a k+1 mesh
if     (nsd==2 && NumElementNodes==(k+1)*(k+2)/2)       || ...
       (nsd==3 && NumElementNodes==(k+1)*(k+2)*(k+3)/6)
  isPostProcess=false;
elseif (nsd==2 && NumElementNodes==((k+1)+1)*((k+1)+2)/2)           || ...
       (nsd==3 && NumElementNodes==((k+1)+1)*((k+1)+2)*((k+1)+3)/6)
  isPostProcess=true;
  k=k+1;
end

% Cell type (Lagrange triangle/tetrahedron)
if nsd==2
  CellType=69*ones(NumElements,1);
elseif nsd==3
  CellType=71*ones(NumElements,1);
end
if strcmp(Parameters.NodesDistribution,'Fekete')
  warning('Paraview assumes uniform nodes distribution only!');
end

% Nodes re-ordering
if nsd==2
  switch k
    case 1
      order=[1,2,3];
    case 2
      order=[1,2,3,4,5,6];
    case 3
      order=[1,2,3,4,5,6,7,8,9,10];
    case 4
      order=[1,2,3,4,5,6,7,8,9,10,11,12,15,13,14];
    case 5
      order=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,16,18,21,17,19];
    case 6
      order=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25,19,22,26,27,20,21,23,24,28];
    case 7
      order=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,30,22,26,31,32,33,23,24,25,27,28,29,36,34,35];
    case 8
      order=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,35,25,30,36,37,38,39,26,27,28,29,31,32,33,34,44,40,42,45,41,43];
    otherwise
      error('Element not yet implemented');
  end
elseif nsd==3
  switch k
    case 1
      order=[1,2,3,4];
    case 2
      order=[1,2,3,4,5,6,7,8,9,10];
    case 3
      order=[1,2,3,4,5,6,7,8,10,9,11,12,13,14,15,16,19,18,20,17];
    case 4
      order=[1,2,3,4,5,6,7,8,9,10,13,12,11,14,15,16,17,18,19,20,21,22,29,30,31,28,26,27,32,33,34,25,23,24,35];
    case 5
      order=[1,2,3,4,5,6,7,8,9,10,11,12,16,15,14,13,17,18,19,20,21,22,23,24,25,26,27,28,41,43,45,42,44,46,39,35,37,40,36,38,47,49,51,48,50,52,33,29,31,34,30,32,53,54,55,56];
    case 6
      order=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,19,18,17,16,15,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,55,58,61,56,57,59,60,62,63,64,51,45,48,52,53,46,47,49,50,54,65,68,71,66,67,69,70,72,73,74,41,35,38,42,43,36,37,39,40,44,75,76,77,78,79,80,81,82,83,84];
    case 7
      order=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,22,21,20,19,18,17,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,71,75,79,72,73,74,76,77,78,80,81,82,83,84,85,64,56,60,65,66,67,57,58,59,61,62,63,70,68,69,86,90,94,87,88,89,91,92,93,95,96,97,98,99,100,49,41,45,50,51,52,42,43,44,46,47,48,55,53,54,101,102,103,104,105,106,107,108,110,109,111,112,113,114,115,116,119,118,120,117];
    case 8
      order=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25,24,23,22,21,20,19,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,89,94,99,90,91,92,93,95,96,97,98,100,101,102,103,104,106,108,105,107,109,78,68,73,79,80,81,82,69,70,71,72,74,75,76,77,87,83,85,88,84,86,110,115,120,111,112,113,114,116,117,118,119,121,122,123,124,125,127,129,126,128,130,57,47,52,58,59,60,61,48,49,50,51,53,54,55,56,66,62,64,67,63,65,131,132,133,134,135,136,137,138,139,140,143,142,141,144,145,146,147,148,149,150,151,152,159,160,161,158,156,157,162,163,164,155,153,154,165];
    otherwise
      error('Element not yet implemented');
  end
end

% Connectivity for paraview
ConnecVTK=[NumElementNodes*ones(NumElements,1),C(:,order)-1];

% Open file
vtkFile=[fileparts(which(mfilename)),'/../output/',Options.SaveResultsFolder,FileName,'.vtk'];
vtkId=fopen(vtkFile,'w');

% Write header
fprintf(vtkId,'# vtk DataFile Version 3.0\n');
fprintf(vtkId,'My problem\n');
fprintf(vtkId,'\nASCII\n');

% Unstructured grid
fprintf(vtkId,'\nDATASET UNSTRUCTURED_GRID\n');

% Write coordinates
fprintf(vtkId,'\nPOINTS %d float\n',NumNodes);
fprintf(vtkId,'%.12f %.12f %.12f\n',X');

% Write connectivity
fprintf(vtkId,'\nCELLS %d %d\n',NumElements,NumElements*(NumElementNodes+1));
fprintf(vtkId,sprintf('%s',['%d',repmat(' %d',1,NumElementNodes),'\n']),ConnecVTK');

% Write cell types
fprintf(vtkId,'\nCELL_TYPES %d\n',NumElements);
fprintf(vtkId,'%d\n',CellType);

% Get data
[PointData,CellData]=Formulation.dataForParaview(Results,Parameters,Mesh,Sizes,isPostProcess);

% Write point data
fprintf(vtkId,'\nPOINT_DATA %d\n',NumNodes);
fprintf(vtkId,'%c',PointData);

% Write cell data
fprintf(vtkId,'\nCELL_DATA %d\n',NumElements);
fprintf(vtkId,'%c',CellData);

% Close file
fclose(vtkId);

% Open paraview
if matchField(Options,'OpenParaview','yes')
  system(['/Applications/ParaView-5.10.0-RC1.app/Contents/MacOS/paraview',...
          ' --data=','''',vtkFile,'''',' &']);
end

end