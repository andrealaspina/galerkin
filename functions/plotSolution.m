function plotSolution(Simulation,Mesh,Faces,Results,Parameters,RefElement,Sizes,Options,Index)
% Plot the real (and imaginary) part of each component of each variable of each discretization

% Axis limits
Limits=minmax(double([Mesh.Nodes]));
Limits=Limits+[-1,+1].*(Limits(:,1)==Limits(:,2));
ThreeD=Sizes(1).NumSpaceDim==3;

for iP=1:length(Options.PlotSolution)
  Variable=Options.PlotSolution{iP};
  Post=strcmp(Variable(end-3:end),'Post');
  NumComp=size(vertcat(Results.(Variable)),2);

  % Abscissa against which the variable is stored
  if     matchField(Results,'Time') && not(matchField(Results,'Frequency'))
    Domain='Time';
  elseif matchField(Results,'Frequency') && not(matchField(Results,'Time'))
    Domain='Frequency';
  elseif find(strcmp(fieldnames(Results),'Time'))<find(strcmp(fieldnames(Results),'Frequency'))
    if not(strcmp(Variable(end-2:end),'FFT'))
      Domain='Time';
    else
      Domain='Frequency';
    end
  elseif find(strcmp(fieldnames(Results),'Frequency'))<find(strcmp(fieldnames(Results),'Time'))
    if not(strcmp(Variable(end-2:end),'IFFT'))
      Domain='Frequency';
    else
      Domain='Time';
    end
  elseif matchField(Results,'Mode')
    Domain='Mode';
  end

  for iC=1:NumComp
    figure('Color','w');
    isComplex=not(isreal(Results(1).(Variable)));
    for iRI=1:(1+isComplex)
      subplot(1,1+isComplex,iRI)
      for iD=1:Simulation.NumDiscretizations
        if not(isempty(Results(iD).(Variable)))
          hold on
          if Post
            MeshPlot=Mesh(iD).Post;
            RefElementPlot=RefElement(iD,iD).Post;
          else
            MeshPlot=Mesh(iD);
            RefElementPlot=RefElement(iD,iD);
          end
          u=Results(iD).(Variable)(:,iC,Index);
          if iRI==1
            u=real(u);
          elseif iRI==2
            u=imag(u);
          end
          if Sizes(iD).NumSpaceDim==2

            % Surface of each element
            for iElem=1:Sizes(iD).NumElements
              Ce=MeshPlot.Elements(:,iElem)';
              Xe=double(MeshPlot.Nodes(1:2,Ce)');
              if Parameters(iD).Degree==0
                ue=repelem(double(u(iElem)),Sizes(iD).NumSpaceDim+1);
              else
                ue=double(u(Ce));
              end
              trisurf(delaunayn(Xe),Xe(:,1),Xe(:,2),ue)
            end
          elseif strcmp(Parameters(iD).DiscretizationType,'HDG')

            % Faces of the sub-elements of each element on the exterior boundary
            for iFaceExt=1:size(Faces(iD,iD).Exterior,1)
              iElem=Faces(iD,iD).Exterior(iFaceExt,1);
              Ce=MeshPlot.Elements(:,iElem)';
              Xe=double(MeshPlot.Nodes(:,Ce)');
              if Parameters(iD).Degree==0
                ue=repelem(double(u(iElem)),Sizes(iD).NumSpaceDim+1);
              else
                ue=double(u(Ce));
              end
              De=delaunayn(Xe);
              for iElemSub=1:size(De,1)
                Xes=Xe(De(iElemSub,:),:);
                ues=ue(De(iElemSub,:),:);
                for iFaceSub=1:Sizes(iD).NumElementFaces
                  Dfs=RefElementPlot.FaceNodesElem(iFaceSub,1:3);
                  trisurf(1:3,Xes(Dfs,1),Xes(Dfs,2),Xes(Dfs,3),ues(Dfs,:))
                end
              end
            end
            view(3)
          else

            % Surface of the whole mesh
            pdeplot3D(Mesh(iD).Nodes,Mesh(iD).Elements,'ColorMapData',u)
            Children=get(gca,'Children'); delete(Children([2,3]));
          end
        end
        hold on
      end
      hold off; box on; axis equal; shading interp; colormap jet; colorbar;
      set(gca,'Visible','on')
      title([repmat('real(',1,isComplex&&iRI==1),...
             repmat('imag(',1,isComplex&&iRI==2),...
             Variable,...
             repmat(sprintf(' (%d)',iC),1,NumComp>1),...
             repmat(')',1,isComplex),...
             sprintf(' at %s %.2e',lower(Domain),Results(1).(Domain)(Index))]);
      xlim(Limits(1,:)); xlabel('x');
      ylim(Limits(2,:)); ylabel('y');
      if ThreeD
        zlim(Limits(3,:)); zlabel('z');
      end
    end
    pause(eps);
  end
end

end
