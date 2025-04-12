function []=plotSolution(...
         Simulation,Mesh,Faces,Results,Parameters,RefElement,Sizes,Options,Index)
         % Plot solution

% Get axis limits
MinMax=minmax(double([Mesh.Nodes]));
XMinMax=MinMax(1,:)+[-1,+1]*(MinMax(1,1)==MinMax(1,2));
YMinMax=MinMax(2,:)+[-1,+1]*(MinMax(2,1)==MinMax(2,2));
if Sizes(1).NumSpaceDim==3
  ZMinMax=MinMax(3,:)+[-1,+1]*(MinMax(3,1)==MinMax(3,2));
end

% Plot real (and imaginary) part of each component of each variable of each discretization
for iP=1:length(Options.PlotSolution)
  if     matchField(Results,'Time') && not(matchField(Results,'Frequency'))
    Domain='Time';
  elseif matchField(Results,'Frequency') && not(matchField(Results,'Time'))
    Domain='Frequency';
  elseif find(strcmp(fieldnames(Results),'Time')) < find(strcmp(fieldnames(Results),'Frequency'))
    if not(strcmp(Options.PlotSolution{iP}(end-2:end),'FFT'))
      Domain='Time';
    else
      Domain='Frequency';
    end
  elseif find(strcmp(fieldnames(Results),'Frequency')) < find(strcmp(fieldnames(Results),'Time'))
    if not(strcmp(Options.PlotSolution{iP}(end-2:end),'IFFT'))
      Domain='Frequency';
    else
      Domain='Time';
    end
  elseif matchField(Results,'Mode')
    Domain='Mode';
  end
  for iC=1:size(vertcat(Results.([Options.PlotSolution{iP}])),2)
    figure('Color','w');
    isComplex=not(isreal(Results(1).([Options.PlotSolution{iP}])));
    for iRI=1:(1+isComplex)
      subplot(1,1+isComplex,iRI)
      for iD=1:Simulation.NumDiscretizations
        if not(isempty(Results(iD).([Options.PlotSolution{iP}])))
          hold on
          if Sizes(iD).NumSpaceDim==2
            for iElem=1:Sizes(iD).NumElements
              if not(strcmp(Options.PlotSolution{iP}(end-3:end),'Post'))
                Ce=Mesh(iD).Elements(:,iElem)';
                Xe=double(Mesh(iD).Nodes(1:2,Ce)');
              else
                Ce=Mesh(iD).Post.Elements(:,iElem)';
                Xe=double(Mesh(iD).Post.Nodes(1:2,Ce)');
              end
              if Parameters(iD).Degree==0
                ue=repelem(double(Results(iD).([Options.PlotSolution{iP}])(iElem,iC,Index)),...
                  Sizes(iD).NumSpaceDim+1);
              else
                ue=double(Results(iD).([Options.PlotSolution{iP}])(Ce,iC,Index));
              end
              if iRI==1
                ue=real(ue);
              elseif iRI==2
                ue=imag(ue);
              end
              De=delaunayn(Xe);
              trisurf(De,Xe(:,1),Xe(:,2),ue)
            end
          elseif Sizes(iD).NumSpaceDim==3
            if strcmp(Parameters(iD).DiscretizationType,'HDG')
              for iFaceExt=1:size(Faces(iD,iD).Exterior,1)
                iElem=Faces(iD,iD).Exterior(iFaceExt,1);
                if not(strcmp(Options.PlotSolution{iP}(end-3:end),'Post'))
                  Ce=Mesh(iD).Elements(:,iElem)';
                  Xe=double(Mesh(iD).Nodes(:,Ce)');
                else
                  Ce=Mesh(iD).Post.Elements(:,iElem)';
                  Xe=double(Mesh(iD).Post.Nodes(:,Ce)');
                end
                if Parameters(iD).Degree==0
                  ue=repelem(double(Results(iD).([Options.PlotSolution{iP}])(iElem,iC,Index)),...
                    Sizes(iD).NumSpaceDim+1);
                else
                  ue=double(Results(iD).([Options.PlotSolution{iP}])(Ce,iC,Index));
                end
                if iRI==1
                  ue=real(ue);
                elseif iRI==2
                  ue=imag(ue);
                end
                De=delaunayn(Xe);
                for iElemSub=1:size(De,1)
                  Ces=De(iElemSub,:);
                  Xes=Xe(Ces,:);
                  ues=ue(Ces,:);
                  for iFaceSub=1:Sizes(iD).NumElementFaces
                    if not(strcmp(Options.PlotSolution{iP}(end-3:end),'Post'))
                      Dfs=RefElement(iD,iD).FaceNodesElem(iFaceSub,1:3);
                    else
                      Dfs=RefElement(iD,iD).Post.FaceNodesElem(iFaceSub,1:3);
                    end
                    Xfs=Xes(Dfs,:);
                    ufs=ues(Dfs,:);
                    trisurf(1:3,Xfs(:,1),Xfs(:,2),Xfs(:,3),ufs)
                  end
                end
              end
              view(3)
            elseif strcmp(Parameters(iD).DiscretizationType,'CG')
              u=Results(iD).([Options.PlotSolution{iP}])(:,iC,Index);
              if iRI==1
                u=real(u);
              elseif iRI==2
                u=imag(u);
              end
              pdeplot3D(Mesh(iD).Nodes,Mesh(iD).Elements,'ColorMapData',u)
              ChildrenAux=get(gca,'Children'); delete(ChildrenAux([2,3]));
            end
          end
        end
        hold on
      end
      hold off
      box on
      axis equal
      shading interp
      colormap jet
      colorbar
      set(gca,'Visible','on')
      title([repmat('real(',isComplex&(iRI==1)),...
             repmat('imag(',isComplex&(iRI==2)),...
             Options.PlotSolution{iP},...
             repmat(sprintf(' (%d)',iC),size(vertcat(Results.([Options.PlotSolution{iP}])),2)>1),...
             repmat(')',isComplex)...              
             sprintf(' at %s %.2e',lower(Domain),Results(1).(Domain)(Index))]);
      xlim(XMinMax); xlabel('x'); ylim(YMinMax); ylabel('y');
      if Sizes(iD).NumSpaceDim==3
        zlim(ZMinMax); zlabel('z');
      end
    end
    pause(eps);
  end
end

end