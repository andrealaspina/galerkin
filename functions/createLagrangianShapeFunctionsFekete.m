% Initialization
close all; clear; clc;

% Input
k=3
nsd=1

% Giacomini's element
if nsd==1
  Ge=load('/Users/andrealaspina/Desktop/HDGlab/dat/refElem/refFaceTRI_opt0.mat'); Ge=Ge.refFace(k,k);
elseif nsd==2
  Ge=load('/Users/andrealaspina/Desktop/HDGlab/dat/refElem/refElemTRI_opt0.mat'); Ge=Ge.refElem(k);
  Gf=load('/Users/andrealaspina/Desktop/HDGlab/dat/refElem/refFaceTRI_opt0.mat'); Gf=Gf.refFace(k,k);
elseif nsd==3
  Ge=load('/Users/andrealaspina/Desktop/HDGlab/dat/refElem/refElemTET_opt0.mat'); Ge=Ge.refElem(k);
  Gf=load('/Users/andrealaspina/Desktop/HDGlab/dat/refElem/refFaceTET_opt0.mat'); Gf=Gf.refFace(k,k);
end

% Fekete's element
if nsd==1
  Fe=load('/Users/andrealaspina/Desktop/HDGlab/dat/refElem/refFaceTRI_opt1.mat'); Fe=Fe.refFace(k,k);
elseif nsd==2
  Fe=load('/Users/andrealaspina/Desktop/HDGlab/dat/refElem/refElemTRI_opt1.mat'); Fe=Fe.refElem(k);
  Ff=load('/Users/andrealaspina/Desktop/HDGlab/dat/refElem/refFaceTRI_opt1.mat'); Ff=Ff.refFace(k,k);
elseif nsd==3
  Fe=load('/Users/andrealaspina/Desktop/HDGlab/dat/refElem/refElemTET_opt1.mat'); Fe=Fe.refElem(k);
  Ff=load('/Users/andrealaspina/Desktop/HDGlab/dat/refElem/refFaceTET_opt1.mat'); Ff=Ff.refFace(k,k);
end

% Andrea's element
A=createReferenceElement(max(nsd,2),k,k,'Uniform');

% Plot
figure('color','w')
subplot(2,3,1)
if nsd==1
  plot(Ge.coordinates([1,end]),[0,0],'k'); hold on
  plot(Ge.coordinates,0,'ko','MarkerSize',10)
  arrayfun(@(i) text(Ge.coordinates(i),0,sprintf('%d',i),'HorizontalAlignment','Center'),1:length(Ge.coordinates))
elseif nsd==2
  plot(Ge.coordinates([1,2,3,1],1),Ge.coordinates([1,2,3,1],2),'k'); hold on
  plot(Ge.coordinates(:,1),Ge.coordinates(:,2),'ko','MarkerSize',10)
  arrayfun(@(i) text(Ge.coordinates(i,1),Ge.coordinates(i,2),sprintf('%d',i),'HorizontalAlignment','Center'),1:size(Ge.coordinates,1))
elseif nsd==3
  plot3(Ge.coordinates([1,2,3,1,4,3,4,2],1),Ge.coordinates([1,2,3,1,4,3,4,2],2),Ge.coordinates([1,2,3,1,4,3,4,2],3),'k'); hold on
  plot3(Ge.coordinates(:,1),Ge.coordinates(:,2),Ge.coordinates(:,3),'ko','MarkerSize',10)
  arrayfun(@(i) text(Ge.coordinates(i,1),Ge.coordinates(i,2),Ge.coordinates(i,3),sprintf('%d',i),'HorizontalAlignment','Center'),1:size(Ge.coordinates,1))
end
axis equal; xlim([-1.2,+1.2]); ylim([-1.2,+1.2]); zlim([-1.2,+1.2]);
title('Giacomini''s element');
subplot(2,3,4)
if nsd==2
  plot(Gf.coordinates([1,end]),[0,0],'k'); hold on
  plot(Gf.coordinates,0,'ko','MarkerSize',10)
  arrayfun(@(i) text(Gf.coordinates(i),0,sprintf('%d',i),'HorizontalAlignment','Center'),1:length(Gf.coordinates))
elseif nsd==3
  plot(Gf.coordinates([1,2,3,1],1),Gf.coordinates([1,2,3,1],2),'k'); hold on
  plot(Gf.coordinates(:,1),Gf.coordinates(:,2),'ko','MarkerSize',10)
  arrayfun(@(i) text(Gf.coordinates(i,1),Gf.coordinates(i,2),sprintf('%d',i),'HorizontalAlignment','Center'),1:size(Gf.coordinates,1))
end
axis equal; xlim([-1.2,+1.2]); ylim([-1.2,+1.2]); zlim([-1.2,+1.2]);
title('Giacomini''s face');
subplot(2,3,2)
if nsd==1
  plot(Fe.coordinates([1,end]),[0,0],'k'); hold on
  plot(Fe.coordinates,0,'ko','MarkerSize',10)
  arrayfun(@(i) text(Fe.coordinates(i),0,sprintf('%d',i),'HorizontalAlignment','Center'),1:length(Fe.coordinates))
elseif nsd==2
  plot(Fe.coordinates([1,2,3,1],1),Fe.coordinates([1,2,3,1],2),'k'); hold on
  plot(Fe.coordinates(:,1),Fe.coordinates(:,2),'ko','MarkerSize',10)
  arrayfun(@(i) text(Fe.coordinates(i,1),Fe.coordinates(i,2),sprintf('%d',i),'HorizontalAlignment','Center'),1:size(Fe.coordinates,1))
elseif nsd==3
  plot3(Fe.coordinates([1,2,3,1,4,3,4,2],1),Fe.coordinates([1,2,3,1,4,3,4,2],2),Fe.coordinates([1,2,3,1,4,3,4,2],3),'k'); hold on
  plot3(Fe.coordinates(:,1),Fe.coordinates(:,2),Fe.coordinates(:,3),'ko','MarkerSize',10)
  arrayfun(@(i) text(Fe.coordinates(i,1),Fe.coordinates(i,2),Fe.coordinates(i,3),sprintf('%d',i),'HorizontalAlignment','Center'),1:size(Fe.coordinates,1))
end
axis equal; xlim([-1.2,+1.2]); ylim([-1.2,+1.2]); zlim([-1.2,+1.2]);
title('Fekete''s element');
subplot(2,3,5)
if nsd==2
  plot(Ff.coordinates([1,end]),[0,0],'k'); hold on
  plot(Ff.coordinates,0,'ko','MarkerSize',10)
  arrayfun(@(i) text(Ff.coordinates(i),0,sprintf('%d',i),'HorizontalAlignment','Center'),1:length(Ff.coordinates))
elseif nsd==3
  plot(Ff.coordinates([1,2,3,1],1),Ff.coordinates([1,2,3,1],2),'k'); hold on
  plot(Ff.coordinates(:,1),Ff.coordinates(:,2),'ko','MarkerSize',10)
  arrayfun(@(i) text(Ff.coordinates(i,1),Ff.coordinates(i,2),sprintf('%d',i),'HorizontalAlignment','Center'),1:size(Ff.coordinates,1))
end
axis equal; xlim([-1.2,+1.2]); ylim([-1.2,+1.2]); zlim([-1.2,+1.2]);
title('Fekete''s face');
subplot(2,3,3)
if nsd==1
  plot(A.NodesCoordFace([1,2]),[1/2,1/2],'k'); hold on
  plot(A.NodesCoordFace,1/2,'ko','MarkerSize',10)
  arrayfun(@(i) text(A.NodesCoordFace(i),1/2,sprintf('%d',i),'HorizontalAlignment','Center'),1:length(A.NodesCoordFace))
elseif nsd==2
  plot(A.NodesCoordElem([1,2,3,1],1),A.NodesCoordElem([1,2,3,1],2),'k'); hold on
  plot(A.NodesCoordElem(:,1),A.NodesCoordElem(:,2),'ko','MarkerSize',10)
  arrayfun(@(i) text(A.NodesCoordElem(i,1),A.NodesCoordElem(i,2),sprintf('%d',i),'HorizontalAlignment','Center'),1:size(A.NodesCoordElem,1))
elseif nsd==3
  plot3(A.NodesCoordElem([1,2,3,1,4,3,4,2],1),A.NodesCoordElem([1,2,3,1,4,3,4,2],2),A.NodesCoordElem([1,2,3,1,4,3,4,2],3),'k'); hold on
  plot3(A.NodesCoordElem(:,1),A.NodesCoordElem(:,2),A.NodesCoordElem(:,3),'ko','MarkerSize',10)
  arrayfun(@(i) text(A.NodesCoordElem(i,1),A.NodesCoordElem(i,2),A.NodesCoordElem(i,3),sprintf('%d',i),'HorizontalAlignment','Center'),1:size(A.NodesCoordElem,1))
end
axis equal; xlim([-0.2,+1.2]); ylim([-0.2,+1.2]); zlim([-0.2,+1.2]);
title('Andrea''s element');
subplot(2,3,6)
if nsd==2
  plot(A.NodesCoordFace([1,2]),[1/2,1/2],'k'); hold on
  plot(A.NodesCoordFace,1/2,'ko','MarkerSize',10)
  arrayfun(@(i) text(A.NodesCoordFace(i),1/2,sprintf('%d',i),'HorizontalAlignment','Center'),1:length(A.NodesCoordFace))
elseif nsd==3
  plot(A.NodesCoordFace([1,2,3,1],1),A.NodesCoordFace([1,2,3,1],2),'k'); hold on
  plot(A.NodesCoordFace(:,1),A.NodesCoordFace(:,2),'ko','MarkerSize',10)
  arrayfun(@(i) text(A.NodesCoordFace(i,1),A.NodesCoordFace(i,2),sprintf('%d',i),'HorizontalAlignment','Center'),1:size(A.NodesCoordFace,1))
end
axis equal; xlim([-0.2,+1.2]); ylim([-0.2,+1.2]); zlim([-0.2,+1.2]);
title('Andrea''s face');
pause(eps);

% Conversion Giacomini=Fekete-->Andrea
if nsd==1
  [~,Oe]=ismembertol(A.NodesCoordFace,(Ge.coordinates'+1)/2,1e-12,'ByRows',true); Oe=Oe';
else
  [~,Oe]=ismembertol(A.NodesCoordElem,(Ge.coordinates+1)/2,1e-12,'ByRows',true); Oe=Oe';
end
if nsd==2
  [~,Of]=ismembertol(A.NodesCoordFace,(Gf.coordinates'+1)/2,1e-12,'ByRows',true); Of=Of';
elseif nsd==3
  [~,Of]=ismembertol(A.NodesCoordFace,(Gf.coordinates+1)/2,1e-12,'ByRows',true); Of=Of';
end

% Element coordinates in [0,1]^nsd
if nsd==1
  Xe=(Fe.coordinates(Oe)'+1)/2;
else
  Xe=(Fe.coordinates(Oe,:)+1)/2;
end

% Face coordinates in [0,1]^(nsd-1)
if nsd==2
  Xf=(Ff.coordinates(Of)'+1)/2;
elseif nsd==3
  Xf=(Ff.coordinates(Of,:)+1)/2;
end

% Define symbolic variables
syms O I x y z

% Symbolic Pascal triangle
if nsd==1
  Pascal=[1,...
          x,...
          x^2,...
          x^3,...
          x^4,...
          x^5,...
          x^6,...
          x^7,...
          x^8];
elseif nsd==2
  Pascal=[1,...
          x,y,...
          x^2,y^2,x*y,...
          x^3,y^3,x^2*y,x*y^2,...
          x^4,y^4,x^3*y,x^2*y^2,x*y^3,...
          x^5,y^5,x^4*y,x^3*y^2,x^2*y^3,x*y^4,...
          x^6,y^6,x^5*y^1,x^4*y^2,x^3*y^3,x^2*y^4,x*y^5,...
          x^7,y^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,...
          x^8,y^8,x^7*y,x^6*y^2,x^5*y^3,x^4*y^4,x^3*y^5,x^2*y^6,x*y^7];
elseif nsd==3
  Pascal=[1,...
          x,y,z,...
          x^2,y^2,z^2,x*y,y*z,x*z,...
          x^3,y^3,z^3,x^2*y,x*y^2,y^2*z,y*z^2,x^2*z,x*z^2,x*y*z,...
          x^4,y^4,z^4,x^3*y,x^2*y^2,x*y^3,y^3*z,y^2*z^2,y*z^3,x^3*z,x^2*z^2,x*z^3,x^2*y*z,x*y^2*z,x*y*z^2,...
          x^5,y^5,z^5,x^4*y,x^3*y^2,x^2*y^3,x*y^4,y^4*z,y^3*z^2,y^2*z^3,y*z^4,x^4*z,x^3*z^2,x^2*z^3,x*z^4,x^3*y*z,x^2*y^2*z,x*y^3*z,x*y^2*z^2,x*y*z^3,x^2*y*z^2,...
          x^6,y^6,z^6,x^5*y^1,x^4*y^2,x^3*y^3,x^2*y^4,x*y^5,y^5*z,y^4*z^2,y^3*z^3,y^2*z^4,y*z^5,x^5*z,x^4*z^2,x^3*z^3,x^2*z^4,x*z^5,x^4*y*z,x^3*y^2*z,x^2*y^3*z,x*y^4*z,x*y^3*z^2,x*y^2*z^3,x*y*z^4,x^2*y*z^3,x^3*y*z^2,x^2*y^2*z^2,...
          x^7,y^7,z^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^6*z,y^5*z^2,y^4*z^3,y^3*z^4,y^2*z^5,y*z^6,x^6*z,x^5*z^2,x^4*z^3,x^3*z^4,x^2*z^5,x*z^6,x^5*y*z,x^4*y^2*z,x^3*y^3*z,x^2*y^4*z,x*y^5*z,x*y^4*z^2,x*y^3*z^3,x*y^2*z^4,x*y*z^5,x^2*y*z^4,x^3*y*z^3,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x^2*y^2*z^3,...
          x^8,y^8,z^8,x^7*y,x^6*y^2,x^5*y^3,x^4*y^4,x^3*y^5,x^2*y^6,x*y^7,y^7*z,y^6*z^2,y^5*z^3,y^4*z^4,y^3*z^5,y^2*z^6,y*z^7,x^7*z,x^6*z^2,x^5*z^3,x^4*z^4,x^3*z^5,x^2*z^6,x*z^7,x^6*y*z,x^5*y^2*z,x^4*y^3*z,x^3*y^4*z,x^2*y^5*z,x*y^6*z,x*y^5*z^2,x*y^4*z^3,x*y^3*z^4,x*y^2*z^5,x*y*z^6,x^2*y*z^5,x^3*y*z^4,x^4*y*z^3,x^5*y*z^2,x^4*y^2*z^2,x^3*y^3*z^2,x^2*y^4*z^2,x^2*y^3*z^3,x^2*y^2*z^4,x^3*y^2*z^3];
end

% Modified Pascal triangle and derivatives
nen=size(Xe,1);
P=Pascal(1:nen); P(P==0)=O;   P(P==1)=I;
Px=diff(P,x);    Px(Px==0)=O; Px(Px==1)=I;
Py=diff(P,y);    Py(Py==0)=O; Py(Py==1)=I;
Pz=diff(P,z);    Pz(Pz==0)=O; Pz(Pz==1)=I;

% % Compute shape functions coefficients
x=mp(Xe(:,1)); if nsd>1; y=mp(Xe(:,2)); end; if nsd>2; z=mp(Xe(:,3)); end; O=0*x; I=1+O; 
C=mp(eval(P))\eye(nen);

% Compute shape functions and derivatives
N=mp(eval(P))*C;
Nx=mp(eval(Px))*C;
Ny=mp(eval(Py))*C;
Nz=mp(eval(Pz))*C;

% Check shape functions at nodes (identity matrix)
error=max(max(abs(eval(P)*C-eye(nen))))

% Create strings
Xs =['[',strjoin(arrayfun(@(i) sprintf(repmat(' %.15f',1,size(Xe,2)),Xe(i,:)),1:size(Xe,1),'UniformOutput',0),';'),']']; Xs(2)=[]
Cs =['[',strjoin(arrayfun(@(i) sprintf(repmat(' %.15f',1,size(C ,2)),C(i ,:)),1:size(C ,1),'UniformOutput',0),';'),']']; Cs(2)=[]
Ps =sprintf('%c',replace(char(P ),{' ','*','/','^',','},{'','.*','./','.^',', '})); Ps( [1:8,end-1:end])=[]
Pxs=sprintf('%c',replace(char(Px),{' ','*','/','^',','},{'','.*','./','.^',', '})); Pxs([1:8,end-1:end])=[]
Pys=sprintf('%c',replace(char(Py),{' ','*','/','^',','},{'','.*','./','.^',', '})); Pys([1:8,end-1:end])=[]
Pzs=sprintf('%c',replace(char(Pz),{' ','*','/','^',','},{'','.*','./','.^',', '})); Pzs([1:8,end-1:end])=[]

% Write to file
fileId=fopen('temp.txt','w');
fprintf(fileId,'X=%s\n\n',Xs);
fprintf(fileId,'C=%s\n\n',Cs);
fprintf(fileId,'P=%s\n\n',Ps);
fprintf(fileId,'Px=%s\n\n',Pxs);
fprintf(fileId,'Py=%s\n\n',Pys);
fprintf(fileId,'Pz=%s\n\n',Pzs);
fclose(fileId);
fprintf('\nCOMPLETED\n');

% Plot shape functions
if nsd==1
  figure('color','w')
  x=(0:0.001:1)'; O=zeros(size(x)); I=ones(size(x));
  plot(x,eval(P)*C);
  hold on;
  plot(Xe,0,'ko');
  x=Xe; O=zeros(size(x)); I=ones(size(x));
  plot(Xe,eval(P)*C,'ro');
  xlim([0,1]);
  syms x O I
end