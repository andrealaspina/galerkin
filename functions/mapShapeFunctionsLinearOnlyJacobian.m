function [N,Nx,Ny,Nz,w]=mapShapeFunctionsLinearOnlyJacobian(...
         Where,RefElement,Xdeg1,nsd)
         % Map shape functions (linear faces)

% Inizialize
Nz=[];

if Where==1
  
  % Element calulcations
  N=RefElement.ShapeFunctionsElem;
  g=RefElement.GaussWeigthsElem;
  if nsd==2
    x1=Xdeg1(1,1); y1=Xdeg1(1,2);
    x2=Xdeg1(2,1); y2=Xdeg1(2,2);
    x3=Xdeg1(3,1); y3=Xdeg1(3,2);
    lf1=sqrt((x1-x2)^2+(y1-y2)^2);
    lf2=sqrt((x2-x3)^2+(y2-y3)^2);
    lf3=sqrt((x3-x1)^2+(y3-y1)^2);
    sp=(lf1+lf2+lf3)/2;
    A=sqrt(sp*(sp-lf1)*(sp-lf2)*(sp-lf3));
    detJ=2*A;
    Nx=[];
    Ny=[];
  elseif nsd==3
    detJ=dot(Xdeg1(2,:)-Xdeg1(1,:),cross(Xdeg1(3,:)-Xdeg1(1,:),Xdeg1(4,:)-Xdeg1(1,:)));
    Nx=[];
    Ny=[];
  end
  w=g.*detJ;
  
else
  
  % Face calulcations
  N=RefElement.ShapeFunctionsFace;
  g=RefElement.GaussWeigthsFace;
  if nsd==2
    x1=Xdeg1(1,1); y1=Xdeg1(1,2);
    x2=Xdeg1(2,1); y2=Xdeg1(2,2);
    L=sqrt((x1-x2)^2+(y1-y2)^2);
    detJ=1*L;
    Nx=+(y2-y1);
    Ny=-(x2-x1);
    Nnorm=sqrt(Nx^2+Ny^2);
    Nx=Nx/Nnorm;
    Ny=Ny/Nnorm;
  elseif nsd==3
    x1=Xdeg1(1,1); y1=Xdeg1(1,2); z1=Xdeg1(1,3);
    x2=Xdeg1(2,1); y2=Xdeg1(2,2); z2=Xdeg1(2,3);
    x3=Xdeg1(3,1); y3=Xdeg1(3,2); z3=Xdeg1(3,3);
    lf1=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
    lf2=sqrt((x2-x3)^2+(y2-y3)^2+(z2-z3)^2);
    lf3=sqrt((x3-x1)^2+(y3-y1)^2+(z3-z1)^2);
    sp=(lf1+lf2+lf3)/2;
    A=sqrt(sp*(sp-lf1)*(sp-lf2)*(sp-lf3));
    detJ=2*A;
    Nx=((y1-y2)*(z1+z2)+(y2-y3)*(z2+z3)+(y3-y1)*(z3+z1));
    Ny=((z1-z2)*(x1+x2)+(z2-z3)*(x2+x3)+(z3-z1)*(x3+x1));
    Nz=((x1-x2)*(y1+y2)+(x2-x3)*(y2+y3)+(x3-x1)*(y3+y1));
    Nnorm=sqrt(Nx^2+Ny^2+Nz^2);
    Nx=Nx/Nnorm;
    Ny=Ny/Nnorm;
    Nz=Nz/Nnorm;
  end
  w=g.*detJ;
  
end