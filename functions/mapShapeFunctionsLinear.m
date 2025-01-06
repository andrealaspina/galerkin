function [N,Nx,Ny,Nz,w]=mapShapeFunctionsLinear(...
         Where,RefElement,Xdeg1,nsd)
         % Map shape functions (linear faces)

% Inizialize
Nz=[];

if Where==1
  
  % Element calulcations
  N=RefElement.ShapeFunctionsElem;
  Nxi=RefElement.ShapeFunctionDerivativesElemXi;
  Neta=RefElement.ShapeFunctionDerivativesElemEta;
  if nsd==3
    Nzeta=RefElement.ShapeFunctionDerivativesElemZeta;
  end
  g=RefElement.GaussWeigthsElem;
  if nsd==2
    Nxideg1= [-1,1,0];
    Netadeg1=[-1,0,1];
    J11=Nxideg1*Xdeg1(:,1);  J12=Nxideg1*Xdeg1(:,2);
    J21=Netadeg1*Xdeg1(:,1); J22=Netadeg1*Xdeg1(:,2);
    detJ=J11*J22-J12*J21;
    invJ11=+J22/detJ;
    invJ12=-J12/detJ;
    invJ21=-J21/detJ;
    invJ22=+J11/detJ;
    
    Nx=invJ11*Nxi+invJ12*Neta;
    Ny=invJ21*Nxi+invJ22*Neta;
  elseif nsd==3
    Nxideg1=  [0,0,-1,1];
    Netadeg1= [1,0,-1,0];
    Nzetadeg1=[0,1,-1,0];
    J11=Nxideg1*Xdeg1(:,1);   J12=Nxideg1*Xdeg1(:,2);   J13=Nxideg1*Xdeg1(:,3);
    J21=Netadeg1*Xdeg1(:,1);  J22=Netadeg1*Xdeg1(:,2);  J23=Netadeg1*Xdeg1(:,3);
    J31=Nzetadeg1*Xdeg1(:,1); J32=Nzetadeg1*Xdeg1(:,2); J33=Nzetadeg1*Xdeg1(:,3);
    detJ=J11*J22*J33-J11*J23*J32-J12*J21*J33...
        +J12*J23*J31+J13*J21*J32-J13*J22*J31;
    invJ11=(J22*J33-J23*J32)/detJ;
    invJ12=(J13*J32-J12*J33)/detJ;
    invJ13=(J12*J23-J13*J22)/detJ;
    invJ21=(J23*J31-J21*J33)/detJ;
    invJ22=(J11*J33-J13*J31)/detJ;
    invJ23=(J13*J21-J11*J23)/detJ;
    invJ31=(J21*J32-J22*J31)/detJ;
    invJ32=(J12*J31-J11*J32)/detJ;
    invJ33=(J11*J22-J12*J21)/detJ;
    
    Nx=invJ11*Nxi+invJ12*Neta+invJ13*Nzeta;
    Ny=invJ21*Nxi+invJ22*Neta+invJ23*Nzeta;
    Nz=invJ31*Nxi+invJ32*Neta+invJ33*Nzeta;
  end
  w=g.*detJ;
  
else
  
  % Face calulcations
  N=RefElement.ShapeFunctionsFace;
  g=RefElement.GaussWeigthsFace;
  if nsd==2
    x1=Xdeg1(1,1); y1=Xdeg1(1,2);
    x2=Xdeg1(2,1); y2=Xdeg1(2,2);
    Nx=+(y2-y1);
    Ny=-(x2-x1);
    Nnorm=sqrt(Nx^2+Ny^2);
    Nx=Nx/Nnorm;
    Ny=Ny/Nnorm;
    L=sqrt((x1-x2)^2+(y1-y2)^2);
    detJ=1*L;
  elseif nsd==3
    x1=Xdeg1(1,1); y1=Xdeg1(1,2); z1=Xdeg1(1,3);
    x2=Xdeg1(2,1); y2=Xdeg1(2,2); z2=Xdeg1(2,3);
    x3=Xdeg1(3,1); y3=Xdeg1(3,2); z3=Xdeg1(3,3);
    lf1=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
    lf2=sqrt((x2-x3)^2+(y2-y3)^2+(z2-z3)^2);
    lf3=sqrt((x3-x1)^2+(y3-y1)^2+(z3-z1)^2);
    sp=(lf1+lf2+lf3)/2;
    Nx=((y1-y2)*(z1+z2)+(y2-y3)*(z2+z3)+(y3-y1)*(z3+z1));
    Ny=((z1-z2)*(x1+x2)+(z2-z3)*(x2+x3)+(z3-z1)*(x3+x1));
    Nz=((x1-x2)*(y1+y2)+(x2-x3)*(y2+y3)+(x3-x1)*(y3+y1));
    Nnorm=sqrt(Nx^2+Ny^2+Nz^2);
    Nx=Nx/Nnorm;
    Ny=Ny/Nnorm;
    Nz=Nz/Nnorm;
    A=sqrt(sp*(sp-lf1)*(sp-lf2)*(sp-lf3));
    detJ=2*A;
  end
  w=g.*detJ;
  
end