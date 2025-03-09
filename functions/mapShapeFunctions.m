function [N,Nx,Ny,Nz,w,Ng,pinvN]=mapShapeFunctions(...
         Where,RefElementGeometry,RefElementApproxim,X,nsd)
         % Map shape functions

% Inizialize
N=[]; Nx=[]; Ny=[]; Nz=[]; w=[]; Ng=[]; pinvN=[];

if Where==1
  
  % Element calculations
  Ng=RefElementGeometry.ShapeFunctionsElem;
  N=RefElementApproxim.ShapeFunctionsElem;
  if size(X,1)==nsd+1
    if nsd==2
      Ngxi= [-1,1,0];
      Ngeta=[-1,0,1];
    elseif nsd==3
      Ngxi=  [0,0,-1,1];
      Ngeta= [1,0,-1,0];
      Ngzeta=[0,1,-1,0];
    end
  else
    Ngxi=RefElementGeometry.ShapeFunctionDerivativesElemXi;
    Ngeta=RefElementGeometry.ShapeFunctionDerivativesElemEta;
    if nsd==3
      Ngzeta=RefElementGeometry.ShapeFunctionDerivativesElemZeta;
    end
  end
  Naxi=RefElementApproxim.ShapeFunctionDerivativesElemXi;
  Naeta=RefElementApproxim.ShapeFunctionDerivativesElemEta;
  if nsd==3
    Nazeta=RefElementApproxim.ShapeFunctionDerivativesElemZeta;
  end
  g=RefElementApproxim.GaussWeigthsElem;
  pinvN=RefElementApproxim.PseudoinverseShapeFunctionsElem;
  if nsd==2
    J11=Ngxi *X(:,1); J12=Ngxi *X(:,2);
    J21=Ngeta*X(:,1); J22=Ngeta*X(:,2);
    detJ=J11.*J22-J12.*J21;
    invJ11=+J22./detJ;
    invJ12=-J12./detJ;
    invJ21=-J21./detJ;
    invJ22=+J11./detJ;
    Nx=invJ11.*Naxi+invJ12.*Naeta;
    Ny=invJ21.*Naxi+invJ22.*Naeta;
  elseif nsd==3
    J11=Ngxi  *X(:,1); J12=Ngxi  *X(:,2); J13=Ngxi  *X(:,3);
    J21=Ngeta *X(:,1); J22=Ngeta *X(:,2); J23=Ngeta *X(:,3);
    J31=Ngzeta*X(:,1); J32=Ngzeta*X(:,2); J33=Ngzeta*X(:,3);
    detJ=J11.*J22.*J33-J11.*J23.*J32-J12.*J21.*J33...
        +J12.*J23.*J31+J13.*J21.*J32-J13.*J22.*J31;
    invJ11=(J22.*J33-J23.*J32)./detJ;
    invJ12=(J13.*J32-J12.*J33)./detJ;
    invJ13=(J12.*J23-J13.*J22)./detJ;
    invJ21=(J23.*J31-J21.*J33)./detJ;
    invJ22=(J11.*J33-J13.*J31)./detJ;
    invJ23=(J13.*J21-J11.*J23)./detJ;
    invJ31=(J21.*J32-J22.*J31)./detJ;
    invJ32=(J12.*J31-J11.*J32)./detJ;
    invJ33=(J11.*J22-J12.*J21)./detJ;
    Nx=invJ11.*Naxi+invJ12.*Naeta+invJ13.*Nazeta;
    Ny=invJ21.*Naxi+invJ22.*Naeta+invJ23.*Nazeta;
    Nz=invJ31.*Naxi+invJ32.*Naeta+invJ33.*Nazeta;
  end
  w=g.*detJ;
  
elseif Where==0
  
  % Face calculations
  Ng=RefElementGeometry.ShapeFunctionsFace;
  N=RefElementApproxim.ShapeFunctionsFace;
  if size(X,1)==nsd
    if nsd==2
      Ngxi= [-1,1];
    elseif nsd==3
      Ngxi= [-1,1,0];
      Ngeta=[-1,0,1];
    end
  else
    Ngxi=RefElementGeometry.ShapeFunctionDerivativesFaceXi;
    if nsd==3
      Ngeta=RefElementGeometry.ShapeFunctionDerivativesFaceEta;
    end
  end
  g=RefElementApproxim.GaussWeigthsFace;
  pinvN=RefElementApproxim.PseudoinverseShapeFunctionsFace;
  if nsd==2
    tx=Ngxi*X(:,1); ty=Ngxi*X(:,2);
    detJ=sqrt(tx.^2+ty.^2);
    Nx=+ty./detJ;
    Ny=-tx./detJ;
  elseif nsd==3
    t1x=Ngxi *X(:,1); t1y=Ngxi *X(:,2); t1z=Ngxi *X(:,3);
    t2x=Ngeta*X(:,1); t2y=Ngeta*X(:,2); t2z=Ngeta*X(:,3);
    tx=t1y.*t2z-t1z.*t2y;
    ty=t1z.*t2x-t1x.*t2z;
    tz=t1x.*t2y-t1y.*t2x;
    detJ=sqrt(tx.^2+ty.^2+tz.^2);
    Nx=tx./detJ;
    Ny=ty./detJ;
    Nz=tz./detJ;
  end
  w=g.*detJ;
  
end

end