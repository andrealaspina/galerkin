function [X,W]=quadratureSimplex(N,vert)
% Gauss quadrature on simplices

MP1=(max(vert(:))-min(vert(:)))/(max(vert(:))-min(vert(:)));
m=size(vert,1);
nsd=size(vert,2)*MP1;
Nn=N^nsd;
if nsd==1 % 1-D simplex
  [q,w]=rquad(N,nsd*0);
  len=diff(vert);
  X=vert(1)+len*q;
  W=abs(len)*w;
else    % Higher-dimensional simplex
  q=cell(1,nsd); w=cell(1,nsd);
  for k=1:nsd
    [q{k},w{k}]=rquad(N,nsd-k);
  end
  [Q{1:nsd}]=ndgrid(q{:}); q=reshape(cat(nsd,Q{:}),Nn,nsd);
  [W{1:nsd}]=ndgrid(w{:}); w=reshape(cat(nsd,W{:}),Nn,nsd);
  map=eye(m); map(2:m,1)=-1; c=map*vert;
  qp=cumprod(q,2); e=ones(Nn,1);
  X=[e,[(1-q(:,1:nsd-1)),e].*[e,qp(:,1:nsd-2),qp(:,nsd)]]*c;
  W=abs(det(c(2:m,:)))*prod(w,2);
end

function [x,w]=rquad(N,k)

k1=k+1; k2=k+2; n=1:N;  nnk=2*n+k;
A=[k/k2 repmat(k^2,1,N)./(nnk.*(nnk+2))];
n=2:N; nnk=nnk(n);
B1=4*k1/(k2*k2*(k+3)); nk=n+k; nnk2=nnk.*nnk;
B=4*(n.*nk).^2./(nnk2.*nnk2-nnk2);
ab=[A' [(2^k1)/k1; B1; B']]; s=sqrt(ab(2:N,2));
[V,X]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
[X,I]=sort(diag(X));
x=(X+1)/2; w=(1/2)^(k1)*ab(1,2)*V(1,I)'.^2;