close all; clear; clc;

k=1
         
%Note: k<=8
if(k>8)
  error('Please enter a value for k less than 9');
end

r=k+1;
X=0:(1/k):1;
syms x y delta L N S real;
delta=1-x-y;
L=sym(ones(1,r));

% The following loop generates the Lagrange coefficients
for i=2:r
  for j=1:i
    if(i~=j)
      vpa(L);
      L(i)=(L(i)*(x-X(j))/(X(i)-X(j)));
    end
  end
end

% The following MATLAB routines generate the shape functions for any ordered polynomial
nodes=(k+1)*(k+2)/2; % total no. of nodes
n=1;
index=zeros(1,3,nodes);
for j=0:k
  for i=0:k-j
    index(:,:,n)=[i,j,k-(i+j)]; % Area coordinates for each nodes
    S(n)=expand(L(i+1)*subs(L(j+1),y)*subs(L(k-(i+j)+1),delta));% Shape functions for all the nodes
    n=n+1;
  end
end

%The follwoing MATLAB routines are used to generate Lagrange Shape functions in Anticlockwise sequence for up to k=8, and
%the same logic can be used to extend the generation of shape functions for any order.
N(1)=S(k+1);
N(2)=S(nodes);
N(3)=S(1);
a=k;t=k+1;
for i=4:k+2 % 1st hypot
  t=t+a;
  N(i)=S(t);
  a=a-1;
end
c=-2;t=nodes;
for i=(k+3):(2*k+1) % 1st vert
  t=t+c ;
  N(i)=S(t) ;
  c=c-1;
end
for i=(2*k+2):(3*k) % 1st horiz
  t=i-2*k;
  N(i)=S(t);
end
a=k;
for i=(3*k+1):(4*k-2) % 2nd hypot
  t=t+a;
  N(i)=S(t);
  a=a-1;
end
c=-4;
for i=(4*k-1):(5*k-5) % 2nd vert
  t=t+c;
  N(i)=S(t);
  c=c-1;
end
for i=(5*k-4):(6*k-9) % 2nd horiz
  t=t+1;
  N(i)=S(t);
end
a=k-1;
for i=(6*k-8):(7*k-14) % 3nd hypot
  t=t+a;
  N(i)=S(t);
  a=a-1;
end
c=-6;
for i=(7*k-13):(8*k-20) % 3rd vert
  t=t+c;
  N(i)=S(t);
  c=c-1;
end
for i=(8*k-19):(9*k-27) % 3rd horiz
  t=t+1;
  N(i)=S(t);
end

% MY CHANGES ---------------------------------------------------------------------------------------
% Adapt nodes ordering
N(1:3)=circshift(N(1:3),1);
N(4:3*k)=circshift(N(4:3*k),k-1);

% Derivatives of the shape functions
NXi=diff(N,x);
NEta=diff(N,y);

% Vectors of ones and zeros
syms I O real

% Create strings
N(abs(N)==1)=I*sign(N(abs(N)==1)); N(N==0)=O;
N=['[',strjoin(arrayfun(@char,N,'uniform',0),', '),']'];
N=strrep(strrep(strrep(strrep(strrep(N,' ',''),'*','.*'),'/','./'),'^','.^'),',',', ')

NXi(abs(NXi)==1)=I*sign(NXi(abs(NXi)==1)); NXi(NXi==0)=O;
NXi=['[',strjoin(arrayfun(@char,NXi,'uniform',0),', '),']'];
NXi=strrep(strrep(strrep(strrep(strrep(NXi,' ',''),'*','.*'),'/','./'),'^','.^'),',',', ')

NEta(abs(NEta)==1)=I*sign(NEta(abs(NEta)==1)); NEta(NEta==0)=O;
NEta=['[',strjoin(arrayfun(@char,NEta,'uniform',0),', '),']'];
NEta=strrep(strrep(strrep(strrep(strrep(NEta,' ',''),'*','.*'),'/','./'),'^','.^'),',',', ')
% --------------------------------------------------------------------------------------------------