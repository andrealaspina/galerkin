close all; clear; clc;

syms mu lambda real positive;
syms Fxx Fxy Fxz real
syms Fyx Fyy Fyz real
syms Fzx Fzy Fzz real

% Input
nsd=3;
Model='LinearElasticity';
Transform='no';

% Deformation gradient
if nsd==2
  F=[Fxx Fxy
     Fyx Fyy];
elseif nsd==3
  F=[Fxx Fxy Fxz
     Fyx Fyy Fyz
     Fzx Fzy Fzz];
end

% Stress
if strcmp(Model,'LinearElasticity')
  epsilon=((F-eye(nsd))+(F-eye(nsd))')/2;
  s_mat=2*mu*epsilon+lambda*trace(epsilon)*eye(nsd);
elseif strcmp(Model,'StVenantKirchhoff')
  E=(F'*F-eye(nsd))/2;
  psi=mu*sum(sum(E.*E))+lambda/2*(trace(E))^2;
  s_mat=mygradient(psi,F);
elseif strcmp(Model,'NeoHook')
  psi=mu/2*(trace(F'*F)-2*log(det(F))-nsd)+lambda/2*(log(det(F)))^2;
  s=mygradient(psi,reshape(F',1,nsd^2));
  s_mat=transpose(reshape(s,nsd,nsd));
end

% Transformation of stress
if strcmp(Transform,'no')
  
elseif strcmp(Transform,'push-forward')
  s_mat=det(F)^(-1)*s_mat*F';
elseif strcmp(Transform,'pull-back')
  s_mat=det(F)*s_mat*(F^(-1))';
end

%s_mat=simplify(s_mat);
s_mat=s_mat;
s=reshape(transpose(s_mat),1,nsd^2);

% Linearization of stress
%dsdF_mat=simplify(mygradient(s,reshape(F',1,nsd^2)));
dsdF_mat=mygradient(s,reshape(F',1,nsd^2));
dsdF=reshape(transpose(dsdF_mat),1,nsd^4);

% Collect terms
syms J real positive;
s=subs(s,det(F),J);
dsdF=subs(dsdF,det(F),J);
if strcmp(Model,'LinearElasticity')
  syms trEps real positive;
  s=subs(s,trace(epsilon),trEps);
  dsdF=subs(dsdF,trace(epsilon),trEps);
elseif strcmp(Model,'StVenantKirchhoff')
  syms trE real positive;
  s=subs(s,trace(E),trE);
  dsdF=subs(dsdF,trace(E),trE);
elseif strcmp(Model,'NeoHook')
  syms logJ real;
  s=subs(s,log(J),logJ);
  dsdF=subs(dsdF,log(J),logJ);
end

% Print terms
fprintf('\n\nCOLLECTED TERMS')
J_str=strrep(strrep(strrep(strrep(strrep(char(det(F)),...
  ' ',''),'*','.*'),'/','./'),'^','.^'),',',', ');
fprintf('\nJ=%s;',J_str)
if strcmp(Model,'LinearElasticity')
  trEps_str=strrep(strrep(strrep(strrep(strrep(char(trace(epsilon)),...
    ' ',''),'*','.*'),'/','./'),'^','.^'),',',', ');
  fprintf('\ntrEps=%s;',strrep(strrep(strrep(strrep(strrep(char(trace(epsilon)),...
    ' ',''),'*','.*'),'/','./'),'^','.^'),',',', '))
elseif strcmp(Model,'StVenantKirchhoff')
  trE_str=strrep(strrep(strrep(strrep(strrep(char(trace(E)),...
    ' ',''),'*','.*'),'/','./'),'^','.^'),',',', ');
  fprintf('\ntrE=%s;',trE_str)
elseif strcmp(Model,'NeoHook')
  logJ_str=strrep(strrep(strrep(strrep(strrep(char(log(J)),...
    ' ',''),'*','.*'),'/','./'),'^','.^'),',',', ');
  fprintf('\nlogJ=%s;',logJ_str)
end

% Simplify and print stress
fprintf('\n\nSTRESS')
for i=1:size(s,2)
  % Simplify
  s(:,i)=simplify(collect(s(:,i),[mu,lambda]),'Steps',1000);
  s(:,i)=subs(s(:,i),det(F),J);
  if strcmp(Model,'LinearElasticity')
    s(:,i)=collect(subs(s(:,i),trace(epsilon),trEps),trEps);
  elseif strcmp(Model,'StVenantKirchhoff')
    s(:,i)=collect(subs(s(:,i),trace(E),trE),trE);
  elseif strcmp(Model,'NeoHook')
    s(:,i)=collect(subs(s(:,i),log(J),logJ),logJ);
  end
  s(:,i)=simplify(s(:,i),'Steps',1000);
 

  % Print
  s_str=strrep(strrep(strrep(strrep(strrep(char(s(:,i)),...
    ' ',''),'*','.*'),'/','./'),'^','.^'),',',', ');
  fprintf('\ns(:,%1d)=%s;',i,s_str)
end

% Simplify and print linearization of stress
fprintf('\n\nLINEARIZATION OF STRESS')
for i=1:size(dsdF,2)
  % Simplify
  dsdF(:,i)=simplify(collect(dsdF(:,i),[mu,lambda]),'Steps',1000);
  dsdF(:,i)=subs(dsdF(:,i),det(F),J);
  if strcmp(Model,'LinearElasticity')
    dsdF(:,i)=collect(subs(dsdF(:,i),trace(epsilon),trEps),trEps);
  elseif strcmp(Model,'StVenantKirchhoff')
    dsdF(:,i)=collect(subs(dsdF(:,i),trace(E),trE),trE);
  elseif strcmp(Model,'NeoHook')
    dsdF(:,i)=collect(subs(dsdF(:,i),log(J),logJ),logJ);
  end
  dsdF(:,i)=simplify(dsdF(:,i),'Steps',1000);

  % Print
  dsdF_str=strrep(strrep(strrep(strrep(strrep(char(dsdF(:,i)),...
    ' ',''),'*','.*'),'/','./'),'^','.^'),',',', ');
  fprintf('\ndsdF(:,%2d)=%s;',i,dsdF_str)
end

fprintf('\n\n')