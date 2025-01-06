close all; clear; clc;

syms mu lambda real positive;
syms Fxx Fxy Fxz real
syms Fyx Fyy Fyz real
syms Fzx Fzy Fzz real

% Input
nsd=2;
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

s_mat=simplify(s_mat);
s=reshape(transpose(s_mat),1,nsd^2);

% Linearization of stress
dsdF_mat=simplify(mygradient(s,reshape(F',1,nsd^2)));
dsdF=reshape(transpose(dsdF_mat),1,nsd^4);

% Print stress
fprintf('\n\nSTRESS')
for i=1:size(s,2)
  s_str=strrep(strrep(strrep(strrep(strrep(char(s(:,i)),...
    ' ',''),'*','.*'),'/','./'),'^','.^'),',',', ');
  fprintf('\ns(:,%d)=%s;',i,s_str)
end

% Print linearization of stress
fprintf('\n\nLINEARIZATION OF STRESS')
for i=1:size(dsdF,2)
  dsdF_str=strrep(strrep(strrep(strrep(strrep(char(dsdF(:,i)),...
    ' ',''),'*','.*'),'/','./'),'^','.^'),',',', ');
  fprintf('\ndsdF(:,%d)=%s;',i,dsdF_str)
end

fprintf('\n\n')