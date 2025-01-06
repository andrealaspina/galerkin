close all; clear; clc;

run check_jacobian_normal_euler_equations.mlx;

%% Print "positive" and "negative" part of the Jacobian matrix

% Print "positive part" of the Jacobian matrix with singularity in n_x
fprintf('\n\nAn+_1')
Anp=A_n_p_1;
for i=1:size(Anp,1)
  for j=1:size(Anp,2)
    Anp_str=char(Anp(i,j));
    Anp_str=strrep(strrep(strrep(strrep(Anp_str,' ',''),'*','.*'),'/','./'),'^','.^');
    Anp_str=strrep(strrep(strrep(Anp_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
    Anp_str=strrep(strrep(strrep(Anp_str,'v_n','vn'),'e_k','ek'),'h_0','ho');
    Anp_str=strrep(strrep(Anp_str,'2.*','2*'),'4.*','4*');
    Anp_str=strrep(strrep(Anp_str,'./2','/2'),'./4','/4');
    Anp_str=strrep(Anp_str,'^2*','^2.*');
    fprintf('\nAnp(:,%d,%d)=%s;',i,j,Anp_str)
  end
end

% Print "negative part" of the Jacobian matrix with singularity in n_x
fprintf('\n\nAn-_1')
Anm=A_n_m_1;
for i=1:size(Anm,1)
  for j=1:size(Anm,2)
    Anm_str=char(Anm(i,j));
    Anm_str=strrep(strrep(strrep(strrep(Anm_str,' ',''),'*','.*'),'/','./'),'^','.^');
    Anm_str=strrep(strrep(strrep(Anm_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
    Anm_str=strrep(strrep(strrep(Anm_str,'v_n','vn'),'e_k','ek'),'h_0','ho');
    Anm_str=strrep(strrep(Anm_str,'2.*','2*'),'4.*','4*');
    Anm_str=strrep(strrep(Anm_str,'./2','/2'),'./4','/4');
    Anm_str=strrep(Anm_str,'^2*','^2.*');
    fprintf('\nAnm(:,%d,%d)=%s;',i,j,Anm_str)
  end
end

if nsd==3
  
  % Print "positive part" of the Jacobian matrix with singularity in n_y
  fprintf('\n\nAn+_2')
  Anp=A_n_p_2;
  for i=1:size(Anp,1)
    for j=1:size(Anp,2)
      Anp_str=char(Anp(i,j));
      Anp_str=strrep(strrep(strrep(strrep(Anp_str,' ',''),'*','.*'),'/','./'),'^','.^');
      Anp_str=strrep(strrep(strrep(Anp_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
      Anp_str=strrep(strrep(strrep(Anp_str,'v_n','vn'),'e_k','ek'),'h_0','ho');
      Anp_str=strrep(strrep(Anp_str,'2.*','2*'),'4.*','4*');
      Anp_str=strrep(strrep(Anp_str,'./2','/2'),'./4','/4');
      Anp_str=strrep(Anp_str,'^2*','^2.*');
      fprintf('\nAnp(:,%d,%d)=%s;',i,j,Anp_str)
    end
  end
  
  % Print "negative part" of the Jacobian matrix with singularity in n_y
  fprintf('\n\nAn-_2')
  Anm=A_n_m_2;
  for i=1:size(Anm,1)
    for j=1:size(Anm,2)
      Anm_str=char(Anm(i,j));
      Anm_str=strrep(strrep(strrep(strrep(Anm_str,' ',''),'*','.*'),'/','./'),'^','.^');
      Anm_str=strrep(strrep(strrep(Anm_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
      Anm_str=strrep(strrep(strrep(Anm_str,'v_n','vn'),'e_k','ek'),'h_0','ho');
      Anm_str=strrep(strrep(Anm_str,'2.*','2*'),'4.*','4*');
      Anm_str=strrep(strrep(Anm_str,'./2','/2'),'./4','/4');
      Anm_str=strrep(Anm_str,'^2*','^2.*');
      fprintf('\nAnm(:,%d,%d)=%s;',i,j,Anm_str)
    end
  end
  
  % Print "positive part" of the Jacobian matrix with singularity in n_z
  fprintf('\n\nAn+_3')
  Anp=A_n_p_3;
  for i=1:size(Anp,1)
    for j=1:size(Anp,2)
      Anp_str=char(Anp(i,j));
      Anp_str=strrep(strrep(strrep(strrep(Anp_str,' ',''),'*','.*'),'/','./'),'^','.^');
      Anp_str=strrep(strrep(strrep(Anp_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
      Anp_str=strrep(strrep(strrep(Anp_str,'v_n','vn'),'e_k','ek'),'h_0','ho');
      Anp_str=strrep(strrep(Anp_str,'2.*','2*'),'4.*','4*');
      Anp_str=strrep(strrep(Anp_str,'./2','/2'),'./4','/4');
      Anp_str=strrep(Anp_str,'^2*','^2.*');
      fprintf('\nAnp(:,%d,%d)=%s;',i,j,Anp_str)
    end
  end
  
  % Print "negative part" of the Jacobian matrix with singularity in n_z
  fprintf('\n\nAn-_3')
  Anm=A_n_m_3;
  for i=1:size(Anm,1)
    for j=1:size(Anm,2)
      Anm_str=char(Anm(i,j));
      Anm_str=strrep(strrep(strrep(strrep(Anm_str,' ',''),'*','.*'),'/','./'),'^','.^');
      Anm_str=strrep(strrep(strrep(Anm_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
      Anm_str=strrep(strrep(strrep(Anm_str,'v_n','vn'),'e_k','ek'),'h_0','ho');
      Anm_str=strrep(strrep(Anm_str,'2.*','2*'),'4.*','4*');
      Anm_str=strrep(strrep(Anm_str,'./2','/2'),'./4','/4');
      Anm_str=strrep(Anm_str,'^2*','^2.*');
      fprintf('\nAnm(:,%d,%d)=%s;',i,j,Anm_str)
    end
  end
  
end

%% Print linearization of "positive" and "negative" part of the Jacobian matrix

% Print "positive part" of the Jacobian matrix with singularity in n_x
fprintf('\n\ndAn+dU_1')
dAnpdU=dA_n_p_1_dU;
for i=1:size(dAnpdU,1)
  for j=1:size(dAnpdU,2)
    for k=1:size(dAnpdU,3)
      dAnpdU_str=char(dAnpdU(i,j,k));
      dAnpdU_str=strrep(strrep(strrep(strrep(dAnpdU_str,' ',''),'*','.*'),'/','./'),'^','.^');
      dAnpdU_str=strrep(strrep(strrep(dAnpdU_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
      dAnpdU_str=strrep(strrep(strrep(dAnpdU_str,'w_x','wx'),'w_y','wy'),'w_z','wz');
      fprintf('\ndAnpdU(:,%d,%d,%d)=%s;',i,j,k,dAnpdU_str)
    end
  end
end

% Print "negative part" of the Jacobian matrix with singularity in n_x
fprintf('\n\ndAn-dU_1')
dAnmdU=dA_n_m_1_dU;
for i=1:size(dAnmdU,1)
  for j=1:size(dAnmdU,2)
    for k=1:size(dAnmdU,3)
      dAnmdU_str=char(dAnmdU(i,j,k));
      dAnmdU_str=strrep(strrep(strrep(strrep(dAnmdU_str,' ',''),'*','.*'),'/','./'),'^','.^');
      dAnmdU_str=strrep(strrep(strrep(dAnmdU_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
      dAnmdU_str=strrep(strrep(strrep(dAnmdU_str,'w_x','wx'),'w_y','wy'),'w_z','wz');
      fprintf('\ndAnmdU(:,%d,%d,%d)=%s;',i,j,k,dAnmdU_str)
    end
  end
end

if nsd==3

  % Print "positive part" of the Jacobian matrix with singularity in n_y
  fprintf('\n\ndAn+dU_2')
  dAnpdU=dA_n_p_2_dU;
  for i=1:size(dAnpdU,1)
    for j=1:size(dAnpdU,2)
      for k=1:size(dAnpdU,3)
        dAnpdU_str=char(dAnpdU(i,j,k));
        dAnpdU_str=strrep(strrep(strrep(strrep(dAnpdU_str,' ',''),'*','.*'),'/','./'),'^','.^');
        dAnpdU_str=strrep(strrep(strrep(dAnpdU_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
        dAnpdU_str=strrep(strrep(strrep(dAnpdU_str,'w_x','wx'),'w_y','wy'),'w_z','wz');
        fprintf('\ndAnpdU(:,%d,%d,%d)=%s;',i,j,k,dAnpdU_str)
      end
    end
  end
  
  % Print "negative part" of the Jacobian matrix with singularity in n_y
  fprintf('\n\ndAn-dU_2')
  dAnmdU=dA_n_m_2_dU;
  for i=1:size(dAnmdU,1)
    for j=1:size(dAnmdU,2)
      for k=1:size(dAnmdU,3)
        dAnmdU_str=char(dAnmdU(i,j,k));
        dAnmdU_str=strrep(strrep(strrep(strrep(dAnmdU_str,' ',''),'*','.*'),'/','./'),'^','.^');
        dAnmdU_str=strrep(strrep(strrep(dAnmdU_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
        dAnmdU_str=strrep(strrep(strrep(dAnmdU_str,'w_x','wx'),'w_y','wy'),'w_z','wz');
        fprintf('\ndAnmdU(:,%d,%d,%d)=%s;',i,j,k,dAnmdU_str)
      end
    end
  end
  
  % Print "positive part" of the Jacobian matrix with singularity in n_z
  fprintf('\n\ndAn+dU_3')
  dAnpdU=dA_n_p_3_dU;
  for i=1:size(dAnpdU,1)
    for j=1:size(dAnpdU,2)
      for k=1:size(dAnpdU,3)
        dAnpdU_str=char(dAnpdU(i,j,k));
        dAnpdU_str=strrep(strrep(strrep(strrep(dAnpdU_str,' ',''),'*','.*'),'/','./'),'^','.^');
        dAnpdU_str=strrep(strrep(strrep(dAnpdU_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
        dAnpdU_str=strrep(strrep(strrep(dAnpdU_str,'w_x','wx'),'w_y','wy'),'w_z','wz');
        fprintf('\ndAnpdU(:,%d,%d,%d)=%s;',i,j,k,dAnpdU_str)
      end
    end
  end
  
  % Print "negative part" of the Jacobian matrix with singularity in n_z
  fprintf('\n\ndAn-dU_3')
  dAnmdU=dA_n_m_3_dU;
  for i=1:size(dAnmdU,1)
    for j=1:size(dAnmdU,2)
      for k=1:size(dAnmdU,3)
        dAnmdU_str=char(dAnmdU(i,j,k));
        dAnmdU_str=strrep(strrep(strrep(strrep(dAnmdU_str,' ',''),'*','.*'),'/','./'),'^','.^');
        dAnmdU_str=strrep(strrep(strrep(dAnmdU_str,'n_x','nx'),'n_y','ny'),'n_z','nz');
        dAnmdU_str=strrep(strrep(strrep(dAnmdU_str,'w_x','wx'),'w_y','wy'),'w_z','wz');
        fprintf('\ndAnmdU(:,%d,%d,%d)=%s;',i,j,k,dAnmdU_str)
      end
    end
  end
  
end

fprintf('\n\n')