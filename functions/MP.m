function [ONE]=MP()
         % Return double(1) or mp(1)

try mp;
  if mp.Digits~=16
    ONE=mp(1);
  else
    ONE=1;
  end
catch
  ONE=1;
end