function [ONE]=MP()
% Return the scalar one in the active precision (double, or mp when Advanpix is in use)

try
  if mp.Digits~=16
    ONE=mp(1);
  else
    ONE=1;
  end
catch
  ONE=1;
end

end
