function [Match]=matchField(Structure,Field,Value)
% Check whether a structure field exists and, optionally, whether it matches a given value

Match=any(strcmp(fieldnames(Structure),Field));
if Match && nargin==3
  if isnumeric(Value)
    Match=isequal(Structure.(Field),Value);
  elseif ischar(Value)
    Match=strcmp(Structure.(Field),Value);
  else
    Match=false;
  end
end

end
