function Match=matchField(Structure,Field,Value)
% Check whether a structure has a given field and, optionally, a given value

Match=any(strcmp(fieldnames(Structure),Field));
if Match && nargin==3
  if ischar(Value)
    Match=strcmp(Structure.(Field),Value);
  else
    Match=isequal(Structure.(Field),Value);
  end
end

end
