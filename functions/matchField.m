function [Match]=matchField(...
         Structure,Field,Value)
         % Check if a structure field exists and potentially matches a given value

Match=false;
switch nargin
  case 2 % Check if a structure field exists
    if any(strcmp(fieldnames(Structure),Field))
      Match=true;
    end
  case 3 % Check if a structure field exists and matches a given value
    if any(strcmp(fieldnames(Structure),Field))
      if     isa(Value,'numeric') && Structure.(Field)==Value
        Match=true;
      elseif isa(Value,'char')    && strcmp(Structure.(Field),Value)
        Match=true;
      end
    end
end

end