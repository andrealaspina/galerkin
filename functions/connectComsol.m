function []=connectComsol(NumProcessors)
         % Connect Comsol

% Add path to Comsol LiveLink for Matlab
addpath /Applications/COMSOL56/Multiphysics/mli
import com.comsol.model.*
import com.comsol.model.util.*

% Get number of processors
if nargin==0
  NumProcessors=1;
end

% Get number of workers
Pool=gcp;
if isempty(Pool)
  Pool=[];
  Pool.NumWorkers=1;
end

% Get current worker
Worker=getCurrentTask();
if isempty(Worker)
  Worker.ID=1;
end

% Get message
Message=getReport(MPHSTART,'basic','hyperlinks','off');

% Do different actions depending on the message
if     strcmp(Message(1:22),'Connection established')
  % Do nothing
elseif strcmp(Message(22:50),'Already connected to a server')
  % Do nothing
elseif strcmp(Message(22:68),'A connection to COMSOL could not be established')
  % Connect given worker
  Port=2035+Worker.ID;
  system(sprintf('/Applications/COMSOL56/Multiphysics/bin/comsol -np %d server -port %d &',...
    NumProcessors,Port));
  pause(Pool.NumWorkers+2);
  mphstart(Port);
end

end

function [Message]=MPHSTART()

% Throw message after trying mphstart()
try mphstart()
  Message=MException('','Connection established');
catch Message
end

end