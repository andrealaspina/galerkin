%% galerkin test 

% Setup
clear; clc; path(pathdef); close('all');
addpath('./');
addpath('./formulations');
addpath('./functions');
addpath('./geometry');
addpath('./input');
addpath('./output');
addpath('./symbolic');
addpath('./tests');
set(0,'DefaultFigureVisible','off');

% Run tests
TestName=dbstack().name;
Files={dir(sprintf('tests/%s*.m',TestName)).name}';
Passed=zeros(size(Files));
CPUTime=zeros(size(Files));
for Test=1:length(Files)
  CPUTime(Test)=cputime;
  run(['tests/',Files{Test}]);
  Passed(Test)=Simulation.TestPassed;
  CPUTime(Test)=cputime-CPUTime(Test);
end

% Print results
clc;
Exitus='passed';
for Test=1:length(Files)
  fprintf('\n%s%s',Files{Test},repmat('.',1,66-length(Files{Test})));
  if Passed(Test)
    fprintf('PASSED in %.1f sec',CPUTime(Test));
  else
    fprintf('FAILED!!!');
    Exitus='failed';
  end
end
fprintf('\n\nTest completed (%s) in %.0f sec (CPU time)\n',Exitus,sum(CPUTime));

% Show and then close figures
set(0,'DefaultFigureVisible','on');
set(findobj('type','figure'),'Visible','on');
for iFig=1:length(findobj('type','figure'))
  pause(0.2); close;
end

% Delete useless variables
clearvars -except Files Passed CPUTime;