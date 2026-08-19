% Run the whole test suite and report the outcome of each test

% Reset the environment and hide the figures while the tests run
init;
set(0,'DefaultFigureVisible','off');

% Run every test whose name starts with the name of this script
Files={dir(sprintf('tests/%s*.m',dbstack().name)).name}';
Passed=zeros(size(Files));
CPUTime=zeros(size(Files));
for Test=1:length(Files)
  CPUTime(Test)=cputime;
  run(['tests/',Files{Test}]);
  Passed(Test)=Simulation.TestPassed;
  CPUTime(Test)=cputime-CPUTime(Test);
end

% Print the results
clc;
for Test=1:length(Files)
  fprintf('\n%s%s',Files{Test},repmat('.',1,66-length(Files{Test})));
  if Passed(Test)
    fprintf('PASSED in %.1f sec',CPUTime(Test));
  else
    fprintf('FAILED!!!');
  end
end
Outcome={'failed','passed'};
fprintf('\n\nTest completed (%s) in %.0f sec (CPU time)\n',Outcome{all(Passed)+1},sum(CPUTime));

% Show and then close the figures
set(0,'DefaultFigureVisible','on');
Figures=findobj('type','figure');
set(Figures,'Visible','on');
for iFig=1:length(Figures)
  pause(0.2); close;
end

clearvars -except Files Passed CPUTime;
