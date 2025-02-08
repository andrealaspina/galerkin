classdef Formulation
  
  properties (Access=private)
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) 1;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) 0;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) 1;
    
    % Number of post-processed components
    NumPostComp=@(NumSpaceDim) 1;

    % Discretization type
    DiscretizationType='';
    
    % Time derivative order
    TimeDerOrder=1;
    
    % Time/frequency domain
    Domain='Time';
    
  end
  
  methods
    
    % Initialize unknowns
    function [Block]=initializeUnknowns(~,~,Block,~,~,~)
      % Input: iD,Block,Parameters,Time,Sizes
    end
    
    % Compute initial conditions
    function [Block]=computeInitialConditions(~,~,Block,~,~,~,~,~,~)
      % Input: iD,Block,Parameters,Mesh,Faces,Time,RefElement,Sizes
    end

    % Evaluate solution at fixed DOFs
    function [Block]=evaluateSolutionFixedDofs(~,~,Block,~,~,~,~)
      % Input: iD,Block,Parameters,Mesh,Time,Sizes
    end
    
    % Build block
    function [Block]=buildBlock(~,~,Block,~,~,~,~,~,~,~,~)
      % Input: iD,Block,Elements,Simulation,Parameters,Mesh,Faces,Time,RefElement,Sizes
    end
    
    % Do coupling
    function [Block]=doCoupling(~,~,~,Block,~,~,~,~,~,~,~,~)
      % Input: iD1,iD2,Block,Elements,Simulation,Parameters,Mesh,Faces,Time,RefElement,Sizes
    end
    
    % Do post-process
    function [Block]=doPostProcess(~,~,~,~,~,~,~,~,~)
      % Input: Elements,Simulation,Parameters,Mesh,Faces,Time,RefElement,Sizes
      Block=[];
    end
    
    % Store results
    function [Results]=storeResults(~,~,~,Results,~,~,~,~,~,~,~)
      % Input: iD,iST,Results,Block,Simulation,Parameters,Mesh,Time,RefElement,Sizes
    end
    
  end
  
end