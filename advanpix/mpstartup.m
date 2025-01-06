function mpstartup()
%MPSTARTUP Set default settings for Multiprecision Computing Toolbox.
%
%   Loads every time on toolbox startup before any other routine runs.
% 
%   Configures toolbox with default settings. 
%   Can be modified to change default precision, guard digits and other options.
% 
%   Does not need to be called explicitly by the user.
%
%   See also mp.Init

%   Copyright (c) 2008 - 2022 Advanpix LLC.    

    % Initialize Multiprecision Computing Toolbox
    mp.Init();

    % Set up quadruple precision by default
    mp.Digits(34);

    % Set up guard digits, 'mp' uses mp.Digits() + mp.GuardDigits() decimal digits in all calculations.
    mp.GuardDigits(0);
    
    % Set up maximum number of threads to use in computations enabled with multi-core parallelism:
    %
    % mp.NumberOfThreads(N) with:
    %
    % N = maxNumCompThreads [default]
    % 
    %    Use MATLAB's default setting. This is probably optimal strategy for most of the users,
    %    especially if user already controls number of computational threads by built-in MATLAB's commands.
    %
    %    Usually MATLAB assigns number of computational threads equal to number
    %    of physical cores of CPU. Also each thread is binned to one physical
    %    core using thread affinity control. This is optimal for most heavy
    %    mathematical computations.
    %
    %    Multi-threading in toolbox/MATLAB is based on OpenMP framework. 
    %    OpenMP allows flexible control on number of threads, affinity, 
    %    timings, thread scheduling, etc. using OS environment variables.
    %    See OpenMP specifications for details: https://www.openmp.org/specifications/
    %
    %    On Windows we rely on Intel OpenMP which allows even more detailed
    %    tuning with (KMP_) environment variables: https://software.intel.com/en-us/node/522775
    %
    %    For example, OpenMP configuration for 16-cores/32-threads Intel Core i9 7960X 
    %    might look like (for Windows):
    %
    %       OMP_NUM_THREADS = 16
    %       KMP_AFFINITY = explicit,proclist=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]
    %    
    %    Here we restrict OpenMP to use only 16 computational threads and attach every thread 
    %    to distinct physical CPU core ( = to use logical cores with even IDs).  
    % 
    % N = 0:
    %    Use all available CPU cores in the system (including logical cores).
    %    Pushes CPU utilization to the limit. This mode is beneficial for computations 
    %    which is massively parallel, e.g. for matrix multiplication.
    %    In other cases speed can actually go down. Please test this mode for your
    %    particular code to see if it provides any speed-up.
    % 
    %    Use with caution since this mode might slow down all other applications running on the computer.
    %    Suitable for users who run one instance of toolbox at a time (no parfor).
    % 
    % N ~= 0:
    %    Toolbox will use exactly N cores (including logical cores). No thread affinity is applied.
    %
    %    This is useful if you run several toolbox instances (e.g. with
    %    parfor) and want to balance the load among all workers.
    %
    %    In this case compute number of threads as:
    %
    %        N = total_number_of_cores / number_of_matlab_workers
    %
    % NOTE. The 'maxNumCompThreads' was declared deprecated starting from R2011 up to R2014. 
    %       But status of this command was restored later and now it is valid command.
    %
    N = 0;    
    if(exist('maxNumCompThreads') ~= 0) %#ok<EXIST>
        warning('off', 'MATLAB:maxNumCompThreads:Deprecated');
        N = maxNumCompThreads;
        warning('on', 'MATLAB:maxNumCompThreads:Deprecated');
    end
    mp.NumberOfThreads(N);  % <-- Change N here if you want to use non-default number of computational threads.
    
    % Auto-extend accuracy of built-in 'double' precision constants (pi, exp, sqrt(2), etc.)
    mp.ExtendConstAccuracy(false);

    % Follow MATLAB's numeric format preferences by default.
    % All digits are shown only in case of long* formats.     
    mp.FollowMatlabNumericFormat(false);

    % Override built-in array-creation routines (zeros, ones, eye, magic, etc.) to
    % produce multiprecision arrays instead of 'double' or 'single'.
    mp.OverrideDoubleBasicArrays(false);
    
    % Implicit arithmetic expansion is enabled by default to match
    % standard behavior of MATLAB after R2016b version.
    mp.EnableImplicitExpansion(true);
end

