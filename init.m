%% galerkin initialization
%
% Clears the workspace and adds the framework folders to the path.

clear; clc; close('all');
path(pathdef);
addpath('./','./formulations','./functions','./geometry','./input','./output','./tests');
set(0,'DefaultFigureVisible','on');
