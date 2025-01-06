function varargout = cell2mat(varargin)
%CELL2MAT Convert the contents of a cell array into a single matrix.
%
%   Identical to built-in function CELL2MAT but enabled with multiprecision arrays support.
%
%   See also MAT2CELL, NUM2CELL

   [varargout{1:nargout}] = mpcell2mat(varargin{:}); 
end
