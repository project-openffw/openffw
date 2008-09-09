function p = closure(p)
% closure algorithm
% Marks further edges such that if one edge of an element
% T is marked also it's reference edge E(T) is marked.

% Copyright 2007 Jan Reininghaus, David Guenther
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ed4e = p.level(end).enum.ed4e;
markedEdges = p.level(end).markedEdges;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

markedEdgesBeforeClosure = markedEdges;
refEd4e = ed4e(:,1);

I =  markedEdges(ed4e(:,2)) | markedEdges(ed4e(:,3));
while nnz(markedEdges(refEd4e(I))) < nnz(I);
   markedEdges(refEd4e(I)) = true;
   I =  markedEdges(ed4e(:,2)) | markedEdges(ed4e(:,3));
end


%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.level(end).markedEdges = markedEdges;
p.level(end).markedEdgesBeforeClosure = markedEdgesBeforeClosure;
