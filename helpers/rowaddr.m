function result = rowaddr(A,I,J,S)
% rowaddr.m pointwise reads entries [I,J] from A if S is omitted, i.e.,
% result = A( I(j,k), J(j,k)).
% If S is submitted, then entries [I,J] from A are pointwise set to S and 
% the new A is returned as a result.

% Copyright 2007 Jan Reininghaus, David Guenther, Andreas Byfut
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

if nargin == 3
    % initialise result
    result = zeros(size(I));
    % if size of A is small enough we can use linear adressing
    if(numel(A) < 2^31)
        for curCol = 1:size(I,2)
            linIndex = sub2ind(size(A),I(:,curCol),J(:,curCol));
            result(:,curCol) = full( A(linIndex) );
        end
    else
        for j = 1:size(I,2)
            dummyI = I(:,j);
            dummyJ = J(:,j);
            dummyres = zeros(size(I,1),1);
            for k = 1:size(I,1)
                dummyres(k) = A(dummyI(k),dummyJ(k));
            end
            result(:,j) = dummyres;
        end
    end
elseif nargin == 4
    % initialise result
    result = A;
    % if size of A is small enough we can use linear adressing
    if(numel(A) < 2^31)
        for curCol = 1:size(I,2)
            linIndex = sub2ind(size(A),I(:,curCol),J(:,curCol));
            result(linIndex) = S(:,curCol);
        end
    else
        for j = 1:size(I,2)
            dummyI = I(:,j);
            dummyJ = J(:,j);
            dummyS = S(:,j);
            for k = 1:size(I,1)
                result(dummyI(k),dummyJ(k)) = dummyS(k);
            end
        end
    end
else
    if nargin < 3  
        error('Not enough parameters submitted,');
    end
end
