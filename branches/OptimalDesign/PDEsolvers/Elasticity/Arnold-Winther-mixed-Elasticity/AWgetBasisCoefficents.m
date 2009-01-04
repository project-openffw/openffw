function p = AWgetBasisCoefficents(p)
% computes the coefficents for the P3-basis-
% functions such that they span the discrete ansatz spaces by
% Arnold and Winther

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


%%%%%%%%%%Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grad4e = p.level(end).enum.grad4e;
normals4e = p.level(end).enum.normals4e;
nrElems = p.level(end).nrElems;
e4ed = p.level(end).enum.e4ed;
ed4e = p.level(end).enum.ed4e;
length4ed = p.level(end).enum.length4ed;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%computation of the coefficents for the basis functions of sigma such that
%the divergence of sigma is in P1

IdentityBig = eye(9);
IdentitySmall = eye(3);
IdentityTriple = [IdentitySmall,IdentitySmall,IdentitySmall];

N = zeros(2,3,3);
Q = zeros(6,3,4);
dummyZeros = zeros(2,3);

coeff1 = zeros(12,9);
coeff1([1 2],[1 2 3]) = 10;
coeff1([5 6],[4 5 6]) = 10;
coeff1([9 10],[7 8 9]) = 10;

coeff2 = zeros(12,9);
coeff2([3 4],[1 2 3])   = -1;
coeff2([7 8],[4 5 6])   = -1;
coeff2([11 12],[7 8 9]) = -1;

rhs = [eye(24);zeros(6,24)];

basisCoefficents = zeros(30,24,nrElems);
matrixC = zeros(30,30,nrElems);

for curElem = 1:nrElems
	curGrad = grad4e(:,:,curElem);
	curEdges = ed4e(curElem,:);
	curNormals = normals4e(:,:,curElem);
	curEdgeLengths = length4ed(curEdges);
	for j = 1:3
		N(:,:,j) = [curNormals(j,1),    0,          curNormals(j,2);
			0,              curNormals(j,2),curNormals(j,1)];
	end

	Nsig = zeros(2,3,3);
	signum = ones(1,3);
	signum(find(curElem == e4ed(curEdges,2))) = -1;
	Nsig(:,:,1) = N(:,:,1)*signum(1);
	Nsig(:,:,2) = N(:,:,2)*signum(2);
	Nsig(:,:,3) = N(:,:,3)*signum(3);

	Q(:,:,1) = [3*(curGrad(1,1)-curGrad(2,1))   , 0, 3*(curGrad(1,2)-curGrad(2,2));
		curGrad(1,1)                    , 0,    curGrad(1,2);
		-curGrad(2,1)                   , 0,    -curGrad(2,2);
		0, 3*(curGrad(1,2)-curGrad(2,2)), 3*(curGrad(1,1)-curGrad(2,1));
		0, curGrad(1,2)                 , curGrad(1,1);
		0, -curGrad(2,2)                ,-curGrad(2,1)  ];

	Q(:,:,2) = [-curGrad(3,1)                   , 0,    -curGrad(3,2);
		3*(curGrad(2,1)-curGrad(3,1))   , 0, 3*(curGrad(2,2)-curGrad(3,2));
		curGrad(2,1)                    , 0,    curGrad(2,2);
		0, -curGrad(3,2)                ,-curGrad(3,1);
		0, 3*(curGrad(2,2)-curGrad(3,2)), 3*(curGrad(2,1)-curGrad(3,1));
		0, curGrad(2,2)                 , curGrad(2,1)  ];

	Q(:,:,3) = [curGrad(3,1)                    , 0,    curGrad(3,2);
		-curGrad(1,1)                   , 0,   -curGrad(1,2);
		3*(curGrad(3,1)-curGrad(1,1))   , 0, 3*(curGrad(3,2)-curGrad(1,2));
		0, curGrad(3,2)                 , curGrad(3,1);
		0, -curGrad(1,2)                ,-curGrad(1,1);
		0, 3*(curGrad(3,2)-curGrad(1,2)), 3*(curGrad(3,1)-curGrad(1,1)) ];

	Q(:,:,4) = [curGrad(3,1)                   , 0,    curGrad(3,2);
		curGrad(1,1)                   , 0,    curGrad(1,2);
		curGrad(2,1)                   , 0,    curGrad(2,2);
		0, curGrad(3,2)                , curGrad(3,1);
		0, curGrad(1,2)                , curGrad(1,1);
		0, curGrad(2,2)                , curGrad(2,1)   ];

	Q = 60*Q;
	Q = max(curEdgeLengths) * Q;

	R = [30*Nsig(:,:,1),   30*Nsig(:,:,1),    dummyZeros;
		-5*N(:,:,1),    5*N(:,:,1),    dummyZeros;
		dummyZeros ,   30*Nsig(:,:,2),    30*Nsig(:,:,2);
		dummyZeros ,   -5*N(:,:,2),     5*N(:,:,2);
		30*Nsig(:,:,3),    dummyZeros,    30*Nsig(:,:,3);
		5*N(:,:,3),    dummyZeros,    -5*N(:,:,3)];

	S = [     Nsig(:,:,1),     dummyZeros,    dummyZeros;
		N(:,:,1),     dummyZeros,    dummyZeros;
		dummyZeros,       Nsig(:,:,2),    dummyZeros;
		dummyZeros,       N(:,:,2),    dummyZeros;
		dummyZeros,     dummyZeros,       Nsig(:,:,3);
		dummyZeros,     dummyZeros,       N(:,:,3)];

	S1 = S.*coeff1;
	S2 = S.*coeff2;

	C = [60*IdentityBig,    zeros(9,21);
		R,                         S1,          S2,         zeros(12,3);
		20*IdentityTriple, 5*IdentityTriple, zeros(3,9),  IdentitySmall;
		zeros(6,18),         Q(:,:,1),Q(:,:,2),Q(:,:,3),       Q(:,:,4)];

	C = 1/60*C;
	matrixC(:,:,curElem) = C;

	basisCoefficents(:,:,curElem) = C\rhs;
end

%%%%%%%%Output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).basisCoefficents = basisCoefficents;
p.level(end).matrixC = matrixC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
