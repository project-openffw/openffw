function p = redGreenBlue(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Db = p.level(end).geom.Db;
Nb = p.level(end).geom.Nb;

DbEd = p.level(end).enum.DbEd;
NbEd = p.level(end).enum.NbEd;
ed4e = p.level(end).enum.ed4e;
midPoint4ed = p.level(end).enum.midPoint4ed;

% refineElemsBisec5 = p.level(end).refineElemsBisec5;
refineEdges = p.level(end).refineEdges;
nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
nextLevel = length(p.level) + 1;
if(nnz(refineEdges) == 0)
% 	p.level(nextLevel).geom = p.level(nextLevel-1).geom;
% 	p.level(nextLevel-1).enum.newNode4ed = [];
    fprintf('No refinement!!! \n');
	return;
end


% store the parent element number for each new element
% upper bound for new nr. of elements (exact for bisec5 refinement)
% parents4e = zeros(6*nrElems,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update C4N									   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create new node numbers from refineEdges
newNode4ed = zeros(1,nrEdges);
newNode4ed( find(refineEdges) ) = (nrNodes+1):(nrNodes+nnz(refineEdges));

% Create coordinates of the new nodes
[dontUse,J,S] = find(newNode4ed);
c4n(S,:) = midPoint4ed(J,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update N4E									   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newNode4e = newNode4ed(ed4e);
unrefinedElems = find( all(newNode4e == 0 ,2) );
refineElems = find( any(newNode4e,2) );
nrMarkedEd4MarkedElems = sum(refineEdges( ed4e(refineElems,:) ),2);


%%%%%%%%%%%%%%
% old element number for unrefined elements
nrUnrefinedElems = length(unrefinedElems);
parents4e(1:nrUnrefinedElems) = unrefinedElems;
%%%%%%%%%%%%%%%

newn4e = n4e(unrefinedElems,:);
% Green refinement
nrGreenElems = 0;
I = find(nrMarkedEd4MarkedElems == 1);
if ~isempty(I)
	gElems = refineElems(I);
	% todo
	[dontUse,dontUse,newN] = find( newNode4e(gElems,:)' );
	newGreenElems = [n4e(gElems,[2 3]) newN;...
                     n4e(gElems,[3 1]) newN];
    nrGreenElems = size(newGreenElems,1);
    dummy = [gElems;gElems];
    parents4e( size(newn4e,1)+1 : size(newn4e,1)+nrGreenElems ) = dummy;
  	newn4e = [newn4e;newGreenElems];
end

% Blue refinement
nrBlueLeftElems = 0;
nrBlueRightElems = 0;

I = find(nrMarkedEd4MarkedElems == 2);
if ~isempty(I)
	bElems = refineElems(I);
	[RefineEdge,dontUse,newN] = find( newNode4e(bElems,:)' );
	RefineEdge = reshape(RefineEdge,2,length(I));
	newN = reshape(newN,2,length(I))';
	% BlueLeft refinement
	bl = find(RefineEdge(2,:) == 3);
	if ~isempty(bl)
		newBlueLeft = [ n4e(bElems(bl),1), newN(bl,1), newN(bl,2);...
                        n4e(bElems(bl),2), n4e(bElems(bl),3), newN(bl,1);...
                        newN(bl,1), n4e(bElems(bl),3),  newN(bl,2)];
        nrBlueLeftElems = size(newBlueLeft,1);
        dummy = [bElems(bl);bElems(bl);bElems(bl)];
        parents4e( size(newn4e,1)+1 : size(newn4e,1)+nrBlueLeftElems ) = dummy;
        newn4e = [newn4e; newBlueLeft]; 
	end
	% BlueRight refinement
	br = find(RefineEdge(2,:) == 2);
	if ~isempty(br)
		newBlueRight = [ newN(br,1), n4e(bElems(br),2), newN(br,2);...
                        n4e(bElems(br),3), newN(br,1), newN(br,2);...
                        n4e(bElems(br),3), n4e(bElems(br),1), newN(br,1)];
        nrBlueRightElems = size(newBlueRight,1);
        dummy = [bElems(br);bElems(br);bElems(br)];
        parents4e( size(newn4e,1)+1 : size(newn4e,1)+nrBlueRightElems ) = dummy;
        newn4e = [newn4e; newBlueRight];
	end
end

% Red refinement (Do not refine elements refine for bisec5 refinement)
% nrMarkedEd4MarkedElems(refineElemsBisec5(refineElems)) = 0;
nrRedElems = 0;
I = find(nrMarkedEd4MarkedElems == 3);
if ~isempty(I)
	rElems = refineElems(I);
	[dontUse,dontUse,newN] = find( newNode4e(rElems,:)' );
	newN = reshape(newN,3,length(I))';
	newRedElems = [newN(:,[2 3 1]);...
		n4e(rElems,1) newN(:,1) newN(:,3);...
		newN(:,1) n4e(rElems,2) newN(:,2);...
		newN(:,3) newN(:,2) n4e(rElems,3)];
    nrRedElems = size(newRedElems,1);
    dummy = [rElems;rElems;rElems;rElems];
    parents4e( size(newn4e,1)+1 : size(newn4e,1)+nrRedElems ) = dummy;
    newn4e = [newn4e; newRedElems];
end

% Bisec5 / Oscillation reduction
% bisec5Elems = find(refineElemsBisec5);
% if ~isempty(bisec5Elems)
% 	newNode4elem = zeros(1,nrElems);
% 	newNrNodes = size(c4n,1);
% 	newNode4elem(bisec5Elems) = (newNrNodes+1):(newNrNodes+nnz(refineElemsBisec5));
% 
% 	% Create coordinates of the new nodes
% 	[I,J,S] = find(newNode4elem);
% 	nonRefEdges = ed4e(bisec5Elems,[2 3]);
% 	midPoints1 = midPoint4ed(nonRefEdges(:,1),:);
% 	midPoints2 = midPoint4ed(nonRefEdges(:,2),:);
% 	c4n(S,:) = (midPoints1 + midPoints2)/2;
% 	
% 	edgeNodes = newNode4e(bisec5Elems,:);
% 	[dontUse,dontUse,innerNode] = find( newNode4elem(bisec5Elems)' );
% 	
% 	newOscElems = [n4e(bisec5Elems,1),edgeNodes(:,1),edgeNodes(:,3);...
% 					edgeNodes(:,1),n4e(bisec5Elems,2),edgeNodes(:,2);...
% 					edgeNodes(:,1),edgeNodes(:,2),innerNode;...
% 					edgeNodes(:,3),edgeNodes(:,1),innerNode;...
% 					n4e(bisec5Elems,3),edgeNodes(:,3),innerNode;...
% 					edgeNodes(:,2),n4e(bisec5Elems,3),innerNode];
% 	newn4e = [newn4e; newOscElems];
% 	nrOscElems = size(newOscElems,1);
% 	dummy = repmat(bisec5Elems,6,1);
% 	parents4e( size(newn4e,1)+1 : size(newn4e,1)+nrOscElems ) = dummy;
% end

% reduce the size of parents4e to the new number of elements
% parents4e = parents4e(1:size(newn4e,1));
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update boundary								   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Db = updateBoundary(Db,DbEd,newNode4ed); 
Nb = updateBoundary(Nb,NbEd,newNode4ed);

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(nextLevel-1).enum.newNode4ed = newNode4ed;
p.level(nextLevel).enum.parents4e = parents4e;
p.level(nextLevel).geom.n4e = newn4e;
p.level(nextLevel).geom.c4n = c4n;
p.level(nextLevel).geom.Db = Db;
p.level(nextLevel).geom.Nb = Nb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newBoundary = updateBoundary(oldB,ed4b,newNode4ed)
if(isempty(oldB))
	newBoundary = [];
else
	unrefinedEd = find(~newNode4ed(ed4b));
	refineEd = find(newNode4ed(ed4b));
	newBoundary = [oldB(unrefinedEd,:);...
		oldB(refineEd,1) newNode4ed(ed4b(refineEd))' ;...
		newNode4ed(ed4b(refineEd))' oldB(refineEd,2)];
end
