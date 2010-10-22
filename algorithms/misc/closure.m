function p = closure(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ed4e = p.level(end).enum.ed4e;
refineEdges = p.level(end).refineEdges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

refineEdges = refineEdges;
refEd4e = ed4e(:,1);

I =  refineEdges(ed4e(:,2)) | refineEdges(ed4e(:,3));
while nnz(refineEdges(refEd4e(I))) < nnz(I);
   refineEdges(refEd4e(I)) = true;
   I =  refineEdges(ed4e(:,2)) | refineEdges(ed4e(:,3));
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).refineEdges = refineEdges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
