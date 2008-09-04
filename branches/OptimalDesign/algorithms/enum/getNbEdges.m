function NbEd  = getNbEdges(Nb,ed4n)

if ~isempty(Nb)
    NbEd = rowaddr(ed4n,Nb(:,1),Nb(:,2));
else
    NbEd = [];
end