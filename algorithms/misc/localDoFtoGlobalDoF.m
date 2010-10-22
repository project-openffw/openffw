function [I,J] = localDoFtoGlobalDoF(dof1,dof2)
dof1 = dof1';
I = repmat(dof1(:),1,size(dof2,2))';
J = repmat(dof2,1,size(dof1,1))';
I = I(:);
J = J(:);