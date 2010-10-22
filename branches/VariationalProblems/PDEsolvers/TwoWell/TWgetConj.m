function p = getConj(p)

p.problem.conjNonLinearFunc = @conjNonLinearFunc;
p.problem.conjNonLinearFuncDer = @conjNonLinearFuncDer;
p.problem.conjNonLinearFuncSecDer = @conjNonLinearFuncSecDer;

%% W***(x)
function z = conjNonLinearFunc(x,y,curElem,lvl,p)
% change Basis: E = x*(1,0) + y*(0,1) = E1 * F1 + E2 * F2
% with F1 = 1/sqrt(13) (3,2) and F2 = 1/sqrt(13) (-2,3)

F2 = p.problem.F2;   %first basis vector F2 equals given constant vector F2
F1 = [0 -1;1 0]*F2;  %F1 is a basis vector, not the constant defined in startTwoWell.m!

% E1 = E*F1
E1 = [x y]*F1;
% E2 = E*F2
E2 = [x y]*F2;

G1_0 = ones(length(x),1);
G2_0 = ones(length(x),1);

% Find G=[G1;G2], such that fun(G)=0
options=optimset('Display','off');
G1 = fsolve(@(G1)Func1(E1,E2,G1),G1_0,options);
G2 = fsolve(@(G2)Func2(E1,E2,G2),G2_0,options);

maxvalue = max(0,G1.^2+G2.^2-1);
z = G1.*E1+ G2.*E2 - maxvalue.^2 - 4*G1.^2;

index1 = find( E1 == 0 );
if z(index1) <= (3/16) * E2(index1).^2
   z(index1) = (3/16) * E2(index1).^2;
end

%z = zeros(length(x),1);

function val = Func1(E1,E2,G1)
val = 1/6*(27*E1+3*(81*E1.^2+192*(1/8*E2+1/72*(81*E2.^2+192*G1.^6+576*G1.^4+576*G1.^2+192).^(1/2)).^2 ...
     -576*(1/8*E2+1/72*(81*E2.^2+192*G1.^6+576*G1.^4+576*G1.^2+192).^(1/2)).^(4/3)+576*(1/8*E2+...
     1/72*(81*E2.^2+192*G1.^6+576*G1.^4+576*G1.^2+192).^(1/2)).^(2/3)-192).^(1/2)).^(1/3) - G1;

function val = Func2(E1,E2,G2) 
val = (1/8*E2+1/72*(81*E2.^2+1/243*(27*E1+3*(81*E1.^2+192*G2.^6-576*G2.^4+576*G2.^2-192).^(1/2)).^2+...
     4/9*(27*E1+3*(81*E1.^2+192*G2.^6-576*G2.^4+576*G2.^2-192).^(1/2)).^(4/3)+16*(27*E1+3*(81*E1.^2+...
     192*G2.^6-576*G2.^4+576*G2.^2-192).^(1/2)).^(2/3)+192).^(1/2)).^(1/3) - G2;

 
%% exact derivative of conjugative nonlinear operator (relaxed)
%% DW***(F)[G] = (W***'(F),G)
function [DiffQuo_xh_y DiffQuo_x_yh] = conjNonLinearFuncDer(x,y,curElem,lvl,p)

h = 10^(-5);
DiffQuo_xh_y = (1/h)*(conjNonLinearFunc(x+h,y,curElem,lvl,p) - conjNonLinearFunc(x,y,curElem,lvl,p));
DiffQuo_x_yh = (1/h)*(conjNonLinearFunc(x,y+h,curElem,lvl,p) - conjNonLinearFunc(x,y,curElem,lvl,p));


%% exact derivative of conjugative nonlinear operator (relaxed)
%% D2W***(F)[G,H] = W***''(F) * (G,H)
function z = conjNonLinearFuncSecDer(x,y,curElem,lvl,p)

h = 10^(-5);
[DiffQuo_xhh_y DiffQuo_xh_yh] = conjNonLinearFuncDer(x+h,y  ,curElem,lvl,p);
[DiffQuo_xh_yh DiffQuo_x_yhh] = conjNonLinearFuncDer(x  ,y+h,curElem,lvl,p);
[DiffQuo_xh_y  DiffQuo_x_yh ] = conjNonLinearFuncDer(x  ,y  ,curElem,lvl,p);


DiffQuo11 = (1/h)*( DiffQuo_xhh_y - DiffQuo_xh_y );
DiffQuo12 = (1/h)*( DiffQuo_xh_yh - DiffQuo_xh_y );
DiffQuo21 = (1/h)*( DiffQuo_xh_yh - DiffQuo_x_yh );
DiffQuo22 = (1/h)*( DiffQuo_x_yhh - DiffQuo_x_yh );

%z = ones(length(x),1);
%z = DiffQuo11.^2 + DiffQuo22.^2;
%z = reshape(z,1,1,[]);

z = [DiffQuo11 DiffQuo12 DiffQuo21 DiffQuo22];
z = reshape(z,[],4)';
z = reshape(z,2,2,[]);