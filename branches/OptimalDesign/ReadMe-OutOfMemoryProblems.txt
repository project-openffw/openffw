Optimal Design:

This package solves a nonlinear optimal design problem. 
In order to do so the Optimization Toolbox of MATLAB is needed. 
Especially the functions 'fsolve.m' and 'nlsq.m' are needed.

Since the nonlinear least square solver nlsq.m of MATLAB does
not use sparse matrices the computation will terminate at 5k DoFs
with the error message "Out of memory."

To fix this problem one needs to modify nlsq.m manually. The command
eye(...) that generates a full unit matrix of specific dimensions needs
to be replaced with speye(...). Now, the generated unit matrix is sparse.
(There are four commands that needs to be replaced).

With this modification calculations up to 100k DoFs are possible
without any problems.

Details to Optimal Design and the implementation can be found in the diploma
thesis (german).

----
David Guenther
