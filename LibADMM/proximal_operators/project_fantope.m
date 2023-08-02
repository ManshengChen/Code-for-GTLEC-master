function X = project_fantope(Q,k)


[U,D] = eig(Q);
Dr = cappedsimplexprojection(diag(D),k);
% Dr = cappedsimplexprojection_matlab(diag(D),k);
X = U*diag(Dr)*U';