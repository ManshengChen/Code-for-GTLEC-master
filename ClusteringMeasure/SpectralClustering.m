
function groups = SpectralClustering(CKSym,n)

warning off;
N = size(CKSym,1);
MAXiter = 1000; % Maximum number of iterations for KMeans 
REPlic = 20; % Number of replications for KMeans



DN = diag( 1./sqrt(sum(CKSym)+eps) );
LapN = speye(N) - DN * CKSym * DN;
[uN,sN,vN] = svd(LapN);
kerN = vN(:,N-n+1:N);
for i = 1:N
    kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
end
groups = kmeans(kerNS,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');