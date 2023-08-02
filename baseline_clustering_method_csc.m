function [AC1,MIhat1]=baseline_clustering_method_csc(W,nClass,gndnew)


% returen the results of SC
% reture the results of k-means

D=sum(W,2);
DD=D.^(-0.5);
L=diag(DD)*W*diag(DD);
n=length(W);

[V,~] = eig(L);
% H = V(:,1421:end);
% H = V(:,end-k+1:end);
H = V(:,1:nClass);
Q=normr(H);
% Q=H;

labelnew = litekmeans(Q,nClass,'Replicates',20);
    %                 gndnew=gnd;
MIhat1 = MutualInfo(gndnew,labelnew);
labelnew = bestMap(gndnew,labelnew);
AC1 = length(find(gndnew == labelnew))/length(gndnew);
disp('sc finished')

