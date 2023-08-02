function [CA]=compute_Av(BP)
% input the base partition, where each column BP is a base partitions
n= size(BP,1); % number of samples
m=size(BP,2); % number of the BPs
% convert the BP into a feature vector
CA=cell(m,1);
for i=1:m
    v=BP(:,i);
    s=zeros(n,max(v));
    for j=1:n
        for k=1:max(v)
            if v(j)==k
                s(j,k)=1;
            end
        end
    end
    CA{i}=s*s';
end
