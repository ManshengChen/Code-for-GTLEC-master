function [newL2] = bestMap(L1,L2)
 

L1 = L1(:);
L2 = L2(:);
if size(L1) ~= size(L2)
    error('size(L1) must == size(L2)');
end

Label1 = unique(L1);
nClass1 = length(Label1);
Label2 = unique(L2);
nClass2 = length(Label2);

nClass = max(nClass1,nClass2);
G = zeros(nClass);
for i=1:nClass1
	for j=1:nClass2
		G(i,j) = length(find(L1 == Label1(i) & L2 == Label2(j)));
	end
end

[c,t] = hungarian(-G);
newL2 = zeros(size(L2));
for i=1:nClass2
    newL2(L2 == Label2(i)) = Label1(c(i));
end
