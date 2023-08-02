function [accuracy, true_labels, CM] = calculateAccuracy(yte, y)

N = length(y);

cluster_names = unique(yte);
accuracy = 0;
maxInd = 1;

perm = perms(unique(y));
[pN pM] = size(perm);

true_labels = y;

for i=1:pN
    flipped_labels = zeros(1,N);
    for cl = 1 : pM
        flipped_labels(yte==cluster_names(cl)) = perm(i,cl);
    end

    testAcc = sum(flipped_labels == y')/N;
    if testAcc > accuracy
        accuracy = testAcc;
        maxInd = i;
        true_labels = flipped_labels;
    end

end

CM = zeros(pM,pM);
for rc = 1 : pM
    for cc = 1 : pM
        CM(rc,cc) = sum( ((y'==rc) .* (true_labels==cc)) );
    end
end