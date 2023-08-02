function result=clustering8(S, cls_num, gt)
[C] = SpectralClustering(S,cls_num);
result = Clustering8Measure(C,gt);
end
