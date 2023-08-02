function [Acc,rand_index,match]=AccMeasure(T,idx)


k=max(T);
n=length(T);
for i=1:k
    temp=find(T==i);
    a{i}=temp; %#ok<AGROW>
end

b1=[];
t1=zeros(1,k);
for i=1:k
    tt1=find(idx==i);
    for j=1:k
       t1(j)=sum(ismember(tt1,a{j}));
    end
    b1=[b1;t1]; %#ok<AGROW>
end
    Members=zeros(1,k); 
    
P = perms((1:k));
    Acc1=0;
for pi=1:size(P,1)
    for ki=1:k
        Members(ki)=b1(P(pi,ki),ki);
    end
    if sum(Members)>Acc1
        match=P(pi,:);
        Acc1=sum(Members);
    end
end

rand_ss1=0;
rand_dd1=0;
for xi=1:n-1
    for xj=xi+1:n
        rand_ss1=rand_ss1+((idx(xi)==idx(xj))&&(T(xi)==T(xj)));
        rand_dd1=rand_dd1+((idx(xi)~=idx(xj))&&(T(xi)~=T(xj)));
    end
end
rand_index=200*(rand_ss1+rand_dd1)/(n*(n-1));
Acc=Acc1/n*100; 
match=[1:k;match];