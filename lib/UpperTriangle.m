function pairs=UpperTriangle(n)
%upper tirangle pairs (above the diagonal)
pairs=zeros(2,CountUpperTriangle(n));
count=1;
for i=1:n
    for  j=1+i:n
        pairs(:,count)=[i;j];
        count=count+1;
    end
end
end