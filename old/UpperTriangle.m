function pairs=UpperTriangle(n)
%upper tirangle pairs (above the diagonal)
pairs=zeros(CountUpperTriangle(n),2);
count=1;
for i=1:n
    for  j=1+i:n
        pairs(count,:)=[i;j];
        count=count+1;
    end
end
end