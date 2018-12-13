function out=CorrRadX(corr,counts_txy)
%pass in a cell array of counts and then calculate the correlations using a
%bin on the fly method
%for the the gxyz we require the difference in the other axes to be less
%than a threshold

updates=100;

%find un-norm g2
shots =size(counts_txy,2);

update_frac=shots/updates;
parfor_progress_imp(updates);



[rad_bins,corr.redges]=histcounts([],corr.redges);
[x_bins,corr.xedges]=histcounts([],corr.xedges);
x_bins=zeros(shots,size(x_bins,2));
rad_bins=zeros(shots,size(rad_bins,2));
pairs_count=zeros(1,shots);
for n=1:shots
    shot_txy=counts_txy{n};
    pairs=UpperTriangle(size(shot_txy,2));
    pairs_count(n)=size(pairs,2);
    %fprintf('shot %3i \n',n)
    delta=zeros(size(pairs,2),3);
    for i=1:size(pairs,2)
        delta(i,:)=shot_txy(pairs(1,i),:)-shot_txy(pairs(2,i),:);
    end
    xmask=abs(delta(:,1))<corr.window(1) & abs(delta(:,3))<corr.window(3);
    x_bins(n,:)=histcounts(delta(2,xmask),corr.xedges);
    rad_bins(n,:)=histcounts(sqrt(sum(delta.^2,1)),corr.redges);
    if rand<(1/update_frac)
        parfor_progress_imp;
    end
end
parfor_progress_imp(0);

out.rad_centers=(corr.redges(2:end)+corr.redges(1:end-1))/2;
out.x_centers=(corr.xedges(2:end)+corr.xedges(1:end-1))/2;
out.pairs=sum(pairs_count);
out.x_bins=sum(x_bins,1);
out.rad_bins=sum(rad_bins,1);


end