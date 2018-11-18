%numerical ghost imaging
%simulate thermal ghost imaging system
%create a thermal distributed cloud and then add in correlations to
%measure the correlation function
%then split the cloud into two and try to do ghost imaging

%to measure the correlations more quickly and save memory overhead will bin on the fly
%myCluster = parcluster('local')
%myCluster.NumWorkers = 50
%saveProfile(myCluster)
%parpool('local',50)

num_list= round(logspace(log10(10),log10(500),50));
pop_study.num=zeros(size(num_list,2),1);
pop_study.maxg2x=zeros(size(num_list,2),1);
pop_study.maxg2rad=zeros(size(num_list,2),1);
figure(1);
for q=1:size(num_list,2)
tic;
num=num_list(q);
%num=80;
shots=3000;
corr_frac=1;
corr_width=0.19;
qe=1;
corr.redges=linspace(0,1,40);
corr.xedges=linspace(-1,1,40);
corr.window=[0.1,0.1,0.1];
norm_chunk_size=50*qe*num*(1+corr_frac);

%generate some correlated data
fprintf('generating test data \n')
counts_txy={};
counts_in_shot=zeros(shots,1);
for n=1:shots
    shot=reshape(randn(num*3,1),[3,num]);
    coor_chance=rand(num,1)<corr_frac;
    if sum(coor_chance)>0
        delpos=reshape(randn(sum(coor_chance)*3,1),[3,sum(coor_chance)])*corr_width;
        shot_with_corr=[shot,shot(:,coor_chance)+delpos];
    else
        shot_with_corr=shot;
    end
    
    det_chance=rand(size(shot_with_corr,2),1)<qe;
    num_det=sum(det_chance);
    if  num_det>1
        det_counts=shot_with_corr(:, det_chance);
        %sort in time to make sure everything is robust to that
        %[values, order] = sort(det_counts(1,:));
        %det_counts=det_counts(:,order);
        counts_txy=[counts_txy,{det_counts}];
    else
        counts_txy=[counts_txy,{[]}];
    end
    counts_in_shot(n)=num_det;
end
total_counts=size([counts_txy{:}],2);
fprintf('total counts over all shots %2.4e \n',total_counts)
fprintf('mean per shot %2.4e \n',total_counts/shots)
fprintf('rough pairs per shot %2.4e \n',CountUpperTriangle(total_counts/shots))
fprintf('rough total pairs %2.4e \n',CountUpperTriangle(total_counts/shots)*shots)

fprintf('calculating in shot correlations \n')
shotscorr=CorrRadX(corr,counts_txy);

%then calculate the normalization
approx_pairs=CountUpperTriangle(total_counts);
fprintf('rough total pairs in norm %2.4e \n',approx_pairs)
accurate_pairs=0;
for n=1:shots
    accurate_pairs=accurate_pairs+(total_counts-counts_in_shot(n))*counts_in_shot(n)*0.5;
end
fprintf('actual total pairs in norm %2.4e \n',accurate_pairs)
fprintf('ratio acc./approx %3.5f \n',accurate_pairs/approx_pairs)


all_counts=[counts_txy{:}];
all_counts=all_counts(:,randperm(total_counts));
n=1;
i=1;
counts_chunked={};
while n<total_counts
    max_index=min([total_counts,n+norm_chunk_size-1]);
    counts_chunked{i}=all_counts(:,n:max_index);
    i=i+1;
    n=n+norm_chunk_size+1;
end

fprintf('calculating in inter-shot correlations \n ')
normcorr=CorrRadX(corr,counts_chunked);

fprintf('sampled %3.3e fraction of norm pairs \n',normcorr.pairs/accurate_pairs)
pair_sample_factor=(normcorr.pairs/shotscorr.pairs);
fprintf('ratio of sampled norm pairs to corr pairs %2.3e \n',pair_sample_factor)



clf
set(gcf,'color','w');
subplot(2,3,1)
plot(shotscorr.rad_centers,shotscorr.rad_bins.*shotscorr.rad_centers.^-2/shotscorr.pairs)
title('In Shot Rad Dist')
xlabel('Radial Seperation')
subplot(2,3,4)
plot(shotscorr.x_centers,shotscorr.x_bins/shotscorr.pairs)
title('In Shot X Dist (windowed)')
xlabel('X Seperation')
subplot(2,3,2)
plot(normcorr.rad_centers,normcorr.rad_bins.*normcorr.rad_centers.^-2/normcorr.pairs)
title('Between Shot Rad Dist')
subplot(2,3,5)
plot(normcorr.x_centers,normcorr.x_bins/normcorr.pairs)
title('Between Shot X Dist (windowed)')
subplot(2,3,6)
xg2=pair_sample_factor*shotscorr.x_bins./normcorr.x_bins;
plot(shotscorr.x_centers,xg2)
title('g2 X')
subplot(2,3,3)
radg2=pair_sample_factor*shotscorr.rad_bins./normcorr.rad_bins;
plot(shotscorr.rad_centers,radg2)
title('g2 Rad')
pause(0.1);
saveas(gcf,['Correlations',num2str(num),'.png'])
saveas(gcf,['Correlations',num2str(num),'.pdf'])
fprintf('max g2 rad %4.2f x %4.2f \n',[max(xg2),max(radg2)])


pop_study.num(q)=num;
pop_study.maxg2x(q)=max(xg2);
pop_study.maxg2rad(q)=max(radg2);
save('popstudy.mat','pop_study')

toc;
end

mode_occ=(corr_width)^3.*pop_study.num;
figure(2)
clf
loglog(mode_occ,pop_study.maxg2rad-1,'ob',...
    mode_occ,pop_study.maxg2x-1,'xk',...
    mode_occ,1./mode_occ,'-r')
axis tight
xlabel('Mode Occupancy')
ylabel('G2 amp')
saveas(gcf,'PopStudy.png')
saveas(gcf,'PopStudy.pdf')