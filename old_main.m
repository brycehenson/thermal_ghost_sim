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
tic;

shot_counts=100;
%num=80;
number_shots=5000;
corr_frac=1;
tof=0.416;
qe=1;
omegax=100*2*pi;
omegay=omegax;
omegaz=omegax;
trun=tc_int*1.5; %set the temp to something abovecritical



%% Setting Up The Enviorment
%add all subfolders to the path
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));


%matlabs zeta is bugging out
zeta3=zeta(3);

omega_bar=(omegax*omegay*omegaz)^(1/3);
omega_mean=(omegax+omegay+omegaz)/3;

tc_non_int= (const.hb*omega_bar*shot_counts^(1/3))/(const.kb*zeta3^(1/3));
abar=sqrt(const.hb/(const.mhe*omega_bar));
tc_int=tc_non_int+(-0.73*shot_counts^(-1/3)*omega_mean/omega_bar-...
    1.33*shot_counts^(1/6)*const.ahe_scat/abar)*tc_non_int;  %pethick eq11 .14



%call the constants function that makes some globals
hebec_constants 


%from Theory for a Hanbury Brown Twiss experiment with a ballistically expanding cloud of cold atoms
corr_widthx=const.hb*tof*sqrt((omegax^2)/(const.kb*trun))/sqrt(const.mhe);
corr_widthy=const.hb*tof*sqrt((omegay^2)/(const.kb*trun))/sqrt(const.mhe);
corr_widthz=const.hb*tof*sqrt((omegaz^2)/(const.kb*trun))/sqrt(const.mhe);
therm_width=tof*sqrt(3/2)*sqrt(2*const.kb*trun/const.mhe);


mode_occ=(1/sqrt(2))*shot_counts*(corr_widthx*corr_widthy*corr_widthz/(therm_width)^3);
corr_opts=[]
corr_opts.window=[corr_widthx,corr_widthy,corr_widthz]/4;

corr_opts.redges=linspace(0,corr_widthx*4,40);
corr_opts.xedges=linspace(-corr_widthx*5,corr_widthx*5,40);

num_corrected=round(shot_counts/(1+corr_frac));
norm_chunk_size=round(4*qe*num_corrected);
fprintf('critical temp                 %2.3e \n',tc_int)
fprintf('thermal width                 %2.3e \n',therm_width)
fprintf('corr widthx                   %2.3e \n',corr_widthx)
fprintf('corr widthy                   %2.3e \n',corr_widthy)
fprintf('corr widthz                   %2.3e \n',corr_widthz)
fprintf('pixels x                      %2.3f \n',therm_width/corr_widthx)
fprintf('total pixels xy               %2.3e \n',(therm_width/corr_widthx)*(therm_width/corr_widthy))
fprintf('mode occupancy                %2.3e\n',mode_occ)
fprintf('predicted corr amp            %2.3e \n',1+1/mode_occ)
%generate some correlated data
fprintf('generating test data \n')
counts_txy={};
counts_in_shot=zeros(number_shots,1);
parfor n=1:number_shots
    shot=reshape(randn(num_corrected*3,1)*therm_width,[3,num_corrected]);
    coor_chance=rand(num_corrected,1)<corr_frac;
    if sum(coor_chance)>0
        delpos=reshape(randn(sum(coor_chance)*3,1),[3,sum(coor_chance)]).*repmat([corr_widthz;corr_widthx;corr_widthy],[1,sum(coor_chance)]);
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
        counts_txy{n}=det_counts;
    else
        counts_txy{n}=[];
    end
    counts_in_shot(n)=num_det;
end
total_counts=size([counts_txy{:}],2);
fprintf('total counts over all shots   %2.3e \n',total_counts)
fprintf('mean per shot                 %2.3e \n',total_counts/number_shots)
fprintf('rough pairs per shot          %2.3e \n',CountUpperTriangle(total_counts/number_shots))
fprintf('rough total pairs             %2.3e \n',CountUpperTriangle(total_counts/number_shots)*number_shots)

%%


fprintf('calculating in shot correlations \n')
shotscorr=CorrRadX(corr,counts_txy);

figure(1);
clf
set(gcf,'color','w');
subplot(2,3,1)
plot(shotscorr_opts.rad_centers,shotscorr_opts.rad_bins.*shotscorr_opts.rad_centers.^-2/shotscorr_opts.pairs)
title('In Shot Rad Dist')
xlabel('Radial Seperation')
subplot(2,3,4)
plot(shotscorr_opts.x_centers,shotscorr_opts.x_bins/shotscorr_opts.pairs)
title('In Shot X Dist (windowed)')
xlabel('X Seperation')


%then calculate the normalization
approx_pairs=CountUpperTriangle(total_counts);
fprintf('rough total pairs in norm     %2.3e \n',approx_pairs)
accurate_pairs=0;
for n=1:number_shots
    accurate_pairs=accurate_pairs+(total_counts-counts_in_shot(n))*counts_in_shot(n)*0.5;
end
fprintf('actual total pairs in norm    %2.3e \n',accurate_pairs)
fprintf('ratio acc./approx             %2.3e \n',accurate_pairs/approx_pairs)


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

fprintf('calculating inter-shot correlations \n ')
normcorr=CorrRadX(corr,counts_chunked);

fprintf('sampled frac of norm pairs    %3.3e\n',normcorr_opts.pairs/accurate_pairs)
pair_sample_factor=(normcorr_opts.pairs/shotscorr_opts.pairs);
fprintf('sampled norm/corr pairs       %2.3e \n',pair_sample_factor)


subplot(2,3,2)
plot(normcorr_opts.rad_centers,normcorr_opts.rad_bins.*normcorr_opts.rad_centers.^-2/normcorr_opts.pairs)
title('Between Shot Rad Dist')
subplot(2,3,5)
plot(normcorr_opts.x_centers,normcorr_opts.x_bins/normcorr_opts.pairs)
title('Between Shot X Dist (windowed)')
subplot(2,3,6)
xg2=pair_sample_factor*shotscorr_opts.x_bins./normcorr_opts.x_bins;
plot(shotscorr_opts.x_centers,xg2)
title('g2 X')
subplot(2,3,3)
radg2=pair_sample_factor*shotscorr_opts.rad_bins./normcorr_opts.rad_bins;
plot(shotscorr_opts.rad_centers,radg2)
title('g2 Rad')
pause(0.1);
saveas(gcf,'Correlations.png')
saveas(gcf,'Correlations.pdf')
savefig('Correlations.fig')
%saveas(gcf,['Correlations',num2str(num),'.png'])
%saveas(gcf,['Correlations',num2str(num),'.pdf'])
fprintf('max g2 rad                    %4.2f \n',max(radg2))
fprintf('       x                      %4.2f \n',max(xg2))

toc;
fprintf('\n\n\n\n')