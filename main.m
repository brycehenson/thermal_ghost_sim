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
constants
num_tot=1000;
%num=80;
shots=50000;
corr_frac=0.010;
tof=0.416;
qe=0.1;
%matlabs zeta is bugging out
zeta3=1.202056903159594285399738161511449990764986292;

omegax=100*2*pi;
omegay=omegax;
omegaz=omegax;
omega_bar=(omegax*omegay*omegaz)^(1/3);
omega_mean=(omegax+omegay+omegaz)/3;

tc_non_int= (hb*omega_bar*num_tot^(1/3))/(kb*zeta3^(1/3));
abar=sqrt(hb/(m*omega_bar));
tc_int=tc_non_int+(-0.73*num_tot^(-1/3)*omega_mean/omega_bar-...
    1.33*num_tot^(1/6)*ahe_scat/abar)*tc_non_int;  %pethick eq11 .14


trun=tc_int*1.5; %set the temp to something above critical

%from Theory for a Hanbury Brown Twiss experiment with a ballistically expanding cloud of cold atoms
corr_widthx=hb*tof*sqrt((omegax^2)/(kb*trun))/sqrt(m);
corr_widthy=hb*tof*sqrt((omegay^2)/(kb*trun))/sqrt(m);
corr_widthz=hb*tof*sqrt((omegaz^2)/(kb*trun))/sqrt(m);
therm_width=tof*sqrt(3/2)*sqrt(2*kb*trun/m);


mode_occ=(1/sqrt(2))*num_tot*(corr_widthx*corr_widthy*corr_widthz/(therm_width)^3);
corr.window=[corr_widthx,corr_widthy,corr_widthz]/4;

corr.redges=linspace(0,corr_widthx*4,40);
corr.xedges=linspace(-corr_widthx*5,corr_widthx*5,40);

num=round(num_tot/(1+corr_frac));
norm_chunk_size=round(4*qe*num);
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
counts_in_shot=zeros(shots,1);
parfor n=1:shots
    shot=reshape(randn(num*3,1)*therm_width,[3,num]);
    coor_chance=rand(num,1)<corr_frac;
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
fprintf('mean per shot                 %2.3e \n',total_counts/shots)
fprintf('rough pairs per shot          %2.3e \n',CountUpperTriangle(total_counts/shots))
fprintf('rough total pairs             %2.3e \n',CountUpperTriangle(total_counts/shots)*shots)

fprintf('calculating in shot correlations \n')
shotscorr=CorrRadX(corr,counts_txy);

figure(1);
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


%then calculate the normalization
approx_pairs=CountUpperTriangle(total_counts);
fprintf('rough total pairs in norm     %2.3e \n',approx_pairs)
accurate_pairs=0;
for n=1:shots
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

fprintf('sampled frac of norm pairs    %3.3e\n',normcorr.pairs/accurate_pairs)
pair_sample_factor=(normcorr.pairs/shotscorr.pairs);
fprintf('sampled norm/corr pairs       %2.3e \n',pair_sample_factor)


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
saveas(gcf,'Correlations.png')
saveas(gcf,'Correlations.pdf')
savefig('Correlations.fig')
%saveas(gcf,['Correlations',num2str(num),'.png'])
%saveas(gcf,['Correlations',num2str(num),'.pdf'])
fprintf('max g2 rad                    %4.2f \n',max(radg2))
fprintf('       x                      %4.2f \n',max(xg2))

toc;
fprintf('\n\n\n\n')