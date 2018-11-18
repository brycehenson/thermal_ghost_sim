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
%num_tot=100;
%shots=30000;
%corr_frac=0.01;
tof=0.416;
num_tot=100;
shots=3000;
corr_frac=1;
qe=1;

omegax=100*2*pi;
omegay=omegax;
omegaz=omegax;
omega_bar=(omegax*omegay*omegaz)^(1/3);
omega_mean=(omegax+omegay+omegaz)/3;

tc_non_int= (hb*omega_bar*num_tot^(1/3))/(kb*zeta(3)^(1/3));
abar=sqrt(hb/(m*omega_bar));
tc_int=tc_non_int+(-0.73*num_tot^(-1/3)*omega_mean/omega_bar-...
    1.33*num_tot^(1/6)*ahe_scat/abar)*tc_non_int;  %pethick eq11 .14


trun=tc_int*2;

%from Theory for a Hanbury Brown Twiss experiment with a ballistically expanding cloud of cold atoms
corr_widthx=hb*tof*sqrt((omegax^2)/(kb*trun))/sqrt(m);
corr_widthy=hb*tof*sqrt((omegay^2)/(kb*trun))/sqrt(m);
corr_widthz=hb*tof*sqrt((omegaz^2)/(kb*trun))/sqrt(m);
therm_width=tof*sqrt(3/2)*sqrt(2*kb*trun/m);


mode_occ=(1/sqrt(2))*num_tot*(corr_widthx*corr_widthy*corr_widthz/(therm_width)^3);
corr.window=[corr_widthx,corr_widthy,corr_widthz]/2;

corr.redges=linspace(0,corr_widthx*4,40);
corr.xedges=linspace(-corr_widthx*5,corr_widthx*5,40);

num=round(num_tot/(1+corr_frac));
corr.chunk_factor=1;
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
for n=1:shots
    shot=reshape(randn(num*3,1)*therm_width,[num,3]);
    coor_chance=rand(num,1)<corr_frac;
    if sum(coor_chance)>0
        delpos=reshape(randn(sum(coor_chance)*3,1),[sum(coor_chance),3]).*repmat([corr_widthz;corr_widthx;corr_widthy]',[sum(coor_chance),1]);
        shot_with_corr=[shot;shot(coor_chance,:)+delpos];
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

CalcCorr(counts_txy,corr)

toc;
