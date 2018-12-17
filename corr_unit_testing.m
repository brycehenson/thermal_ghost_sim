%corr unit testing
% main questions
% should the in shot pairs be counted in the number of possible normalization pairs

%% Setting Up The Enviorment
%add all subfolders to the path
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));


%% fake dataset
%fake_data=fake_cl_corr(100,1);
%fake_data=fake_therm_cl_corr_small_hot(5000,1)
%fake_data=fake_therm_cl_corr_small_cold(1000,1);
%fake_data=fake_super_cl_corr(100,1);
%fake_data=fake_cl_corr_medium_hot(100,0.1); %proof of principle for ghost imaging
fake_data=fake_cl_corr(100,0.1);

%look at the histogram of the density
rmax=0.2;
rmin=1e-5;
corr_opts.redges=sqrt(linspace(rmin^2,rmax^2,500));
%corr_opts.redges=linspace(rmin,rmax,500);

rad_volume=(4/3)*pi*(corr_opts.redges(2:end).^3-corr_opts.redges(1:end-1).^3);
all_txy=vertcat(fake_data.counts_txy{:});
rad=sqrt(sum(all_txy.^2,2));
rad_centers=(corr_opts.redges(2:end)+corr_opts.redges(1:end-1))/2;
[rad_bins,corr_opts.redges]=histcounts(rad,corr_opts.redges);
figure(3)
clf
plot(rad_centers,rad_bins./rad_volume)



%% Test the new any_g2_type function
%fake_data=fake_cl_corr(100,1);
fake_data=fake_cl_corr_medium_hot(5000,0.1)
%%
corr_opts=[];
corr_opts.type='1d_cart_cl';
corr_opts.one_d_dimension=1;
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*5e-2;
one_d_range=0.05;
corr_opts.one_d_edges=linspace(-one_d_range,one_d_range,300);

corr_opts.low_mem=nan;
corr_opts.plots=true;
corr_opts.norm_samp_factor=2;
corr_opts.attenuate_counts=1;
corr_opts.do_pre_mask=true;
corr_opts.sorted_dir=1;
corr_opts.sort_norm=true;


corr_opts.one_d_smoothing=nan;
%improved code  
%28.4 s with premask & sort chunks 
%30.12  with premask & no sort chunks 
%28.57 new with no premask
%28.61 s old with no premask

% old code 27.137 s
tic
out=calc_any_g2_type(corr_opts,fake_data);
toc


%%

corr_opts.type='radial_cl';
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*5e-2; %only used for prewindow
rmax=0.05;
rmin=1e-4;
corr_opts.redges=sqrt(linspace(rmin^2,rmax^2,600));
%corr_opts.redges=linspace(rmin,rmax,500);
corr_opts.rad_smoothing=0;

corr_opts.attenuate_counts=1;
corr_opts.one_d_smoothing=nan;
%improved code  
%28.4 s with premask & sort chunks 
%30.12  with premask & no sort chunks 
%28.57 new with no premask
%28.61 s old with no premask

% old code 27.137 s
tic
out=calc_any_g2_type(corr_opts,fake_data);
toc



%% vvvvvvvOLD STUFF vvvvvv


corr_opts=[];
corr_opts.one_d_dimension=1;
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*5e-2;
one_d_range=0.6;
corr_opts.one_d_edges=linspace(-one_d_range,one_d_range,80);
corr_opts.one_d_smoothing=0.025;

corr_opts.do_three_d_corr=false;
three_d_num_edge=50;
corr_opts.three_d_edges={linspace(corr_opts.one_d_window(1,1),corr_opts.one_d_window(1,2),three_d_num_edge),...
                        linspace(corr_opts.one_d_window(2,1),corr_opts.one_d_window(2,2),three_d_num_edge),...
                        linspace(corr_opts.one_d_window(3,1),corr_opts.one_d_window(3,2),three_d_num_edge)};

rmax=0.6;
rmin=1e-3;
%corr_opts.redges=linspace(rmin,rmax,40);
corr_opts.redges=sqrt(linspace(rmin^2,rmax^2,500));
corr_opts.rad_smoothing=0.025;


fprintf('calculating in shot correlations \n')
tic
shotscorr=corr_multi_low_mem(corr_opts,fake_data.counts_txy);
toc

sfigure(1);
clf
set(gcf,'color','w');
subplot(2,3,1)
plot(shotscorr.rad_centers,shotscorr.rad_corr_density,'.k-','MarkerSize',10)
title('In Shot Rad Dist')
xlabel('Radial Seperation')
ylabel('G^2 coincedence density')
subplot(2,3,4)
plot(shotscorr.x_centers,shotscorr.one_d_corr_density,'.k-','MarkerSize',10)
title('In Shot X Dist (windowed)')
ylabel('G^2 coincedence density')
xlabel('X Seperation')



%% chunk up the data for normalization
norm_samp_factor=20; %should be about the correlation amp for equal noise contibution from in-shot & between shot
norm_chunk_size=round(mean(fake_data.num_counts)*sqrt(norm_samp_factor));

all_counts=[vertcat(fake_data.counts_txy{:})];

all_counts=all_counts(randperm(total_counts),:);
counts_chunked={};
iimax=ceil(total_counts/norm_chunk_size);
for ii=1:iimax
    if ii~=iimax
        counts_chunked{ii}=all_counts((ii-1)*norm_chunk_size+1:ii*norm_chunk_size,:);
    else
        counts_chunked{ii}=all_counts((ii-1)*norm_chunk_size+1:end,:);
    end
end
if size(vertcat(counts_chunked{:}),1)~=total_counts
    warning('lost counts')
end

%%
fprintf('calculating inter-shot correlations \n ')
normcorr=corr_multi_low_mem(corr_opts,counts_chunked);

%%
subplot(2,3,2)
plot(normcorr.rad_centers,normcorr.rad_corr_density,'.k-','MarkerSize',10)
title('Between Shot Rad Dist')
xlabel('Radial Seperation')
ylabel('G^2 coincedence density')
subplot(2,3,5)
plot(normcorr.x_centers,normcorr.one_d_corr_density,'.k-','MarkerSize',10)
title('Between Shot X Dist (windowed)')
ylabel('G^2 coincedence density')
xlabel('X Seperation')
subplot(2,3,6)
xg2=shotscorr.one_d_corr_density./normcorr.one_d_corr_density;
plot(shotscorr.x_centers,xg2,'.k-','MarkerSize',10)
title('Norm. Corr.')
ylabel('g^{(2)} (X)')
xlabel('X Seperation')
subplot(2,3,3)
radg2=shotscorr.rad_corr_density./normcorr.rad_corr_density; 
plot(shotscorr.rad_centers,radg2,'.k-','MarkerSize',10)
ylabel('g^{(2)} (R)')
title('Norm. Corr.')
xlabel('Radial Seperation')
pause(0.1);
%saveas(gcf,'Correlations.png')
%saveas(gcf,'Correlations.pdf')
%savefig('Correlations.fig')
%saveas(gcf,['Correlations',num2str(num),'.png'])
%saveas(gcf,['Correlations',num2str(num),'.pdf'])
fprintf('max g2 rad       %4.2f \n',max(radg2))
fprintf('       x         %4.2f \n',max(xg2))
toc;



