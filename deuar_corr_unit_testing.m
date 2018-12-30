%deuar correlator unit testing

%% Setting Up The Enviorment
%add all subfolders to the path
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));

%%
fake_data=make_bloby_data(300,1);

% Verify blob
sfigure(2);
cla
set(gcf,'color','w')
subplot(1,2,1)
scatter3(fake_data.counts_txy{1}(:,1),fake_data.counts_txy{1}(:,2),fake_data.counts_txy{1}(:,3))
title('single shot')
%%
corr_opts=[];
corr_opts.cl_or_bb=false;
corr_opts.progress_updates=10;
corr_opts.diff_bins=linspace(1e-3,0.4,1e2);
corr_opts.diffs_to_store=120;

deuar_corr=corr_radial_deuar(corr_opts,fake_data.counts_txy);
sfigure(1);
set(gcf,'color','w')
subplot(3,1,1)
plot(deuar_corr.mean_nth_rad_diff_reg(2:end))
ylabel('regularized distance')
xlabel('nth nearest point')
%%
norm_samp_factor=1;
counts_chunked=chunk_data(fake_data,norm_samp_factor,NaN);
deuar_norm=corr_radial_deuar(corr_opts,counts_chunked);

%%
sfigure(2);
subplot(1,2,2)
scatter3(counts_chunked{1}(:,1),counts_chunked{1}(:,2),counts_chunked{1}(:,3))
title('Cross shots')
%%
sfigure(1);
subplot(3,1,1)
plot(deuar_corr.mean_nth_rad_diff_reg(2:end))
ylabel('regularized distance')
xlabel('nth nearest point')
pause(1e-3)
sfigure(1)
subplot(3,1,2)
plot(deuar_norm.mean_nth_rad_diff_reg(2:end))
ylabel('regularized distance')
xlabel('nth nearest point')
subplot(3,1,3)
plot(deuar_norm.mean_nth_rad_diff_reg(2:end)./deuar_corr.mean_nth_rad_diff_reg(2:end))
xlabel('nth nearest point')
ylabel('normalized distance')


%%
sfigure(3)
set(gcf,'color','w')
subplot(3,1,1)
imagesc(deuar_corr.diff_bin_censx(2:end),deuar_corr.diff_bin_censy,deuar_corr.distance_hist(2:end,:)')
ylabel('distance')
xlabel('j^{th} nearest point')
set(gca,'YDir','normal')
colormap(viridis)
subplot(3,1,2)
imagesc(deuar_norm.diff_bin_censx(2:end),deuar_norm.diff_bin_censy,deuar_norm.distance_hist(2:end,:)')
ylabel('distance')
xlabel('j^{th} nearest point')
set(gca,'YDir','normal')
colormap(viridis)
subplot(3,1,3)
start_idx=2;
dist_hist_norm=(deuar_corr.distance_hist(start_idx:end,:)-deuar_norm.distance_hist(start_idx:end,:));
imagesc(deuar_norm.diff_bin_censx(start_idx:end),deuar_norm.diff_bin_censy,dist_hist_norm')
ylabel('distance')
xlabel('j^{th} nearest point')
set(gca,'YDir','normal')
colormap(viridis)

%%
sfigure(4)
clf
hold on
nth_points=30
colors=viridis(nth_points); 
for ii=2:nth_points
    plot(deuar_norm.diff_bin_censy,deuar_norm.distance_hist(ii,:),'color',colors(ii,:))
end
hold off

%% try some fitting
ii=3;
xdat=deuar_norm.diff_bin_censy;
ydat=deuar_norm.distance_hist(ii,:);
% modelfun=@(b,x) b(1)*(1/(4*b(3)))*(x.^(b(5))).*(sech(((x.^(b(4)))-b(2))/b(3))).^2;
% param_guess=[10,0.2,0.2,1,1];
% modelfun=@(b,x) b(1)*(1/(4*b(3)))*(x.^2).*(sech(((x.^2)-b(2))/b(3))).^2 ;
% param_guess=[10,0.2,0.2];
% modelfun=@(b,x) (b(1)*b(2)/b(3))*sqrt(pi/2)*exp((1/2)*((b(2)/b(3))^2)-((x-b(4))/b(3))).*...
%     erfc((1/sqrt(2))*(b(2)/sqrt(2)-((x-b(4))/b(3))));
% param_guess=[1e-14,0.3,0.04,0.18];
modelfun=@(b,x) (b(1)*b(2)/b(3))*sqrt(pi/2)*exp((1/2)*((b(2)/b(3))^2)-((x.^b(5)-b(4))/b(3))).*...
    erfc((1/sqrt(2))*(b(2)/sqrt(2)-((x.^b(5)-b(4))/b(3))));
param_guess=[1e-14,0.3,0.04,0.18,1];
modelfun=@(b,x) b(1).*exp(-b(2)*(x.^3)).*b(2).*x.^2;
param_guess=[1e4,600];
opts = statset('MaxIter',1e2);
fit_model=fitnlm(xdat,ydat,modelfun,param_guess,'Options',opts)
xsamp=linspace(min(xdat),max(xdat),1e3)';
ysamp=predict(fit_model,xsamp);
plot(xsamp,ysamp)
hold on
plot(xdat,ydat,'.')
hold off
set(gca,'yscale','log')
%set(gca,'xscale','log')
set(gcf,'color','w')
xlabel('Distance (m)')
ylabel('Density')
title(sprintf('j=%u Neighbour Distance',ii-1))

