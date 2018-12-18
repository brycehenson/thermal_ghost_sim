%deuar correlator unit testing

%% Setting Up The Enviorment
%add all subfolders to the path
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));

%%
fake_data=make_bloby_data(100,1);

% Verify blob
sfigure(2);
subplot(1,2,1)
scatter3(fake_data.counts_txy{1}(:,1),fake_data.counts_txy{1}(:,2),fake_data.counts_txy{1}(:,3))

%%
corr_opts=[];
corr_opts.cl_or_bb=false;
corr_opts.progress_updates=10;
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

%%
sfigure(1);
subplot(3,1,1)
plot(deuar_corr.mean_nth_rad_diff_norm(2:end))
ylabel('regularized distance')
xlabel('nth nearest point')
pause(1e-3)
sfigure(1)
subplot(3,1,2)
plot(deuar_norm.mean_nth_rad_diff_norm(2:end))
ylabel('regularized distance')
xlabel('nth nearest point')
subplot(3,1,3)
plot(deuar_norm.mean_nth_rad_diff_reg(2:end)./deuar_corr.mean_nth_rad_diff_reg(2:end))
xlabel('nth nearest point')
ylabel('normalized distance')