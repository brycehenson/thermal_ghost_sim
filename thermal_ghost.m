
addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path

%% fake dataset
%fake_data=fake_cl_corr(100,1);
%fake_data=fake_therm_cl_corr_small_hot(5000,1)
%fake_data=fake_therm_cl_corr_small_cold(5000,1);
%fake_data=fake_super_cl_corr(1000,1);
%fake_data=fake_cl_corr_medium_hot(5000,0.1);
data_faker=@() fake_super_cl_corr(5000,1);

fake_data=data_faker();

total_counts=sum(fake_data.num_counts);
%% Spliting data
fprintf('spliting data\n')
data_up_down={};
for ii=1:size(fake_data.counts_txy,2)
    shot_txy=fake_data.counts_txy{ii};
    mask=rand(size(shot_txy,1),1)>0.5;
    data_up_down.up_counts{ii}=shot_txy(mask,:);
    data_up_down.down_counts{ii}=shot_txy(~mask,:);
end

%% make windowing text function
fprintf('calculating text mask\n')
windowing_function=text_to_window_fun([],0);
%%
fprintf('applying text mask\n')
for ii=1:size(data_up_down.up_counts,2)
    down_counts=data_up_down.down_counts{ii};
    mask=windowing_function(down_counts);
    %simple mask
    %mask=(down_counts(:,2)<0 & down_counts(:,3)<0) | (down_counts(:,2)>0 & down_counts(:,3)>0 & down_counts(:,3)<0.5)  ;
    data_up_down.down_counts_mask{ii}=down_counts(mask,:);
end


%% Verify the mask has been applied
x_edges=linspace(-1,1,50);
y_edges=linspace(-1,1,50);

all_mask_counts=vertcat(data_up_down.down_counts_mask{:});
bins=histcn(all_mask_counts(:,2:3),x_edges,y_edges);
sfigure(1)
clf
subplot(1,3,1)
imagesc(bins')
set(gca,'YDir','normal')
title('Mask port all counts')
colormap(viridis)
pause(1e-3)
%% Correlate

t_edges=linspace(-1,1,300);
x_edges=linspace(-1,1,50);
y_edges=linspace(-1,1,50);

conditional_hist_mask=histcn([nan,nan,nan],t_edges,x_edges,y_edges);
conditional_hist_mask=squeeze(conditional_hist_mask(1,:,:));
conditional_hist_nomask_faked=conditional_hist_mask;
num_x_bins=size(conditional_hist_mask,2);
num_y_bins=size(conditional_hist_mask,3);

x_cen=(x_edges(2:end)+x_edges(1:end-1))/2;
y_cen=(y_edges(2:end)+y_edges(1:end-1))/2;
iimax=size(data_up_down.up_counts,2);
fprintf('Correlating shots with mask %04u:%04u',iimax,0)
for ii=1:iimax%loop over every shot
    
    %spatialy resolved detection
    up_counts=data_up_down.up_counts{ii};
    bins_up=histcn(up_counts,t_edges,x_edges,y_edges);
    
    %condition the up 3d hist on the counts in the down in the single time dimension
    down_counts=data_up_down.down_counts_mask{ii};
    bins_down=histcn(down_counts(:,1),t_edges);
    %bins_down=bins_down*0+1;
    mask_down=bins_down>0;
    condition_mask=repmat(mask_down,[1,num_x_bins,num_y_bins]);
    conditional_hist_inc=bins_up.*condition_mask;
    conditional_hist_inc=squeeze(sum(conditional_hist_inc,1));
    conditional_hist_mask=conditional_hist_mask+conditional_hist_inc;
    
    down_counts=data_up_down.down_counts{ii};
    bins_down=histcn(down_counts(:,1),t_edges);
    %bins_down=bins_down*0+1;
    mask_down=bins_down>0;
    condition_mask=repmat(mask_down,[1,num_x_bins,num_y_bins]);
    conditional_hist_inc=bins_up.*condition_mask;
    conditional_hist_inc=squeeze(sum(conditional_hist_inc,1));
    conditional_hist_nomask_faked=conditional_hist_nomask_faked+conditional_hist_inc;
    if mod(ii,10)==0, fprintf('\b\b\b\b%04u',ii), end
end
fprintf('...Done\n')


fprintf('total counts in conditional hist %u\n',sum(conditional_hist_mask(:)))
sfigure(2)
clf
set(gcf,'color','w')
subplot(2,2,1)
imagesc(x_cen,y_cen,conditional_hist_mask')
set(gca,'YDir','normal')
colormap(viridis)
subplot(2,2,2)
surf(conditional_hist_mask')
colormap(viridis)
pause(1e-3)
%set(gca,'XDir','reverse')

%%
sfigure(1)
conditional_hist_mask_cheating=imgaussfilt(conditional_hist_mask,1)./imgaussfilt(conditional_hist_nomask_faked,1);
subplot(1,3,2)
imagesc(x_cen,y_cen,conditional_hist_mask_cheating')
set(gca,'YDir','normal')
title('same data norm')
colormap(viridis)
subplot(1,3,3)
surf(conditional_hist_mask_cheating')
colormap(viridis)
pause(1e-3)


%% Normalize without the mask in
fake_data=data_faker();
total_counts=sum(fake_data.num_counts);
%%
data_up_down={};
for ii=1:size(fake_data.counts_txy,2)
    shot_txy=fake_data.counts_txy{ii};
    mask=rand(size(shot_txy,1),1)>0.5;
    data_up_down.up_counts{ii}=shot_txy(mask,:);
    data_up_down.down_counts{ii}=shot_txy(~mask,:);
end

%% Correlate without the mask
conditional_hist_nomask_real=histcn([nan,nan,nan],t_edges,x_edges,y_edges);
conditional_hist_nomask_real=squeeze(conditional_hist_nomask_real(1,:,:));
conditional_hist_nomask_real=conditional_hist_nomask_real;
num_x_bins=size(conditional_hist_nomask_real,2);
num_y_bins=size(conditional_hist_nomask_real,3);

x_cen=(x_edges(2:end)+x_edges(1:end-1))/2;
y_cen=(y_edges(2:end)+y_edges(1:end-1))/2;
iimax=size(data_up_down.up_counts,2);
fprintf('Correlating shots without mask %04u:%04u',iimax,0)
for ii=1:size(data_up_down.up_counts,2)
    %spatialy resolved detection
    up_counts=data_up_down.up_counts{ii};
    bins_up=histcn(up_counts,t_edges,x_edges,y_edges);
    
    %condition the up 3d hist on the counts in the down in the single time dimension
    down_counts=data_up_down.down_counts{ii};
    bins_down=histcn(down_counts(:,1),t_edges);
    %bins_down=bins_down*0+1;
    mask_down=bins_down>0;
    condition_mask=repmat(mask_down,[1,num_x_bins,num_y_bins]);
    conditional_hist_inc=bins_up.*condition_mask;
    conditional_hist_inc=squeeze(sum(conditional_hist_inc,1));
    conditional_hist_nomask_real=conditional_hist_nomask_real+conditional_hist_inc;
    if mod(ii,10)==0, fprintf('\b\b\b\b%04u',ii), end
end
fprintf('...Done\n')
%%

%conditional_hist_nomask_faked
sfigure(2)

cl=[0.25,0.325];
conditional_hist_mask_norm=imgaussfilt(conditional_hist_mask,1)./imgaussfilt(conditional_hist_nomask_real,1);
subplot(2,2,3)
imagesc(x_cen,y_cen,conditional_hist_mask_norm')
set(gca,'YDir','normal')
%caxis(cl)
colormap(viridis)
subplot(2,2,4)
surf(x_cen,y_cen,conditional_hist_mask_norm')
colormap(viridis)
%caxis(cl)
pause(1e-3)
%subplot(2,2,3)
%surf(x_cen,y_cen,conditional_hist_mask_norm)


