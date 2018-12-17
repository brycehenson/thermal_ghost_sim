function out=corr_1d_cart_low_mem_basic(corr_opts,counts_txy)
%a low memory implementation of the one D correlation
%Inputs
%   counts_txy- as a cell array of zxy counts size=[num_counts,3]
%               must be ordered in z for pre window optimzation
%   corr_opts.one_d_window=[[tmin,tmax];[xmin,xmax];[ymin,ymax]];
%   corr_opts.cl_or_bb
%   corr_opts.one_d_edges
%   corr_opts.one_d_dimension
%   corr_opts.attenuate_counts- value 0-1 the ammount of data 
%   corr_opts.do_pre_mask- logical, use fast sorted mask to window to a range of 
%                           one_d_window in the sorted axis arround each count
%       corr_opts.sorted_dir - must be passed when corr_opts.do_pre_mask used
%pass in a cell array of counts and then calculate the correlations using a
%bin on the fly method
%for the the gxyz we require the difference in the other axes to be less
%than a threshold

%improvements
% implement pre delta windowing for time sorted data using binary search
% input checking
% kill counts, to speed things up when there is a lot of data remove some fraction

if ~isfield(corr_opts,'attenuate_counts')
    corr_opts.attenuate_counts=1;
else
    if corr_opts.attenuate_counts>1 || corr_opts.attenuate_counts<0
        error('invalid attenuation value, pass a value between 0 and 1\n')
    end
end

if ~isfield(corr_opts,'cl_or_bb')
    error('corr_opts.cl_or_bb not specified!do you want co-linear or back-to-back? ')
end

if isfield(corr_opts,'do_pre_mask') && ~isfield(corr_opts,'sorted_dir')
        error('you must pass corr_opts.sorted_dir with corr_opts.do_pre_mask')
elseif ~isfield(corr_opts,'do_pre_mask')
    corr_opts.do_pre_mask=false;
end

%for the dimension of interest we set the window to be the redion specified in corr_opts.one_d_edges
%this simplifies the code and can speed up the historgraming
corr_opts.one_d_window(corr_opts.one_d_dimension,:)=[min(corr_opts.one_d_edges),max(corr_opts.one_d_edges)];

updates=corr_opts.progress_updates; %number of updates in the progress bar to give the user, can slow thigs down if too high

%find un-norm g2
shots =size(counts_txy,2);
update_interval=ceil(shots/updates);
parfor_progress_imp(ceil(shots/update_interval));

[one_d_bins,corr_opts.one_d_edges]=histcounts([],corr_opts.one_d_edges);
one_d_bins=zeros(shots,size(one_d_bins,2));

pairs_count=zeros(1,shots);
for shotnum=1:shots
    shot_txy=counts_txy{shotnum};
    num_counts_shot=size(shot_txy,1);
    if corr_opts.attenuate_counts~=1 %randomly keep corr_opts.attenuate_counts fraction of the data
        mask=rand(num_counts_shot,1)<corr_opts.attenuate_counts;
        shot_txy=shot_txy(mask,:);
        num_counts_shot=size(shot_txy,1); %recalulate the number of counts
    end
    if num_counts_shot<2
        warning('%u counts input\n',num_counts_shot)
    end
    % full number of pairs in the shot
    pairs_count(shotnum)=num_counts_shot^2 -num_counts_shot;
    for ii=1:num_counts_shot
        %if each shot is sorted in time then this will only produce counts that are one direction in time
        %pre mask optimzation here
        if ~corr_opts.cl_or_bb %colinear is when cl_or_bb=false
                delta=shot_txy(ii,:)-shot_txy(ii+1:num_counts_shot,:);
        else %back to back correlation
                delta=shot_txy(ii,:)+shot_txy(ii+1:num_counts_shot,:);
        end
        
        %one dimensional cartesian correlations
        %mask in 2d for the 1d correlation, this is pretty messy but it avoids another comparison
        switch corr_opts.one_d_dimension
            case 1
                one_d_mask_pos=delta(:,2)>corr_opts.one_d_window(2,1) & delta(:,2)<corr_opts.one_d_window(2,2) &...
                delta(:,3)>corr_opts.one_d_window(3,1) & delta(:,3)<corr_opts.one_d_window(3,2);
                one_d_mask_neg=-delta(:,2)>corr_opts.one_d_window(2,1) & -delta(:,2)<corr_opts.one_d_window(2,2) &...
                -delta(:,3)>corr_opts.one_d_window(3,1) & -delta(:,3)<corr_opts.one_d_window(3,2);
            case 2 
                one_d_mask_pos=delta(:,1)>corr_opts.one_d_window(1,1) & delta(:,1)<corr_opts.one_d_window(1,2) &...
                delta(:,3)>corr_opts.one_d_window(3,1) & delta(:,3)<corr_opts.one_d_window(3,2);
                one_d_mask_neg=-delta(:,1)>corr_opts.one_d_window(1,1) & -delta(:,1)<corr_opts.one_d_window(1,2) &...
                -delta(:,3)>corr_opts.one_d_window(3,1) & -delta(:,3)<corr_opts.one_d_window(3,2);
            case 3
                one_d_mask_pos=delta(:,1)>corr_opts.one_d_window(1,1) & delta(:,1)<corr_opts.one_d_window(1,2) &...
                delta(:,2)>corr_opts.one_d_window(2,1) & delta(:,2)<corr_opts.one_d_window(2,2);
                one_d_mask_neg=-delta(:,1)>corr_opts.one_d_window(1,1) & -delta(:,1)<corr_opts.one_d_window(1,2) &...
                -delta(:,2)>corr_opts.one_d_window(2,1) & -delta(:,2)<corr_opts.one_d_window(2,2);
        end
        
        %to be strictly accurate we must calaulate things symetricaly
        one_d_bins(shotnum,:)=one_d_bins(shotnum,:)+histcounts(delta(one_d_mask_pos,corr_opts.one_d_dimension),corr_opts.one_d_edges);
        one_d_bins(shotnum,:)=one_d_bins(shotnum,:)+histcounts(-delta(one_d_mask_neg,corr_opts.one_d_dimension),corr_opts.one_d_edges);
    end
    if mod(shotnum,update_interval)==0
        parfor_progress_imp;
    end
end
parfor_progress_imp(0);

out.x_centers=(corr_opts.one_d_edges(2:end)+corr_opts.one_d_edges(1:end-1))/2;
out.pairs=sum(pairs_count);

sub_index=[1,2,3];
sub_index(corr_opts.one_d_dimension)=[];
one_d_volume=corr_opts.one_d_window(sub_index,1)-corr_opts.one_d_window(sub_index,2);
one_d_volume=prod(one_d_volume)*(corr_opts.one_d_edges(2:end)-corr_opts.one_d_edges(1:end-1));
out.one_d_corr_density_raw=sum(one_d_bins,1)./(one_d_volume.*out.pairs);
%smooth the correlation function for better normalization
if ~(isnan(corr_opts.one_d_smoothing) || corr_opts.one_d_smoothing==0)
    out.one_d_corr_density=gaussfilt(out.x_centers,out.one_d_corr_density_raw,corr_opts.one_d_smoothing);
else
    out.one_d_corr_density=out.one_d_corr_density_raw;
end


end