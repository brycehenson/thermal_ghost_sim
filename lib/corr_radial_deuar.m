function out=corr_radial_deuar(corr_opts,counts_txy)
% deuar_corr - %a low/high memory implementation of the deuar correlator
% Notes:
% This is a 'non traditional' correlaor that is used to see how 'bloby' the data is. Uses the distribution of the
% ordered differeces. eg the mean distance of the 4th nearest data point.
%
% Usage:
% How to interpret the results of this correlator
% ? How to interpret the width, can this be converted to blob width.
% ? can it be used to find gaps in data
%
% Syntax:  out=corr_1d_cart(corr_opts,counts_txy)
%
% Inputs:
%	counts_txy                  - as a cell array of zxy counts size=[num_counts,3]
%               must be ordered in z for pre window optimzation
%   corr_opts                   - structure of input options
%           .one_d_smoothing    - [float] gaussian smoothing for the  coincidence density output use nan or 0 for no
%                                 smoothing
%           .attenuate_counts   - value 0-1 the ammount of data to keep, like simulating a lower QE
%           .do_pre_mask        - logical, use fast sorted mask to window to a range of 
%                                 one_d_window in the sorted axis arround each count
%               .sorted_dir     - must be passed when corr_opts.do_pre_mask used, usualy 1 for time
%           .low_mem            - use lowmem(true) or highmem(false) options, use nan to set dynamicaly using the
%                                 remaining memory on your computer.
%
% Outputs:
%   out                         - output structure
%       .one_d_corr_density     - (optional smoothed) coincidence density
%       .one_d_corr_density_raw - raw coincidence density
%       .x_centers              - centers of the histogram vector
%       .pairs                  - total number of pairs
%
% Example: 
%    [TO DO]
%
% Other m-files required: fast_sorted_mask
% Subfunctions: none
% MAT-files required: none
% See also: corr_unit_testing, calc_any_g2_type, import_mcp_tdc_data
%
% Known BUGS/ Possible Improvements
%   - decrease memory requirements
%   - implement standard deviation
%   - implemnt histogram of the nth distance
%   - understand how to normalize with different atom numbers
%   - implement low memory option with mean/histogram on the fly
%     - need to figure out how to normalize shots first
%     - do after the high memory option is working

%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-12-18

%------------- BEGIN CODE --------------

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

if ~isfield(corr_opts,'progress_updates')
    corr_opts.progress_updates=50;
end

num_counts=cellfun(@(x)size(x,1),counts_txy);
%if ~isfield(corr_opts,'low_mem') || isnan(corr_opts.low_mem)
%TO BE USED WHEN LOW MEM IS IMPLEMENTED   
%     mem_temp=memory;
%     max_arr_size=floor(mem_temp.MaxPossibleArrayBytes/(8*2)); %divide by 8 for double array size and 2 as a saftey factor
%     %premasking can dramaticaly reduce the number of pairs that are stored in memory
%     %this might be a factor of ~1/100
%     if corr_opts.do_pre_mask
%         premask_factor=1e-2; 
%     else
%         premask_factor=1;
%     end
%     %use the corr_opts.norm_samp_factor so that both corr/uncorr parts are calculated with the same function
%     corr_opts.low_mem=((max(num_counts)*corr_opts.norm_samp_factor*premask_factor)^2)>max_arr_size;
% 
%     if corr_opts.low_mem,fprintf('auto detection:using the low memory option\n'), end
%     if ~corr_opts.low_mem,fprintf('auto detection:using the high memory option\n'), end
%end
corr_opts.low_mem=false;

shots =size(counts_txy,2);
updates=corr_opts.progress_updates; %number of updates in the progress bar to give the user, can slow thigs down if too high
update_interval=ceil(shots/updates);
parfor_progress_imp(ceil(shots/update_interval));

pairs_count=zeros(1,shots);
delta_multiplier=[1,-1];
delta_multiplier=delta_multiplier(1+corr_opts.cl_or_bb); %gives 1 when cl and -1 when bb
       
ordered_delta=cell(shots,1); %each cell should countain an num_counts by num_counts array with the ordered radial differences

    
if corr_opts.low_mem %calculate with the low memory mode
    error('low mem option not yet implemented')
    % will require some more thinking , high memory is easier to implement 
    
else%calculate with the high memory mode
    for shotnum=1:shots
        shot_txy=counts_txy{shotnum};
        num_counts_shot=num_counts(shotnum);
        ordered_delta_shot=nan(num_counts_shot,num_counts_shot);
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
            if corr_opts.do_pre_mask  %pre mask optimzation using sortd masking
                error('do_pre_mask is not implemented yet')
%                 temp_1d_diff=shot_txy(ii+1:num_counts_shot,corr_opts.sorted_dir)...
%                     -delta_multiplier*shot_txy(ii,corr_opts.sorted_dir);
%                 mask_idx=fast_sorted_mask(...
%                     temp_1d_diff,...
%                     corr_opts.one_d_window(corr_opts.sorted_dir,1),...
%                     corr_opts.one_d_window(corr_opts.sorted_dir,2));
%                 delta=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+mask_idx(1):ii+mask_idx(2),:);
            else
                %take the difference to all counts (including itself) %use [1:ii-1,ii+1:num_counts_shot] to exclude
                delta=shot_txy(ii,:)-delta_multiplier*shot_txy;
            end
            rad_delta=sqrt(sum(delta.^2,2));
            rad_delta=sort(rad_delta);
            ordered_delta_shot(ii,:)=rad_delta';
        end
        if mod(shotnum,update_interval)==0
            parfor_progress_imp;
        end
        ordered_delta{shotnum}=ordered_delta_shot;
    end%loop over shots
parfor_progress_imp(0);

%%
diffs_to_store=5e2;
mean_nth_rad_diff=nan(diffs_to_store,shots);
for shotnum=1:shots
    ordered_delt_shot=ordered_delta{shotnum};
    diffs_avail=size(ordered_delt_shot,2);
    padd_size=max(0,diffs_to_store-diffs_avail);
    max_idx=min(diffs_to_store,diffs_avail);
    %padd up to the requested size
    mean_nth_rad_diff(:,shotnum)=cat(2,nanmean(ordered_delt_shot(:,1:max_idx),1),nan(1,padd_size));
end

mean_nth_rad_diff=nanmean(mean_nth_rad_diff,2);
radial_dep=((1:size(mean_nth_rad_diff,1)).*(3/(4*pi))).^(1/3);

mean_nth_rad_diff_reg=mean_nth_rad_diff./radial_dep';

%%


end 



out.mean_nth_rad_diff_raw=mean_nth_rad_diff;
out.mean_nth_rad_diff_reg=mean_nth_rad_diff_reg;


end


% %%you can check that the sorted mask this is equivelent to the below brute approach
% delta_brute=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+1:num_counts_shot,:); 
% min_lim=corr_opts.one_d_window(corr_opts.sorted_dir,1);
% max_lim=corr_opts.one_d_window(corr_opts.sorted_dir,2);
% mask=delta_brute(:,corr_opts.sorted_dir)>min_lim & delta_brute(:,corr_opts.sorted_dir)<max_lim;
% delta_brute =delta_brute(mask,:);
% if ~isequal(delta_brute,delta)
%     fprintf('warning mask is not giving correct answer')
% end

% out.x_centers=(corr_opts.one_d_edges(2:end)+corr_opts.one_d_edges(1:end-1))/2;
% out.pairs=sum(pairs_count);
% 
% sub_index=[1,2,3];
% sub_index(corr_opts.one_d_dimension)=[];
% one_d_volume=corr_opts.one_d_window(sub_index,1)-corr_opts.one_d_window(sub_index,2);
% one_d_volume=prod(one_d_volume)*(corr_opts.one_d_edges(2:end)-corr_opts.one_d_edges(1:end-1));
% out.one_d_corr_density_raw=sum(one_d_bins,1)./(one_d_volume.*out.pairs);
% %smooth the correlation function for better normalization
% if ~(isnan(corr_opts.one_d_smoothing) || corr_opts.one_d_smoothing==0)
%     out.one_d_corr_density=gaussfilt(out.x_centers,out.one_d_corr_density_raw,corr_opts.one_d_smoothing);
% else
%     out.one_d_corr_density=out.one_d_corr_density_raw;
% end
