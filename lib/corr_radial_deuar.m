function out=corr_radial_deuar(corr_opts,counts_txy)
% deuar_corr - %a low/high memory implementation of the deuar correlator
% Notes:
% This is a 'non traditional' correlaor that is used to see how 'bloby' the data is
%
% Usage:s
% pass in a cell array of counts and then calculate the correlations using a
% bin on the fly method (low mem) or a bin all diferences (high mem)
% differences can be prewindowed using fast_sorted_mask (optional) in both cases to reduce memory & comp requirements.
% differences are then windowed in two dimensions and the correlation calculated in the remaining axis.
% conforms to output formats of import_mcp_tdc_data


% corr_1d_cart - %a low/high memory implementation of the one Deuar correlator
% pass in a cell array of counts and then calculate the correlations using a
% difference on the fly method (low mem) or a bin all diferences (high mem)
% differences can be prewindowed using fast_sorted_mask (optional) in both cases to reduce memory & comp requirements.
% differences are then windowed in two dimensions and the correlation calculated in the remaining axis.
% conforms to output formats of import_mcp_tdc_data
%
% Syntax:  out=corr_1d_cart(corr_opts,counts_txy)
%
% Inputs:
%	counts_txy                  - as a cell array of zxy counts size=[num_counts,3]
%               must be ordered in z for pre window optimzation
%   corr_opts                   - structure of input options
%           .one_d_window       - [[tmin,tmax];[xmin,xmax];[ymin,ymax]];
%           .cl_or_bb           - colinear(false) or bb (true)
%           .one_d_edges        - edges of histogram bins for the correlation
%           .one_d_dimension    - dimension to calculate correlation
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
%   - input checking
%	- document output better
%	- more outputs
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

num_counts=cellfun(@(x)size(x,1),counts_txy);
if ~isfield(corr_opts,'low_mem') || isnan(corr_opts.low_mem)
    mem_temp=memory;
    max_arr_size=floor(mem_temp.MaxPossibleArrayBytes/(8*2)); %divide by 8 for double array size and 2 as a saftey factor
    %premasking can dramaticaly reduce the number of pairs that are stored in memory
    %this might be a factor of ~1/100
    if corr_opts.do_pre_mask
        premask_factor=1e-2; 
    else
        premask_factor=1;
    end
    %use the corr_opts.norm_samp_factor so that both corr/uncorr parts are calculated with the same function
    corr_opts.low_mem=((max(num_counts)*corr_opts.norm_samp_factor*premask_factor)^2)>max_arr_size;

    if corr_opts.low_mem,fprintf('auto detection:using the low memory option\n'), end
    if ~corr_opts.low_mem,fprintf('auto detection:using the high memory option\n'), end
end

%for the dimension of interest we set the window to be the region specified in corr_opts.one_d_edges
%this simplifies the code and can speed up the historgraming
corr_opts.one_d_window(corr_opts.one_d_dimension,:)=[min(corr_opts.one_d_edges),max(corr_opts.one_d_edges)];

updates=corr_opts.progress_updates; %number of updates in the progress bar to give the user, can slow thigs down if too high


shots =size(counts_txy,2);
update_interval=ceil(shots/updates);
parfor_progress_imp(ceil(shots/update_interval));

pairs_count=zeros(1,shots);
delta_multiplier=[1,-1];
delta_multiplier=delta_multiplier(1+corr_opts.cl_or_bb); %gives 1 when cl and -1 when bb
       
%build a queue of the dimensions that should be masked in before the histogram
%should not mask in corr_opts.one_d_dimension & corr_opts.sorted_dir (if corr_opts.do_pre_mask )
mask_dimensions=[1,2,3];
mask_dimensions_logic=[true,true,true]; %logic to avoid clashes when corr_opts.one_d_dimension==corr_opts.sorted_dir
mask_dimensions_logic(corr_opts.one_d_dimension)=false;
if corr_opts.do_pre_mask, mask_dimensions_logic(corr_opts.sorted_dir)=false; end
mask_dimensions=mask_dimensions(mask_dimensions_logic);

[one_d_bins,corr_opts.one_d_edges]=histcounts([],corr_opts.one_d_edges);
one_d_bins=zeros(shots,size(one_d_bins,2));
    
if corr_opts.low_mem %calculate with the low memory mode
    for shotnum=1:shots
        shot_txy=counts_txy{shotnum};
        num_counts_shot=num_counts(shotnum);
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
        for ii=1:num_counts_shot-1
            %if each shot is sorted in time then this will only produce counts that are one direction in time
            if corr_opts.do_pre_mask  %pre mask optimzation using sortd masking
                temp_1d_diff=shot_txy(ii+1:num_counts_shot,corr_opts.sorted_dir)...
                    -delta_multiplier*shot_txy(ii,corr_opts.sorted_dir);
                mask_idx=fast_sorted_mask(...
                    temp_1d_diff,...
                    corr_opts.one_d_window(corr_opts.sorted_dir,1),...
                    corr_opts.one_d_window(corr_opts.sorted_dir,2));
                delta=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+mask_idx(1):ii+mask_idx(2),:);
            else
                delta=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+1:num_counts_shot,:);
            end
            
            one_d_mask_pos=true(size(delta,1),1); %initalize the mask
            one_d_mask_neg=one_d_mask_pos;
            %run the mask in all the queued dimensions
            for mask_dim=mask_dimensions
                one_d_mask_pos=one_d_mask_pos & delta(:,mask_dim)>corr_opts.one_d_window(mask_dim,1) ...
                    & delta(:,mask_dim)<corr_opts.one_d_window(mask_dim,2);
                one_d_mask_neg=one_d_mask_neg & -delta(:,mask_dim)>corr_opts.one_d_window(mask_dim,1) ...
                    & -delta(:,mask_dim)<corr_opts.one_d_window(mask_dim,2);
            end

            %to be strictly accurate we must calaulate things symetricaly
            one_d_bins(shotnum,:)=one_d_bins(shotnum,:)+histcounts(delta(one_d_mask_pos,corr_opts.one_d_dimension),corr_opts.one_d_edges);
            one_d_bins(shotnum,:)=one_d_bins(shotnum,:)+histcounts(-delta(one_d_mask_neg,corr_opts.one_d_dimension),corr_opts.one_d_edges);
        end
        if mod(shotnum,update_interval)==0
            parfor_progress_imp;
        end
    end%loop over shots
    
    
    
else%calculate with the high memory mode
    for shotnum=1:shots
        shot_txy=counts_txy{shotnum};
        num_counts_shot=num_counts(shotnum);
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
        upper_triangle_size=pairs_count(shotnum)/2;
        delta_idx=1;
        delta=nan(upper_triangle_size,3);
        for ii=1:num_counts_shot-1
            %if each shot is sorted in time then this will only produce counts that are one direction in time
            %pre mask optimzation here
            if corr_opts.do_pre_mask  %pre mask optimzation using sortd masking
                    temp_1d_diff=shot_txy(ii+1:num_counts_shot,corr_opts.sorted_dir)...
                        -delta_multiplier*shot_txy(ii,corr_opts.sorted_dir);
                    mask_idx=fast_sorted_mask(...
                        temp_1d_diff,...
                        corr_opts.one_d_window(corr_opts.sorted_dir,1),...
                        corr_opts.one_d_window(corr_opts.sorted_dir,2));
                    delta_inc=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+mask_idx(1):ii+mask_idx(2),:);
                else
                    delta_inc=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+1:num_counts_shot,:);
            end
            delta_inc_size=size(delta_inc,1);
            delta(delta_idx:delta_idx+delta_inc_size-1,:)=delta_inc;
            delta_idx=delta_idx+delta_inc_size; %increment the index
        end
        %one dimensional cartesian correlations, mask in the other direction (if not already done by the premask)
        %initalize the mask
        one_d_mask_pos=true(size(delta,1),1);
        one_d_mask_neg=one_d_mask_pos;
        %run the mask in all the queued dimensions
        for mask_dim=mask_dimensions
            one_d_mask_pos=one_d_mask_pos & delta(:,mask_dim)>corr_opts.one_d_window(mask_dim,1) ...
                & delta(:,mask_dim)<corr_opts.one_d_window(mask_dim,2);
            one_d_mask_neg=one_d_mask_neg & -delta(:,mask_dim)>corr_opts.one_d_window(mask_dim,1) ...
                & -delta(:,mask_dim)<corr_opts.one_d_window(mask_dim,2);
        end
        %to be strictly accurate we must calaulate things symetricaly
        one_d_bins(shotnum,:)=one_d_bins(shotnum,:)+histcounts(delta(one_d_mask_pos,corr_opts.one_d_dimension),corr_opts.one_d_edges);
        one_d_bins(shotnum,:)=one_d_bins(shotnum,:)+histcounts(-delta(one_d_mask_neg,corr_opts.one_d_dimension),corr_opts.one_d_edges);
        if mod(shotnum,update_interval)==0
            parfor_progress_imp;
        end
    end%loop over shots
end %done calculating with either high or low mem
    
    
   
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


% %%you can check that the sorted mask this is equivelent to the below brute approach
% delta_brute=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+1:num_counts_shot,:); 
% min_lim=corr_opts.one_d_window(corr_opts.sorted_dir,1);
% max_lim=corr_opts.one_d_window(corr_opts.sorted_dir,2);
% mask=delta_brute(:,corr_opts.sorted_dir)>min_lim & delta_brute(:,corr_opts.sorted_dir)<max_lim;
% delta_brute =delta_brute(mask,:);
% if ~isequal(delta_brute,delta)
%     fprintf('warning mask is not giving correct answer')
% end