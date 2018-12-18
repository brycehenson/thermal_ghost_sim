function out=calc_any_g2_type(corr_opts,data)
%caucluates normalized g2 functions for a lot of different cases
%uses a lot of correlators from https://github.com/spicydonkey/correlation-funs (with some modifications)
%input
%counts_txy- cell array (n_shots*[1 or 2]) of txy data (n counts * 3)
%            txy must be sorted in time for prewindowing option
%            cell array can be n_shots*2 for two port xcorr
%corr_opts.type- what type of g2 fo you want to calculate
%               '1d_cart_cl' %window out the other 2 axes and then do a 1d correlation in the remaining
%               '1d_cart_bb' %differences are cacluated as the vector sum instead of difference
%               'radial_cl'  %Euclidean distance/L2 norm of the differences
%               'radial_bb'  %differences are calculates as the vector sum
%               '3d_cart_cl' 
%               '3d_cart_bb'
%               '3d_cart_cl'

%output
%norm g2 amplitude
%norm g2 vector with cordinates

%improvements
%- fit g2 amplitude
%- implement all the options
%- think more carefully about how the coicidences should be normalized
%  - results should be invariant under
%  - QE change
%  - Change in how the shot is cut up
% - implement bootstraping
%   - most usefull at the shot level (more practical to implement & more informative)
% - different smoothing for G2(in-shot) G2(between-shot)



% set the number of updates in a smart way
% should input check everything used here
if ~isfield(corr_opts,'progress_updates') || isnan(corr_opts.progress_updates)
    update_time=2;
    pairs_per_sec=5e7*(1+corr_opts.do_pre_mask*10);
    dyn_updates=round(update_time*size(data.counts_txy,2)*...
                               mean((corr_opts.attenuate_counts*data.num_counts).^2)/(pairs_per_sec));
    corr_opts.progress_updates=min(100,max(5,dyn_updates));
end

if isfield(corr_opts,'norm_samp_factor')
    if corr_opts.norm_samp_factor<0.01 || corr_opts.norm_samp_factor>200
        error('corr_opts.norm_samp_factor exceeds limits');
    end
else
    corr_opts.norm_samp_factor=3;
end

if ~isfield(corr_opts,'sort_norm')
    corr_opts.sort_norm=false;
end

if isequal(corr_opts.type,'1d_cart_cl')  || isequal(corr_opts.type,'1d_cart_bb') 
    
    if isequal(corr_opts.type,'1d_cart_cl') 
        corr_opts.cl_or_bb=false;
    elseif isequal(corr_opts.type,'1d_cart_bb') 
          corr_opts.cl_or_bb=true;
    end
    shotscorr=corr_1d_cart(corr_opts,data.counts_txy);
    %shotscorr_high=corr_1d_cart_high_mem(corr_opts,data.counts_txy); 
    if corr_opts.plots
        sfigure(1);
        clf
        set(gcf,'color','w');
        subplot(1,3,1)
        plot(shotscorr.x_centers,shotscorr.one_d_corr_density,'.k-','MarkerSize',10)
        title('In Shot X Dist (windowed)')
        ylabel('G^2 coincedence density')
        xlabel('X Seperation')
        pause(1e-6);
    end
    norm_sort_dir=corr_opts.sorted_dir;
    if ~corr_opts.sort_norm,norm_sort_dir=nan; end
    counts_chunked=chunk_data(data,corr_opts.norm_samp_factor,norm_sort_dir);
    corr_opts.do_pre_mask=corr_opts.sort_norm; %can only do premask if data is sorted
    fprintf('calculating inter-shot correlations \n')
    normcorr=corr_1d_cart(corr_opts,counts_chunked);
    if corr_opts.plots
        subplot(1,3,2)
        plot(normcorr.x_centers,normcorr.one_d_corr_density,'.k-','MarkerSize',10)
        title('Between Shot X Dist (windowed)')
        ylabel('G^2 coincedence density')
        xlabel('X Seperation')
        subplot(1,3,3)
        xg2=shotscorr.one_d_corr_density./normcorr.one_d_corr_density;
        plot(shotscorr.x_centers,xg2,'.k-','MarkerSize',10)
        title('Norm. Corr.')
        ylabel('g^{(2)} (X)')
        xlabel('X Seperation')
        pause(1e-6);
    end
    fprintf('g2 peak amplitude         %4.2f \n',max(xg2))
    out.in_shot_corr.x_centers=shotscorr.x_centers;
    out.in_shot_corr.one_d_corr_density=shotscorr.one_d_corr_density;
    out.between_shot_corr.x_centers=normcorr.x_centers;
    out.between_shot_corr.one_d_corr_density=normcorr.one_d_corr_density;
    out.norm_g2.x_centers=shotscorr.x_centers;
    out.norm_g2.g2_amp=xg2;
        

elseif isequal(corr_opts.type,'radial_cl')  || isequal(corr_opts.type,'radial_bb')
    if isequal(corr_opts.type,'radial_cl') 
        corr_opts.cl_or_bb=false;
    elseif isequal(corr_opts.type,'radial_bb') 
          corr_opts.cl_or_bb=true;
    end
    shotscorr=corr_radial(corr_opts,data.counts_txy);

    if corr_opts.plots
        sfigure(1);
        clf
        set(gcf,'color','w');
        subplot(1,3,1)
        plot(shotscorr.rad_centers,shotscorr.rad_corr_density,'.k-','MarkerSize',10)
        title('In Shot X Dist (windowed)')
        ylabel('G^2 coincedence density')
        xlabel('X Seperation')
        pause(1e-6);
    end
    norm_sort_dir=corr_opts.sorted_dir;
    if ~corr_opts.sort_norm,norm_sort_dir=nan; end
    
    counts_chunked=chunk_data(data,corr_opts.norm_samp_factor,norm_sort_dir);
    corr_opts.do_pre_mask=corr_opts.sort_norm; %can only do premask if data is sorted
    
    fprintf('calculating inter-shot correlations \n')
    normcorr=corr_radial(corr_opts,counts_chunked);
    if corr_opts.plots
        subplot(1,3,2)
        plot(normcorr.rad_centers,normcorr.rad_corr_density,'.k-','MarkerSize',10)
        title('Between Shot X Dist (windowed)')
        ylabel('G^2 coincedence density')
        xlabel('X Seperation')
        subplot(1,3,3)
        xg2=shotscorr.rad_corr_density./normcorr.rad_corr_density;
        plot(shotscorr.rad_centers,xg2,'.k-','MarkerSize',10)
        title('Norm. Corr.')
        ylabel('g^{(2)} (X)')
        xlabel('X Seperation')
        pause(1e-6);
    end
    fprintf('g2 peak amplitude         %4.2f \n',max(xg2))
    out.in_shot_corr.rad_centers=shotscorr.rad_centers;
    out.in_shot_corr.one_d_corr_density=shotscorr.rad_corr_density;
    out.between_shot_corr.rad_centers=normcorr.rad_centers;
    out.between_shot_corr.one_d_corr_density=normcorr.rad_corr_density;
    out.norm_g2.rad_centers=shotscorr.rad_centers;
    out.norm_g2.g2_amp=xg2;
    
end


end

