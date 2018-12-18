function counts_chunked=chunk_data(data,norm_samp_factor,sort_dir)
    %data chunking function for normalization
    %see calc_any_g2_type for usage
    
    %input
    %norm_samp_factor-should be about the correlation amp for equal noise contibution from in-shot & between shot
    %data.txy
    %data.data.num_counts
    %to do
    %   -documentation
    
    total_counts=sum(data.num_counts);
    norm_chunk_size=round(mean(data.num_counts)*sqrt(norm_samp_factor));
    all_counts=vertcat(data.counts_txy{:});
    all_counts=all_counts(randperm(total_counts),:);
    iimax=ceil(total_counts/norm_chunk_size);
    counts_chunked=cell(1,iimax);
    for ii=1:iimax
        min_idx=(ii-1)*norm_chunk_size+1;
        if ii==iimax
            max_idx=total_counts;
        else
            max_idx=ii*norm_chunk_size;
        end
        tmp_data=all_counts(min_idx:max_idx,:);
        if ~isnan(sort_dir)
            [~,order]=sort(tmp_data(:,sort_dir));
            tmp_data=tmp_data(order,:);
        end
        counts_chunked{ii}=tmp_data;
    end
    if size(vertcat(counts_chunked{:}),1)~=total_counts
        warning('lost counts')
    end
end

