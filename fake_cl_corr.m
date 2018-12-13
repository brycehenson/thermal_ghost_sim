function fake_data=fake_cl_corr(number_shots,qe)
% test function on that makes some colinearly correlated data
%return data in the same structure that import_mcp_tdc_data produces

shot_counts=1000; %number of counts in each shot
corr_frac=1;
therm_width=1; %isotropic in 3d
corr_widthx=0.1;
corr_widthy=0.1;
corr_widthz=0.1;

fake_data=[];
num_corrected=round(shot_counts/(1+corr_frac));
for n=1:number_shots
    shot=reshape(randn(num_corrected*3,1)*therm_width,[num_corrected,3]);
    coor_chance=rand(num_corrected,1)<corr_frac;
    if sum(coor_chance)>0
        delpos=reshape(randn(sum(coor_chance)*3,1),[sum(coor_chance),3]).*repmat([corr_widthz,corr_widthx,corr_widthy],[sum(coor_chance),1]);
        shot_with_corr=[shot;shot(coor_chance,:)+delpos];
    else
        shot_with_corr=shot;
    end
    
    det_chance=rand(size(shot_with_corr,1),1)<qe;
    num_det=sum(det_chance);
    if  num_det>1
        det_counts=shot_with_corr(det_chance,:);
        %sort in time to make sure everything is robust to that
        [~,sort_idx]=sort(det_counts(:,1));
        det_counts=det_counts(sort_idx,:);
        counts_txy{n}=det_counts;
    else
        counts_txy{n}=[];
    end
    num_counts(n)=num_det;
end

total_counts=size([counts_txy{:}],1);
fprintf('total counts over all shots   %2.3e \n',total_counts)
fprintf('mean per shot                 %2.3e \n',total_counts/number_shots)
fprintf('rough pairs per shot          %2.3e \n',CountUpperTriangle(total_counts/number_shots))
fprintf('rough total pairs             %2.3e \n',CountUpperTriangle(total_counts/number_shots)*number_shots)

fake_data.counts_txy=counts_txy;
fake_data.num_counts=num_counts;

end