function fake_data=make_bloby_data(number_shots,qe)

shot_counts=2000; %number of counts in each shot
therm_width=0.5; %isotropic in 3d, standard deviation of the distribution
blob_width=1e-2;
blob_number=50;%floor(shot_counts/10);
blob_dist_width=0.25; %how wide the distribution of blob centers should be
use_gaussian_background_dist=false;

counts_txy=[];
fake_data=[];
for n=1:number_shots
    sampled_shot_counts=poissrnd(shot_counts);
    %sampled_shot_counts=shot_counts;
    if use_gaussian_background_dist
        shot=mvnrnd(zeros(3,1),diag([1,1,1])*therm_width,poissrnd(shot_counts));
    else
        shot=(rand(sampled_shot_counts,3)-0.5)*4*therm_width;
    end

    if blob_number~=0
        blob_cen=mvnrnd(zeros(3,1),diag([1,1,1])*blob_dist_width,1);
        blob=blob_cen+mvnrnd(zeros(3,1),diag([1,1,1])*blob_width,poissrnd(blob_number));
        shot_with_blob=[shot;blob];
    else
        shot_with_blob=shot;
    end
    
    det_chance=rand(size(shot_with_blob,1),1)<qe;
    num_det=sum(det_chance);
    if  num_det>1
        det_counts=shot_with_blob(det_chance,:);
        %sort in time to make sure everything is robust to that
        [~,sort_idx]=sort(det_counts(:,1));
        det_counts=det_counts(sort_idx,:);
        counts_txy{n}=det_counts;
    else
        counts_txy{n}={};
    end
    num_counts(n)=num_det;
end

fake_data.counts_txy=counts_txy;
fake_data.num_counts=num_counts;

end