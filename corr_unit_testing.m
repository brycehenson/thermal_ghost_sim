%corr unit testing
%first fake dataset
fake_data=fake_correlations_one;

%%
sfigure(1);
clf


corr_opts=[];
corr_opts.window=[1,1,1]*1e-1;
corr_opts.redges=linspace(0,0.5,40);
corr_opts.xedges=linspace(-0.5,0.5,40);

fprintf('calculating in shot correlations \n')
shotscorr=CorrRadX(corr_opts,fake_data.counts_txy);


set(gcf,'color','w');
subplot(2,3,1)
plot(shotscorr.rad_centers,shotscorr.rad_bins.*shotscorr.rad_centers.^-2/shotscorr.pairs)
title('In Shot Rad Dist')
xlabel('Radial Seperation')
subplot(2,3,4)
plot(shotscorr.x_centers,shotscorr.x_bins/shotscorr.pairs)
title('In Shot X Dist (windowed)')
xlabel('X Seperation')




%%
%then calculate the normalization
approx_pairs=CountUpperTriangle(total_counts);
fprintf('rough total pairs in norm     %2.3e \n',approx_pairs)
accurate_pairs=0;
for n=1:number_shots
    accurate_pairs=accurate_pairs+(total_counts-counts_in_shot(n))*counts_in_shot(n)*0.5;
end
fprintf('actual total pairs in norm    %2.3e \n',accurate_pairs)
fprintf('ratio acc./approx             %2.3e \n',accurate_pairs/approx_pairs)


all_counts=[counts_txy{:}];
all_counts=all_counts(:,randperm(total_counts));
n=1;
i=1;
counts_chunked={};
while n<total_counts
    max_index=min([total_counts,n+norm_chunk_size-1]);
    counts_chunked{i}=all_counts(:,n:max_index);
    i=i+1;
    n=n+norm_chunk_size+1;
end

fprintf('calculating inter-shot correlations \n ')
normcorr=CorrRadX(corr,counts_chunked);

fprintf('sampled frac of norm pairs    %3.3e\n',normcorr.pairs/accurate_pairs)
pair_sample_factor=(normcorr.pairs/shotscorr.pairs);
fprintf('sampled norm/corr pairs       %2.3e \n',pair_sample_factor)


subplot(2,3,2)
plot(normcorr.rad_centers,normcorr.rad_bins.*normcorr.rad_centers.^-2/normcorr.pairs)
title('Between Shot Rad Dist')
subplot(2,3,5)
plot(normcorr.x_centers,normcorr.x_bins/normcorr.pairs)
title('Between Shot X Dist (windowed)')
subplot(2,3,6)
xg2=pair_sample_factor*shotscorr.x_bins./normcorr.x_bins;
plot(shotscorr.x_centers,xg2)
title('g2 X')
subplot(2,3,3)
radg2=pair_sample_factor*shotscorr.rad_bins./normcorr.rad_bins;
plot(shotscorr.rad_centers,radg2)
title('g2 Rad')
pause(0.1);
saveas(gcf,'Correlations.png')
saveas(gcf,'Correlations.pdf')
savefig('Correlations.fig')
%saveas(gcf,['Correlations',num2str(num),'.png'])
%saveas(gcf,['Correlations',num2str(num),'.pdf'])
fprintf('max g2 rad                    %4.2f \n',max(radg2))
fprintf('       x                      %4.2f \n',max(xg2))

toc;
fprintf('\n\n\n\n')