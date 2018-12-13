%corr unit testing
% main questions
% should the in shot pairs be counted in the number of possible normalization pairs


%first fake dataset
fake_data=fake_cl_corr(100,1);
total_counts=sum(fake_data.num_counts);
%%
sfigure(1);
clf

corr_opts=[];
corr_opts.dimension=1;
corr_opts.window=[[-1,1];[-1,1];[-1,1]]*1e-2;
rmax=0.8;
rmin=1e-4;
%corr_opts.redges=linspace(rmin,rmax,40);
corr_opts.redges=sqrt(linspace(rmin^2,rmax^2,120));
xrange=0.8;
corr_opts.xedges=linspace(-xrange,xrange,40);
corr_opts.rad_smoothing=0.05;
corr_opts.one_d_smoothing=0.05;

fprintf('calculating in shot correlations \n')
tic
shotscorr=CorrRadX(corr_opts,fake_data.counts_txy);
toc

set(gcf,'color','w');
subplot(2,3,1)
plot(shotscorr.rad_centers,shotscorr.rad_corr_density)
title('In Shot Rad Dist')
xlabel('Radial Seperation')
ylabel('G^2 coincedence density')
subplot(2,3,4)
plot(shotscorr.x_centers,shotscorr.one_d_corr_density)
title('In Shot X Dist (windowed)')
ylabel('G^2 coincedence density')
xlabel('X Seperation')



%% chunk up the data for normalization
norm_chunk_size=mean(fake_data.num_counts)*4;

all_counts=[vertcat(fake_data.counts_txy{:})];

all_counts=all_counts(randperm(total_counts),:);
n=1;
i=1;
counts_chunked={};
while n<total_counts
    max_index=min([total_counts,n+norm_chunk_size-1]);
    counts_chunked{i}=all_counts(n:max_index,:);
    i=i+1;
    n=n+norm_chunk_size;
end
if size(vertcat(counts_chunked{:}),1)~=total_counts
    error('lost counts')
end

%%
fprintf('calculating inter-shot correlations \n ')
normcorr=CorrRadX(corr_opts,counts_chunked);

%%
subplot(2,3,2)
plot(normcorr.rad_centers,normcorr.rad_corr_density)
title('Between Shot Rad Dist')
xlabel('Radial Seperation')
ylabel('G^2 coincedence density')
subplot(2,3,5)
plot(normcorr.x_centers,normcorr.one_d_corr_density)
title('Between Shot X Dist (windowed)')
ylabel('G^2 coincedence density')
xlabel('X Seperation')
subplot(2,3,6)
xg2=shotscorr.one_d_corr_density./normcorr.one_d_corr_density;
plot(shotscorr.x_centers,xg2)
title('g2 X')
xlabel('X Seperation')
subplot(2,3,3)
radg2=shotscorr.rad_corr_density./normcorr.rad_corr_density; 
plot(shotscorr.rad_centers,radg2)
title('g2 Rad')
xlabel('Radial Seperation')
pause(0.1);
%saveas(gcf,'Correlations.png')
%saveas(gcf,'Correlations.pdf')
%savefig('Correlations.fig')
%saveas(gcf,['Correlations',num2str(num),'.png'])
%saveas(gcf,['Correlations',num2str(num),'.pdf'])
fprintf('max g2 rad                    %4.2f \n',max(radg2))
fprintf('       x                      %4.2f \n',max(xg2))

toc;
fprintf('\n\n\n\n')