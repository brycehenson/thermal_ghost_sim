
arg={};
a=rand(1e4,1,arg{:}); %
tic
numbins=1e3;
edges=linspace(-1,1,numbins); %gpuArray.
l=size(a,1);
bins=ones(1,numbins,arg{:}); 
for ii=1:l
    disances=a(ii)-a(ii+1:end);
    bins=bins+histc(disances,edges);
end
toc

plot(edges,bins)



%%
a=gpuArray.rand(1e4,1); %
tic
numbins=1e3;
edges=gpuArray.linspace(-1,1,numbins); %.
bins=gpuArray.ones(1,numbins); 
for ii=1:size(a,1)
    bins=bins+histc(a(ii)-a(:),edges);
end
toc

plot(edges,bins)



%%
numbins=1+1e7;
numcounts=1e6;


a=rand(numcounts,1); %
tic
edges=linspace(0,1,numbins); %.
bins=[histc(a,edges);zeros(numbins,1)]; %
bfft_cpu=fft(bins);
xc=ifft(bfft_cpu.*bfft_cpu);
toc
plot(xc)


%%
tic
edges=linspace(0,1,numbins); %.
bins=[histc(a,edges)]; %
xc=xcorr(bins,bins,10);
toc
plot(xc)
%%

a=gpuArray.rand(numcounts,1); %
tic
edges=gpuArray.linspace(0,1,numbins); %.
bins=[histc(a,edges)]; %
xc=xcorr(bins,bins,10);
toc
plot(xc)


%%

a=gpuArray.rand(numcounts,1); %
tic
edges=gpuArray.linspace(0,1,numbins); %.
bins=[histc(a,edges);zeros(numbins,1)]; %
bfft_gpu=fft(bins);
xc=ifft(bfft_gpu.*bfft_gpu);

toc
plot(real(xc))


%yay a speedup!



%%

x_cpu=single(linspace(0,1,10));
fft_cpu=fft(x_cpu);

x_gpu=gpuArray.linspace(0,1,10)
fft_gpu=fft(x_gpu);


(fft_cpu-fft_gpu)./mean([fft_gpu;fft_cpu])


