function windowing_function=text_to_window_fun(input,verbose)
%inputs
%input domain
%text

%to do
%implement changeable domain
%smart hole detection

input.text='He*';
input.domain=[-1,1];
%develop a windowing function from some text

%%

text_image=ones(180,350);
text_image=insertText(text_image,[-50,-80],input.text,'FontSize',200,'BoxOpacity',0,'Font','Arial' );
text_image=imresize(text_image,10,'bicubic');
text_image = imbinarize(text_image);
text_image= squeeze(text_image(:,:,1));
text_image=imcomplement(text_image);

[B,L] = bwboundaries(text_image); %why does the output not label holes!

if verbose>1
    figure(3)
    set(gcf,'color','w')
    clf
    subplot(1,2,1)
    imagesc(text_image)
    colormap('gray')
    %%
    subplot(1,2,2)
    imshow(label2rgb(L, @jet, [.5 .5 .5]))
    hold on
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
    end


    %%
    figure(4)
    subplot(1,3,1)
    boundary = B{4};
    plot(polyshape(boundary(:,2), boundary(:,1)))
    %plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2)


    subplot(1,3,2)
    poly_simple=[];
    [poly_simple(:,2),poly_simple(:,1)]=reducem(boundary(:,2), boundary(:,1),5);
    plot(polyshape(poly_simple(:,2), poly_simple(:,1)))
end



%% Make a cell array of the simplified polydons and scale to the unit interval
text_segments= {};
text_segments.holes=logical([0,0,0,1]);
text_segments.polygons={};
image_size=size(text_image);
image_size(1)=-image_size(1);
image_shift=[0.46,-0.5];
image_rescale=[2,2];
for ii=1:size(B,1)
    boundary = B{ii};
    poly_simple=[];
    [poly_simple(:,2),poly_simple(:,1)]=reducem(boundary(:,2), boundary(:,1),5);
    %plot(polyshape(poly_simple(:,2), poly_simple(:,1)))
    poly_simple=poly_simple./repmat(image_size,[size(poly_simple,1),1]);
     poly_simple=poly_simple+repmat(image_shift,[size(poly_simple,1),1]);
     poly_simple=poly_simple.*repmat(image_rescale,[size(poly_simple,1),1]);
    text_segments.polygons{ii}=poly_simple;
end


windowing_function=@(txy_in) window_poly_fun(txy_in,text_segments);

if verbose>1
    figure(5)
    clf
    hold on
    for ii=1:numel(text_segments.polygons)
        plot(polyshape(text_segments.polygons{ii}(:,2), text_segments.polygons{ii}(:,1)))
    end
    hold off


    %% Test the method
    txy_rand=(rand(1e6,3)-0.5)*2;
    
    mask=windowing_function(txy_rand) ;%& ~mask_in_holes(2,:);
    plot(txy_rand(mask,2),txy_rand(mask,3),'.')
end

end

function mask=window_poly_fun(txy_in,text_segments)
    mask_in_holes=zeros(size(text_segments.holes,2),size(txy_in,1));
    for kk=1:numel(text_segments.polygons)
        polygons=text_segments.polygons{kk};
        mask_in_holes(kk,:)=inpolygon(txy_in(:,2),txy_in(:,3),polygons(:,2),polygons(:,1));
    end
    %find the case that the count is in any of the segments but in none of the holes
    mask=any(mask_in_holes(~text_segments.holes,:),1) & ~any(mask_in_holes(text_segments.holes,:),1) ;%& ~mask_in_holes(2,:);
end


