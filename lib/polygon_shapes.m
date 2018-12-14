
pgon = polyshape([0  2 2],[2 0 0 2]);
plot(pgon)


%%
figure(3)
set(gcf,'color','w')
clf
subplot(1,2,1)
text_image=ones(180,350);
text_image=insertText(text_image,[-50,-80],'He*','FontSize',200,'BoxOpacity',0,'Font','Arial' );
text_image = imbinarize(text_image);
text_image= squeeze(text_image(:,:,1));
text_image=imcomplement(text_image);
imagesc(text_image)
colormap('gray')
%%
subplot(1,2,2)

[B,L] = bwboundaries(text_image,'noholes');
imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

%%
figure(4)
boundary = B{3};
plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2)



%%
%   [[82.4288  109.5379];
%   [81.4925  252.9781];
%   [102.0918  254.1443];
%   [101.1554  184.1735];
%   [172.3165  184.1735];
%   [173.2528  254.1443];
%   [192.9157  254.1443];
%   [191.9794  109.5379];
%   [169.5075  108.3717]
%   [170.4438  169.0131];
%   [100.2191  170.1793];
%   [101.1554  107.2055];