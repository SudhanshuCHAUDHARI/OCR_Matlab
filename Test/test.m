I = imread('TEST_2.jpg');
Ibw = imbinarize(I);
Ibw = imfill(Ibw,'holes');
%Ilabel = bwlabel(Ibw);
stat = regionprops(~Ibw,'centroid');
imshow(I); hold on;
for x = 1: numel(stat)
      plot(stat(x).Centroid(1),stat(x).Centroid(2),'*');
end