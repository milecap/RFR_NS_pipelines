% Author: Milena Capiglioni, University of Bern
% Contact: milena.capiglioni@extern.insel.ch
% Last update: Nov.2023

function [min_s,max_s,aux_max,norm_min,norm_max] = im_overlay_RFR(Yoff,contrast,mask,min_s,max_s,aux_max,norm_min,norm_max,ax) 

gray_im = Yoff;
over_im = contrast;
alpha_mask = mask;
n_colors = 256;

%turn gray_im image to RGB
max_c = max(double(gray_im),[],'all');
c = gray(256);
RGB_ana = ind2rgb(floor(256*mat2gray(gray_im,[0 max_c+200])),c);
imshow(RGB_ana,[],'Parent',ax);
            
%turn over_im image to RGB
over_im = imgaussfilt(over_im,1);
over_im = mask.*over_im;

if isempty(aux_max)
    aux_max = max(abs(over_im),[],'all');
end
over_im = over_im./aux_max;

% if isempty(norm_min)
%     norm_min = min(over_im,[],'all');
% end
% over_im(over_im < 0) = over_im(over_im < 0)./abs(norm_min);
% 
% if isempty(norm_max)
%     norm_max = max(over_im,[],'all');
% end
% over_im(over_im > 0) = over_im(over_im > 0)./norm_max;

if isempty(max_s)
    max_s = max(over_im,[],'all');
end
if isempty(min_s)
    min_s = min(over_im,[],'all');
end

% over_im(isnan(over_im)) = 0;

[c,num,typ,scheme] = brewermap(n_colors,'RdYlBu');
c = flipud(c);

RGB_im = ind2rgb(floor(n_colors*mat2gray(over_im,[-1 1])),c); %
hold(ax,'on')
h = imshow(RGB_im,'Parent',ax);

colormap(ax,c);
colorbar(ax)
caxis(ax,[min_s max_s]);

hold(ax,'off')

set(h,'AlphaData', (~isnan(alpha_mask))) %.*(over_im<0.01)