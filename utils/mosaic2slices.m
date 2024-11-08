%% This function transform a 2D mosaic matrix into a 3D one with im_size,im_size,slices

function out_sl = mosaic2slices(mat,im_size, slices)
    num = (size(mat,2)/im_size);
    aux = mat(:,:);
    out = permute(reshape(aux',[size(mat,2),size(mat,2)/num,num]),[2,1,3]); %change the numbers accordingly
    out_2 = permute(reshape(out,[im_size,im_size,num^2]),[1,2,3]);
    out_sl = out_2(:,:,1:slices);
end
