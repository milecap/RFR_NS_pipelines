function out = slices2mosaic(mat,fillValue)

    if nargin < 2
        fillValue = 0;  % Default to filling with zeros
    end

    x=size(mat,1);
    y=size(mat,2);
    z=size(mat,3);

    %t = size(mat,4);

    %   find the side length of the mosaic image (number of images)
    m_side = ceil(sqrt(z));

    %   create the 2D mosaic matrix
    if strcmpi(fillValue, 'nan')
        out = nan(m_side*x, m_side*y);
    else
        out = zeros(m_side*x,m_side*y);
    end

    for j = 1:z
        row = (floor((j-1)/m_side)+1);
        col = j-floor((j-1)/m_side)*(m_side);
        out((row-1)*x+1:row*x,(col-1)*y+1:col*y) = mat(:,:,j);
    end     
end