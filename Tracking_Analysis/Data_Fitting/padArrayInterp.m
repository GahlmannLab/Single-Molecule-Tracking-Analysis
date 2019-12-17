function [xf yf] = padArrayInterp(xi,yi,xq)
% Julian Rocha

%This function pads the data so that the interpolation can sample the full
%range appropriately.   

    maxX = max(xq);
    padStep = xq(2)/5;
    padSize1 = floor(min(xi)/padStep);
    padStart1 = padSize1*padStep;
    
    yi = padarray(yi,padSize1,0,'pre');
    for j = 1:padSize1
    xi = padarray(xi,1,padStart1,'pre');
    padStart1 = padStart1 - padStep;
    end
    
    padSize2 = floor((maxX-max(xi))/padStep);
    padStart2 = maxX-padSize2*padStep;
    if padSize2 < 0 
        padSize2 = 1;
        padStart2 = maxX;
    end
    
    yi = padarray(yi,padSize2,1,'post');
    for j = 1:padSize2
    xi = padarray(xi,1,padStart2,'post');
    padStart2 = padStart2 + padStep;
    end

    xf = xi;
    yf = yi;
    
end
