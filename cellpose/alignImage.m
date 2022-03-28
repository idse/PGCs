function aligned = alignImage(im1, im2, shiftyx)
%after mapping points, overlay the (shifted) images and show
%the link between fixed and live points


mn = min(size(im1), size(im2));
m = mn(1); n = mn(2);
im2 = im2(1:mn(1), 1:mn(2));

shifty = shiftyx(1);
shiftx = shiftyx(2);

xinrange = max(1,1+shiftx):min(m,m+shiftx);
xoutrange = max(1,1-shiftx):min(m,m-shiftx);
yinrange = max(1,1+shifty):min(n,n+shifty);
youtrange = max(1,1-shifty):min(n,n-shifty);

aligned = zeros(size(im1),class(im2));
aligned(xoutrange, youtrange) = im2(xinrange, yinrange);

end