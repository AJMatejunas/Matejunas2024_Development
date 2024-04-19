function [ PHI_X, PHI_Y, modX, modY ] = LSA( im, g, p)
% Local spectrum analysis (LSA)
% usage:
% 1) input:
% im: grid image (double)
% g: analysis window (cf build_window.m)
% p: grid pattern pitch (in pixels)
% 2) output
% PHI_X: phase along x-direction
% PHI_Y: phase along y-direction 
% modX: modulus of the LSA along x-direction
% modY: modulus of the LSA along y-direction

size_x=size(im,1);
size_y=size(im,2);
sizepad=size(im)+size(g)-1;
border=(size(g)-1+mod(size(g)-1,2))/2; 
fc=2*pi/p; 
[XX,YY]=meshgrid(0:size_y-1,0:size_x-1);

% x-direction
ima = im.*exp(-1i*fc*XX);
fft2ima=fft2(ima,sizepad(1),sizepad(2));
wfx = ifft2(fft2ima.*fft2(g,sizepad(1),sizepad(2)));
wfx=wfx(border(1)+1:border(1)+size_x,border(2)+1:border(2)+size_y);

modX=abs(wfx);

PHI_X = angle(wfx);

% y-direction
% the window has to be properly rotated to deal with the triangular-rectangular window
if (size(g,1)~=size(g,2))
    gg=g;
    g=zeros(max(size(g)+1));
    g(1:end-1,2:end-1)=gg;
end
ima = im.*exp(-1i*fc*YY);
fft2ima=fft2(ima,sizepad(1),sizepad(2));
wfy = ifft2(fft2ima.*fft2(g',sizepad(1),sizepad(2)));
wfy=wfy(border(1)+1:border(1)+size_x,border(2)+1:border(2)+size_y);

modY=abs(wfy);

PHI_Y=angle(wfy);


end

