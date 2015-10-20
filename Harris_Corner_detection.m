function[R,result]=Harris_Corner_detection(imagex,w)
% To Detect Corner interest point in an image
% input -- imagex (input image)
%          w (window size)
if(ndims(imagex)==3)
image=im2double(rgb2gray(imagex));
else
image=im2double(imagex);
end

nx=size(image);

if(mod((w-1),2)~=0)
    error('size of patch window is not good enough');
end

X=zeros(nx(1)+((w-1)),nx(2)+((w-1)));
X(((w-1)/2)+1:nx(1)+((w-1)/2),((w-1)/2)+1:nx(2)+((w-1)/2))=image;
image=X;
clear X;



%using sobel operators
for i=((w-1)/2)+1:nx(1)+((w-1)/2)
for j=((w-1)/2)+1:nx(2)+((w-1)/2)
Gx(i-((w-1)/2),j-((w-1)/2))=image(i,j)+0-image(i,j+2)+2*(image(i+1,j)+0-image(i+1,j+2))+image(i+2,j)+0-image(i+2,j+2);
Gy(i-((w-1)/2),j-((w-1)/2))=image(i,j)+2*image(i,j+1)+image(i,j+2)+0+0+0-image(i+2,j)-2*image(i+2,j+1)-image(i+2,j+2);
end
end


Ix1=Gx.^2;
Ix2=Gy.^2;
Ix3=Gx.*Gy;


% Gaussian mask with sigma 1.4
gauss_mask=Gaussian_template(w,1.4);

% apply Gaussian mask on image for Smoothing 
IX2=filter2(gauss_mask,Ix1);
IY2=filter2(gauss_mask,Ix2);
IXY=filter2(gauss_mask,Ix3);

nx=size(IX2);
Rmax=0;
for i=1:nx(1)
for j=1:nx(2)
N=im2double([IX2(i,j) IXY(i,j);IXY(i,j) IY2(i,j)]);
R(i,j)=det(N)-0.01*(trace(N))^2;
if(R(i,j)>Rmax)
Rmax=R(i,j);
end
end
end


nx=size(R);
for i = 2:nx(1)-1
for j = 2:nx(2)-1
if(R(i,j)>0.1*Rmax &&R(i,j)> R(i-1,j-1) && R(i,j) > R(i-1,j) && R(i,j) > R(i-1,j+1) && R(i,j) > R(i,j-1) && R(i,j) > R(i,j+1) && R(i,j) > R(i+1,j-1) && R(i,j) > R(i+1,j) &&(( R(i,j) > R(i+1,j+1))|| R(i,j)<0))
result(i,j) = 1;
end;
end;
end;
[posc, posr] = find(result==1);
imshow(imagex);
hold on;
plot(posr,posc,'r.');
end
