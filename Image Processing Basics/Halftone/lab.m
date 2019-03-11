%% We read a figure, produce thresholded image
clear
clc
X=imread('house.tif');
Y= zeros(size(X));
for i =1:256
    for j = 1:384
        if X(i,j) > 127,
            Y(i,j) = 255;
        end
    end
end

imwrite(Y,'threshold.tiff','TIFF');

%% halftoned image using different sized Bayer index matrices
XX = im2double(X);
YY = im2double(Y);

RMSE = sqrt(sum(sum((XX-YY).^2,1),2)/256/384);

fid = fidelity(X,Y);

XX = im2double(X);

XXg = (XX).^2.2;


I2=[1,2;3,0];
I4 = [ 4*I2 + 1, 4*I2 + 2; 4*I2 + 3, 4*I2];
I8 = [ 4*I4 + 1, 4*I4 + 2; 4*I4 + 3, 4*I4];

I2N=repmat(I2,128,192);
I4N=repmat(I4,64,96);
I8N=repmat(I8,32,48);

T2N=255*(I2N+0.5)/256/384;
T4N=255*(I4N+0.5)/256/384;
T8N=255*(I8N+0.5)/256/384;

X2=zeros(size(X));
X4=zeros(size(X));
X8=zeros(size(X));

for i =1:256
    for j = 1:384
        if XXg(i,j) > T2N(i,j),
            X2(i,j) = 255;
        end
        if XXg(i,j) > T4N(i,j),
            X4(i,j) = 255;
        end
        if XXg(i,j) > T8N(i,j),
            X8(i,j) = 255;
        end
    end
end

imwrite(X2,'X2.tiff','TIFF');
imwrite(X4,'X4.tiff','TIFF');
imwrite(X8,'X8.tiff','TIFF');

fid2 = fidelity(X,X2);
RMSE2 = sqrt(sum(sum((XX-X2).^2,1),2)/256/384);
fid4 = fidelity(X,X4);
RMSE4 = sqrt(sum(sum((XX-X4).^2,1),2)/256/384);
fid8 = fidelity(X,X8);
RMSE8 = sqrt(sum(sum((XX-X8).^2,1),2)/256/384);


%% halftoned image using error diffusion
[d1,d2]=size(X);
X_=padarray(X,[1,1]);
X_=reshape(X_,[1,(d1+2)*(d2+2)]);
Z=zeros(size(X_));
E=zeros(size(X_));
O=zeros(size(X_));
for i =d2+4:(d1+1)*(d2+2)
    if X_(i) > 127,
       Z(i) = 255;
    else Z(i) = 0;
    end
    E(i)= X_(i)-Z(i);
    X_(i+1)= X_(i+1)...
        +1/16*E(i-387)...
        +5/16*E(i-386)...
        +3/16*E(i-385)...
        +7/16*E(i-1);    
end

O = reshape(X_,[d1+2,d2+2]);
O = O(2:257,2:385);
imwrite(O,'output.tiff','TIFF');
XX = im2double(X);
OO = im2double(O);
RMSEO = sqrt(sum(sum((XX-OO).^2,1),2)/256/384);
fid_ = fidelity(X,O);