%% IO with 100/4
%% Preparations
clear;
close all;
%% parameters
itermax = 50;
alpha = 1;
phasequantizer = 4;
magnitudequantizer = 100;

lena512 = imread('lena512.bmp');
X = im2double(lena512);
i = sqrt(-1);
Sample = 1.2;
loose = 0;
[imn1, imn2] = size(X);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image1   = zeros(imn11,imn22);   % image with zero paddings
start1   = round(1 + (imn11-imn1)/2);   % position of the true image
end1     = start1 + imn1 - 1;
start2   = round(1 + (imn22-imn2)/2);
end2     = start2 + imn2 - 1;
image1(start1:end2,start2:end2) = X;

%% build quandrant put into first quardrant
Sample = 2;
loose = 0;
[imn1, imn2] = size(image1);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image   = zeros(imn11,imn22);                 % image with zero paddings
image(start1:end2,start2:end2) = X;

imagetype = -2;
image     = image.*exp(i*unifrnd(0,2*pi,imn11,imn22));
% complex images have random phases in [0,2*pi]


S1start   = start1 - loose;% In the case where the image has a loose support,
S1end     = end2   + loose;% only a loose position of the true image is given
S2start   = start2 - loose;
S2end     = end2   + loose;

Index  = [imn1 , imn2 ; imn11 , imn22 ; S1start , S1end; S2start , S2end];
% data to input

X0  = generatestart(imagetype,Index);

%% IO  
N11          = Index(2,1); %size after zeropadding
N22          = Index(2,2); %size after zeropadding
N1start   = start1; %index where the true image start
N1end    = end2;
N2start  = start2;
N2end    = end2;

%% initiation
Xk           = X0;      % Xk
XkP          = X0;      % Pf{Xk}
iterres      = [];      % record residual at each iteration
iterchange   = [];      % record ||XkP{k}-XkP{k-1}||/||XkP{k-1}||

Y_abs = abs(fft2(image));

%% IO iteration
for iter = 1 : itermax  %iterate over IO
    TemMatrix = XkP;     % record XkP input from the previous iteration
    %% compute Pf{Xk}
    XkP   = Xk; %new input modified
    FXk  = ifft2(XkP);
    FXkP = Y_abs.*exp(i*angle(FXk)); %*(FXk./abs(FXk));
    FXkP_abs= floor(magnitudequantizer*abs(FXkP)./(max(max(abs(FXkP)))))...
      ./magnitudequantizer*max(max(abs(FXkP)));
    FXkP_angle = floor(phasequantizer*angle(FXkP)./(2*pi))...
      ./phasequantizer*2*pi;
    FXkP_quant= FXkP_abs.*exp(i.*FXkP_angle);
    XkP  = fft2(FXkP_quant); %new output
    %% calculate delta_g
    termA =  abs(image).*(XkP./abs(XkP)); termA(isnan(termA)) = 0;
    termB =  abs(image).*(XkP./abs(XkP)); termB(isnan(termA)) = 0;
    delta_g = termA -XkP+termA-termB;
    %% apply IO step
    new     = zeros(N11,N22);
    %% apply IO on zero paddings
    new(:,1:N2start-1)              = Xk(:,1:N2start-1) + ...
    alpha * delta_g(:,1:N2start-1);
    new(:,N2end+1:N22)              = Xk(:,N2end+1:N22) + ...
    alpha * delta_g(:,N2end+1:N22);
    new(1:N1start-1,N2start:N2end)  = Xk(1:N1start-1,N2start:N2end)...
    + alpha * delta_g(1:N1start-1,N2start:N2end);
    new(N1end+1:N11,N2start:N2end)  = Xk(N1end+1:N11,N2start:N2end)...
    + alpha * delta_g(N1end+1:N11,N2start:N2end);
    %% apply IO in support
    if imagetype == -2
        new(N1start:N1end,N2start:N2end) = ...
            Xk(N1start:N1end,N2start:N2end)+ ...
            alpha * delta_g(N1start:N1end,N2start:N2end);
        Xk             = new;
    end
    %% record residual
    Xtem = zeros(N11,N22);
    Xtem(N1start:N1end,N2start:N2end) = ...
    projection(XkP(N1start:N1end,N2start:N2end),imagetype);
    intensityerr = norm(abs(fft2(Xtem))-Y_abs,'fro')^2;
    intensityerr = sqrt(intensityerr)/norm(Y_abs,'fro');
    iterres(iter)  = intensityerr;
    %% record the change of XkP
    iterchange(iter) = norm(TemMatrix-XkP,'fro')/norm(TemMatrix,'fro');
    %% display data
    fprintf('%6.0f residual = %6.4f p change = %6.6f p \n',iter,...
        100*intensityerr,100*iterchange(iter))
    Xrec   = XkP;
    figure(1)
    subplot(1,2,1)
    imagesc(abs(Xrec(start1:end1,start2:end2)));
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram reconstructed image')
    subplot(1,2,2)
    imagesc(abs(FXkP_quant))
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram')
    drawnow;
end



%% plot the magnitudes of the recovered image

figure(2)
subplot(1,2,1)
imagesc(abs(Xrec));
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title('Hologram reconstructed image')
subplot(1,2,2)
imagesc(abs(FXkP_quant))
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title(['Hologram with ' num2str(magnitudequantizer)...
    ' discrete magnitude and ' num2str(phasequantizer) ' discrete phase'])
savefig(['h_full_' num2str(itermax) 'IO_'... 
    num2str(magnitudequantizer) '_' num2str(phasequantizer) '.fig'])

iterres_100_4_IO = iterres;
save(['h_' num2str(itermax) 'IO_' num2str(magnitudequantizer)...
    '_' num2str(phasequantizer) '_residue.mat'],'iterres_100_4_IO')

%% IO First with 256/10
%% Preparations
clear;
close all;
%% parameters
itermax = 50;
alpha = 1;
phasequantizer = 10;
magnitudequantizer = 256;

lena512 = imread('lena512.bmp');
X = im2double(lena512);
i = sqrt(-1);
Sample = 1.2;
loose = 0;
[imn1, imn2] = size(X);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image1   = zeros(imn11,imn22);   % image with zero paddings
start1   = round(1 + (imn11-imn1)/2);   % position of the true image
end1     = start1 + imn1 - 1;
start2   = round(1 + (imn22-imn2)/2);
end2     = start2 + imn2 - 1;
image1(start1:end2,start2:end2) = X;

%% build quandrant put into first quardrant
Sample = 2;
loose = 0;
[imn1, imn2] = size(image1);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image   = zeros(imn11,imn22);                 % image with zero paddings
image(start1:end2,start2:end2) = X;

imagetype = -2;
image     = image.*exp(i*unifrnd(0,2*pi,imn11,imn22));
% complex images have random phases in [0,2*pi]


S1start   = start1 - loose;% In the case where the image has a loose support,
S1end     = end2   + loose;% only a loose position of the true image is given
S2start   = start2 - loose;
S2end     = end2   + loose;

Index  = [imn1 , imn2 ; imn11 , imn22 ; S1start , S1end; S2start , S2end];
% data to input

X0  = generatestart(imagetype,Index);

%% IO  
N11          = Index(2,1); %size after zeropadding
N22          = Index(2,2); %size after zeropadding
N1start   = start1; %index where the true image start
N1end    = end2;
N2start  = start2;
N2end    = end2;

%% initiation
Xk           = X0;      % Xk
XkP          = X0;      % Pf{Xk}
iterres      = [];      % record residual at each iteration
iterchange   = [];      % record ||XkP{k}-XkP{k-1}||/||XkP{k-1}||

Y_abs = abs(fft2(image));

%% IO iteration
for iter = 1 : itermax  %iterate over IO
    TemMatrix = XkP;     % record XkP input from the previous iteration
    %% compute Pf{Xk}
    XkP   = Xk; %new input modified
    FXk  = fft2(XkP);
    FXkP = Y_abs.*exp(i*angle(FXk)); %*(FXk./abs(FXk));
    FXkP_abs= floor(magnitudequantizer*abs(FXkP)./(max(max(abs(FXkP)))))...
      ./magnitudequantizer*max(max(abs(FXkP)));
    FXkP_angle = floor(phasequantizer*angle(FXkP)./(2*pi))...
      ./phasequantizer*2*pi;
    FXkP_quant= FXkP_abs.*exp(i.*FXkP_angle);
    XkP  = ifft2(FXkP_quant); %new output
    %% calculate delta_g
    termA =  abs(image).*(XkP./abs(XkP)); termA(isnan(termA)) = 0;
    termB =  abs(image).*(XkP./abs(XkP)); termB(isnan(termA)) = 0;
    delta_g = termA -XkP+termA-termB;
    %% apply IO step
    new     = zeros(N11,N22);
    %% apply IO on zero paddings
    new(:,1:N2start-1)              = Xk(:,1:N2start-1) + ...
    alpha * delta_g(:,1:N2start-1);
    new(:,N2end+1:N22)              = Xk(:,N2end+1:N22) + ...
    alpha * delta_g(:,N2end+1:N22);
    new(1:N1start-1,N2start:N2end)  = Xk(1:N1start-1,N2start:N2end)...
    + alpha * delta_g(1:N1start-1,N2start:N2end);
    new(N1end+1:N11,N2start:N2end)  = Xk(N1end+1:N11,N2start:N2end)...
    + alpha * delta_g(N1end+1:N11,N2start:N2end);
    %% apply IO in support
    if imagetype == -2
        new(N1start:N1end,N2start:N2end) = ...
            Xk(N1start:N1end,N2start:N2end)+ ...
            alpha * delta_g(N1start:N1end,N2start:N2end);
        Xk             = new;
    end
    %% record residual
    Xtem = zeros(N11,N22);
    Xtem(N1start:N1end,N2start:N2end) = ...
    projection(XkP(N1start:N1end,N2start:N2end),imagetype);
    intensityerr = norm(abs(fft2(Xtem))-Y_abs,'fro')^2;
    intensityerr = sqrt(intensityerr)/norm(Y_abs,'fro');
    iterres(iter)  = intensityerr;
    %% record the change of XkP
    iterchange(iter) = norm(TemMatrix-XkP,'fro')/norm(TemMatrix,'fro');
    %% display data
    fprintf('%6.0f residual = %6.4f p change = %6.6f p \n',iter,...
    100*intensityerr,100*iterchange(iter))
    Xrec   = XkP;
    figure(1)
    subplot(1,2,1)
    imagesc(abs(Xrec(start1:end1,start2:end2)));
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram reconstructed image')
    subplot(1,2,2)
    imagesc(abs(FXkP_quant))
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram')
    drawnow;
end

%% plot the magnitudes of the recovered image

figure(2)
subplot(1,2,1)
imagesc(abs(Xrec));
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title('Hologram reconstructed image')
subplot(1,2,2)
imagesc(abs(FXkP_quant))
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title(['Hologram with ' num2str(magnitudequantizer)...
    ' discrete magnitude and ' num2str(phasequantizer) ' discrete phase'])
savefig(['h_full_' num2str(itermax) 'IO_'... 
    num2str(magnitudequantizer) '_' num2str(phasequantizer) '.fig'])

iterres_256_10_IO = iterres;
save(['h_' num2str(itermax) 'IO_' num2str(magnitudequantizer)...
    '_' num2str(phasequantizer) '_residue.mat'],'iterres_256_10_IO')

%% OO with 100/4
%% Preparations
clear;
close all;
%% parameters
itermax = 50;
alpha = 1;
phasequantizer = 4;
magnitudequantizer = 100;

lena512 = imread('lena512.bmp');
X = im2double(lena512);
i = sqrt(-1);
Sample = 1.2;
loose = 0;
[imn1, imn2] = size(X);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image1   = zeros(imn11,imn22);   % image with zero paddings
start1   = round(1 + (imn11-imn1)/2);   % position of the true image
end1     = start1 + imn1 - 1;
start2   = round(1 + (imn22-imn2)/2);
end2     = start2 + imn2 - 1;
image1(start1:end2,start2:end2) = X;

%% build quandrant put into first quardrant
Sample = 2;
loose = 0;
[imn1, imn2] = size(image1);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image   = zeros(imn11,imn22);                 % image with zero paddings
image(start1:end2,start2:end2) = X;

imagetype = -2;
image     = image.*exp(i*unifrnd(0,2*pi,imn11,imn22));
% complex images have random phases in [0,2*pi]


S1start   = start1 - loose;% In the case where the image has a loose support,
S1end     = end2   + loose;% only a loose position of the true image is given
S2start   = start2 - loose;
S2end     = end2   + loose;

Index  = [imn1 , imn2 ; imn11 , imn22 ; S1start , S1end; S2start , S2end];
% data to input

X0  = generatestart(imagetype,Index);

%% OO  
N11          = Index(2,1); %size after zeropadding
N22          = Index(2,2); %size after zeropadding
N1start   = start1; %index where the true image start
N1end    = end2;
N2start  = start2;
N2end    = end2;

%% initiation
Xk           = X0;      % Xk
XkP          = X0;      % Pf{Xk}
iterres      = [];      % record residual at each iteration
iterchange   = [];      % record ||XkP{k}-XkP{k-1}||/||XkP{k-1}||

Y_abs = abs(fft2(image));


%% OO iteration
for iter = 1 : itermax  %iterate over OO
    TemMatrix = XkP;     % record XkP input from the previous iteration
    %% compute Pf{Xk}
    XkP   = Xk; %new input modified
    FXk  = fft2(XkP);
    FXkP = Y_abs.*exp(i*angle(FXk)); %*(FXk./abs(FXk));
    FXkP_abs= floor(magnitudequantizer*abs(FXkP)./(max(max(abs(FXkP)))))...
      ./magnitudequantizer*max(max(abs(FXkP)));
    FXkP_angle = floor(phasequantizer*angle(FXkP)./(2*pi))...
      ./phasequantizer*2*pi;
    FXkP_quant= FXkP_abs.*exp(i.*FXkP_angle);
    XkP  = ifft2(FXkP_quant); %new output
    %% calculate delta_g
    termA =  abs(image).*(XkP./abs(XkP)); termA(isnan(termA)) = 0;
    termB =  abs(image).*(XkP./abs(XkP)); termB(isnan(termA)) = 0;
    delta_g = termA -XkP+termA-termB;
    %% apply IO step
    new     = zeros(N11,N22);
    %% apply IO on zero paddings
    new(:,1:N2start-1)              = XkP(:,1:N2start-1) + ...
    alpha * delta_g(:,1:N2start-1);
    new(:,N2end+1:N22)              = XkP(:,N2end+1:N22) + ...
    alpha * delta_g(:,N2end+1:N22);
    new(1:N1start-1,N2start:N2end)  = XkP(1:N1start-1,N2start:N2end)...
    + alpha * delta_g(1:N1start-1,N2start:N2end);
    new(N1end+1:N11,N2start:N2end)  = XkP(N1end+1:N11,N2start:N2end)...
    + alpha * delta_g(N1end+1:N11,N2start:N2end);
    %% apply IO in support
    if imagetype == -2
        new(N1start:N1end,N2start:N2end) = ...
            XkP(N1start:N1end,N2start:N2end)+ ...
            alpha * delta_g(N1start:N1end,N2start:N2end);
        Xk             = new;
    end
    %% record residual
    Xtem = zeros(N11,N22);
    Xtem(N1start:N1end,N2start:N2end) = ...
    projection(XkP(N1start:N1end,N2start:N2end),imagetype);
    intensityerr = norm(abs(fft2(Xtem))-Y_abs,'fro')^2;
    intensityerr = sqrt(intensityerr)/norm(Y_abs,'fro');
    iterres(iter)  = intensityerr;
    %% record the change of XkP
    iterchange(iter) = norm(TemMatrix-XkP,'fro')/norm(TemMatrix,'fro');
    %% display data
    fprintf('%6.0f residual = %6.4f p change = %6.6f p \n',iter,...
        100*intensityerr,100*iterchange(iter))
    Xrec   = XkP;
    figure(1)
    subplot(1,2,1)
    imagesc(abs(Xrec(start1:end1,start2:end2)));
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram reconstructed image')
    subplot(1,2,2)
    imagesc(abs(FXkP_quant))
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram')
    drawnow;
end


%% plot the magnitudes of the recovered image
figure(2)
subplot(1,2,1)
imagesc(abs(Xrec));
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title('Hologram reconstructed image')
subplot(1,2,2)
imagesc(abs(FXkP_quant))
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title(['Hologram with ' num2str(magnitudequantizer)...
    ' discrete magnitude and ' num2str(phasequantizer) ' discrete phase'])
savefig(['h_full_' num2str(itermax) 'OO_'... 
    num2str(magnitudequantizer) '_' num2str(phasequantizer) '.fig'])

iterres_100_4_OO = iterres;
save(['h_' num2str(itermax) 'OO_' num2str(magnitudequantizer)...
    '_' num2str(phasequantizer) '_residue.mat'],'iterres_100_4_OO')

%% OO with 256/10
%% Preparations
clear;
close all;
%% parameters
itermax = 50;
alpha = 1;
phasequantizer = 10;
magnitudequantizer = 256;

lena512 = imread('lena512.bmp');
X = im2double(lena512);
i = sqrt(-1);
Sample = 1.2;
loose = 0;
[imn1, imn2] = size(X);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image1   = zeros(imn11,imn22);   % image with zero paddings
start1   = round(1 + (imn11-imn1)/2);   % position of the true image
end1     = start1 + imn1 - 1;
start2   = round(1 + (imn22-imn2)/2);
end2     = start2 + imn2 - 1;
image1(start1:end2,start2:end2) = X;

%% build quandrant put into first quardrant
Sample = 2;
loose = 0;
[imn1, imn2] = size(image1);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image   = zeros(imn11,imn22);                 % image with zero paddings
image(start1:end2,start2:end2) = X;

imagetype = -2;
image     = image.*exp(i*unifrnd(0,2*pi,imn11,imn22));
% complex images have random phases in [0,2*pi]


S1start   = start1 - loose;% In the case where the image has a loose support,
S1end     = end2   + loose;% only a loose position of the true image is given
S2start   = start2 - loose;
S2end     = end2   + loose;

Index  = [imn1 , imn2 ; imn11 , imn22 ; S1start , S1end; S2start , S2end];
% data to input

X0  = generatestart(imagetype,Index);

%% OO  
N11          = Index(2,1); %size after zeropadding
N22          = Index(2,2); %size after zeropadding
N1start   = start1; %index where the true image start
N1end    = end2;
N2start  = start2;
N2end    = end2;

%% initiation
Xk           = X0;      % Xk
XkP          = X0;      % Pf{Xk}
iterres      = [];      % record residual at each iteration
iterchange   = [];      % record ||XkP{k}-XkP{k-1}||/||XkP{k-1}||

Y_abs = abs(fft2(image));

%% OO iteration
for iter = 1 : itermax  %iterate over OO
    TemMatrix = XkP;     % record XkP input from the previous iteration
    %% compute Pf{Xk}
    XkP   = Xk; %new input modified
    FXk  = fft2(XkP);
    FXkP = Y_abs.*exp(i*angle(FXk)); %*(FXk./abs(FXk));
    FXkP_abs= floor(magnitudequantizer*abs(FXkP)./(max(max(abs(FXkP)))))...
      ./magnitudequantizer*max(max(abs(FXkP)));
    FXkP_angle = floor(phasequantizer*angle(FXkP)./(2*pi))...
      ./phasequantizer*2*pi;
    FXkP_quant= FXkP_abs.*exp(i.*FXkP_angle);
    XkP  = ifft2(FXkP_quant); %new output
    %% calculate delta_g
    termA =  abs(image).*(XkP./abs(XkP)); termA(isnan(termA)) = 0;
    termB =  abs(image).*(XkP./abs(XkP)); termB(isnan(termA)) = 0;
    delta_g = termA -XkP+termA-termB;
    %% apply IO step
    new     = zeros(N11,N22);
    %% apply IO on zero paddings
    new(:,1:N2start-1)              = XkP(:,1:N2start-1) + ...
    alpha * delta_g(:,1:N2start-1);
    new(:,N2end+1:N22)              = XkP(:,N2end+1:N22) + ...
    alpha * delta_g(:,N2end+1:N22);
    new(1:N1start-1,N2start:N2end)  = XkP(1:N1start-1,N2start:N2end)...
    + alpha * delta_g(1:N1start-1,N2start:N2end);
    new(N1end+1:N11,N2start:N2end)  = XkP(N1end+1:N11,N2start:N2end)...
    + alpha * delta_g(N1end+1:N11,N2start:N2end);
    %% apply IO in support
    if imagetype == -2
        new(N1start:N1end,N2start:N2end) = ...
            XkP(N1start:N1end,N2start:N2end)+ ...
            alpha * delta_g(N1start:N1end,N2start:N2end);
        Xk             = new;
    end
    %% record residual
    Xtem = zeros(N11,N22);
    Xtem(N1start:N1end,N2start:N2end) = ...
    projection(XkP(N1start:N1end,N2start:N2end),imagetype);
    intensityerr = norm(abs(fft2(Xtem))-Y_abs,'fro')^2;
    intensityerr = sqrt(intensityerr)/norm(Y_abs,'fro');
    iterres(iter)  = intensityerr;
    %% record the change of XkP
    iterchange(iter) = norm(TemMatrix-XkP,'fro')/norm(TemMatrix,'fro');
    %% display data
    fprintf('%6.0f residual = %6.4f p change = %6.6f p \n',iter,...
    100*intensityerr,100*iterchange(iter))
    Xrec   = XkP;
    figure(1)
    subplot(1,2,1)
    imagesc(abs(Xrec(start1:end1,start2:end2)));
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram reconstructed image')
    subplot(1,2,2)
    imagesc(abs(FXkP_quant))
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram')
    drawnow;
end

%% plot the magnitudes of the recovered image

figure(2)
subplot(1,2,1)
imagesc(abs(Xrec));
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title('Hologram reconstructed image')
subplot(1,2,2)
imagesc(abs(FXkP_quant))
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title(['Hologram with ' num2str(magnitudequantizer)...
    ' discrete magnitude and ' num2str(phasequantizer) ' discrete phase'])
savefig(['h_full_' num2str(itermax) 'OO_'... 
    num2str(magnitudequantizer) '_' num2str(phasequantizer) '.fig'])

iterres_256_10_OO = iterres;

save(['h_' num2str(itermax) 'OO_' num2str(magnitudequantizer)...
    '_' num2str(phasequantizer) '_residue.mat'],'iterres_256_10_OO')
%% MIO 100/4
%% Preparations
clear;
close all;
%% parameters
itermax = 50;
alpha = 1;
beta = 1.2;
phasequantizer = 4;
magnitudequantizer = 100;

lena512 = imread('lena512.bmp');
X = im2double(lena512);
i = sqrt(-1);
Sample = 1.2;
loose = 0;
[imn1, imn2] = size(X);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image1   = zeros(imn11,imn22);   % image with zero paddings
start1   = round(1 + (imn11-imn1)/2);   % position of the true image
end1     = start1 + imn1 - 1;
start2   = round(1 + (imn22-imn2)/2);
end2     = start2 + imn2 - 1;
image1(start1:end2,start2:end2) = X;

%% build quandrant put into first quardrant
Sample = 2;
loose = 0;
[imn1, imn2] = size(image1);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image   = zeros(imn11,imn22);                 % image with zero paddings
image(start1:end2,start2:end2) = X;

imagetype = -2;
image     = image.*exp(i*unifrnd(0,2*pi,imn11,imn22));
% complex images have random phases in [0,2*pi]


S1start   = start1 - loose;% In the case where the image has a loose support,
S1end     = end2   + loose;% only a loose position of the true image is given
S2start   = start2 - loose;
S2end     = end2   + loose;

Index  = [imn1 , imn2 ; imn11 , imn22 ; S1start , S1end; S2start , S2end];
% data to input

X0  = generatestart(imagetype,Index);

%% MIO  
N11          = Index(2,1); %size after zeropadding
N22          = Index(2,2); %size after zeropadding
N1start   = start1; %index where the true image start
N1end    = end2;
N2start  = start2;
N2end    = end2;

%% initiation
Xk           = X0;      % Xk
XkP          = X0;      % Pf{Xk}
iterres      = [];      % record residual at each iteration
iterchange   = [];      % record ||XkP{k}-XkP{k-1}||/||XkP{k-1}||

Y_abs = abs(fft2(image));


%% MIO iteration
for iter = 1 : itermax  %iterate over MIO
    TemMatrix = XkP;     % record XkP input from the previous iteration
    %% compute Pf{Xk}
    XkP   = Xk; %new input modified
    FXk  = fft2(XkP);
    FXkP = Y_abs.*exp(i*angle(FXk)); %*(FXk./abs(FXk));
    FXkP_abs= floor(magnitudequantizer*abs(FXkP)./(max(max(abs(FXkP)))))...
      ./magnitudequantizer*max(max(abs(FXkP)));
    FXkP_angle = floor(phasequantizer*angle(FXkP)./(2*pi))...
      ./phasequantizer*2*pi;
    FXkP_quant= FXkP_abs.*exp(i.*FXkP_angle);
    XkP  = ifft2(FXkP_quant); %new output
    %% calculate delta_g
    termA =  abs(image).*(XkP./abs(XkP)); termA(isnan(termA)) = 0;
    termB =  abs(image).*(XkP./abs(XkP)); termB(isnan(termA)) = 0;
    delta_g = termA -XkP+termA-termB;
    %% apply IO step
    new     = zeros(N11,N22);
    %% apply IO on zero paddings
    new(:,1:N2start-1)              = Xk(:,1:N2start-1) +...
        alpha/beta^(iter-1) * delta_g(:,1:N2start-1);
    new(:,N2end+1:N22)              = Xk(:,N2end+1:N22) +...
        alpha/beta^(iter-1) * delta_g(:,N2end+1:N22);
    new(1:N1start-1,N2start:N2end)  = Xk(1:N1start-1,N2start:N2end) +...
        alpha/beta^(iter-1) * delta_g(1:N1start-1,N2start:N2end);
    new(N1end+1:N11,N2start:N2end)  = Xk(N1end+1:N11,N2start:N2end) +...
        alpha/beta^(iter-1) * delta_g(N1end+1:N11,N2start:N2end);
    %% apply IO in support
    if imagetype == -2
        new(N1start:N1end,N2start:N2end) = ...
            Xk(N1start:N1end,N2start:N2end) + alpha/beta^(iter-1) * ...
            delta_g(N1start:N1end,N2start:N2end);
        Xk             = new;
    end
    %% record residual
    Xtem = zeros(N11,N22);
    Xtem(N1start:N1end,N2start:N2end) = ...
    projection(XkP(N1start:N1end,N2start:N2end),imagetype);
    intensityerr = norm(abs(fft2(Xtem))-Y_abs,'fro')^2;
    intensityerr = sqrt(intensityerr)/norm(Y_abs,'fro');
    iterres(iter)  = intensityerr;
    %% record the change of XkP
    iterchange(iter) = norm(TemMatrix-XkP,'fro')/norm(TemMatrix,'fro');
    %% display data
    fprintf('%6.0f residual = %6.4f p change = %6.6f p \n',iter,...
    100*intensityerr,100*iterchange(iter))
    Xrec   = XkP;
    figure(1)
    subplot(1,2,1)
    imagesc(abs(Xrec(start1:end1,start2:end2)));
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram reconstructed image')
    subplot(1,2,2)
    imagesc(abs(FXkP_quant))
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram')
    drawnow;
end


%% plot the magnitudes of the recovered image
figure(2)
subplot(1,2,1)
imagesc(abs(Xrec));
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title('Hologram reconstructed image')
subplot(1,2,2)
imagesc(abs(FXkP_quant))
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title(['Hologram with ' num2str(magnitudequantizer)...
    ' discrete magnitude and ' num2str(phasequantizer) ' discrete phase'])
savefig(['h_full_' num2str(itermax) 'MIO_'... 
    num2str(magnitudequantizer) '_' num2str(phasequantizer) '.fig'])

iterres_100_4_MIO = iterres;
save(['h_' num2str(itermax) 'MIO_' num2str(magnitudequantizer)...
    '_' num2str(phasequantizer) '_residue.mat'],'iterres_100_4_MIO')

%% MIO 256/10
%% Preparations
clear;
close all;
%% parameters
itermax = 50;
alpha = 1;
beta = 1.2;
phasequantizer = 10;
magnitudequantizer = 256;

lena512 = imread('lena512.bmp');
X = im2double(lena512);
i = sqrt(-1);
Sample = 1.2;
loose = 0;
[imn1, imn2] = size(X);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image1   = zeros(imn11,imn22);   % image with zero paddings
start1   = round(1 + (imn11-imn1)/2);   % position of the true image
end1     = start1 + imn1 - 1;
start2   = round(1 + (imn22-imn2)/2);
end2     = start2 + imn2 - 1;
image1(start1:end2,start2:end2) = X;

%% build quandrant put into first quardrant
Sample = 2;
loose = 0;
[imn1, imn2] = size(image1);
imn11       = round(Sample * (imn1+2*loose));
% size of the image with zero paddings
imn22       = round(Sample * (imn2+2*loose));
image   = zeros(imn11,imn22);                 % image with zero paddings
image(start1:end2,start2:end2) = X;

imagetype = -2;
image     = image.*exp(i*unifrnd(0,2*pi,imn11,imn22));
% complex images have random phases in [0,2*pi]


S1start   = start1 - loose;% In the case where the image has a loose support,
S1end     = end2   + loose;% only a loose position of the true image is given
S2start   = start2 - loose;
S2end     = end2   + loose;

Index  = [imn1 , imn2 ; imn11 , imn22 ; S1start , S1end; S2start , S2end];
% data to input

X0  = generatestart(imagetype,Index);

%% MIO  
N11          = Index(2,1); %size after zeropadding
N22          = Index(2,2); %size after zeropadding
N1start   = start1; %index where the true image start
N1end    = end2;
N2start  = start2;
N2end    = end2;

%% initiation
Xk           = X0;      % Xk
XkP          = X0;      % Pf{Xk}
iterres      = [];      % record residual at each iteration
iterchange   = [];      % record ||XkP{k}-XkP{k-1}||/||XkP{k-1}||

Y_abs = abs(fft2(image));

%% MIO iteration
for iter = 1 : itermax  %iterate over MIO
    TemMatrix = XkP;     % record XkP input from the previous iteration
    %% compute Pf{Xk}
    XkP   = Xk; %new input modified
    FXk  = fft2(XkP);
    FXkP = Y_abs.*exp(i*angle(FXk)); %*(FXk./abs(FXk));
    FXkP_abs= floor(magnitudequantizer*abs(FXkP)./(max(max(abs(FXkP)))))...
      ./magnitudequantizer*max(max(abs(FXkP)));
    FXkP_angle = floor(phasequantizer*angle(FXkP)./(2*pi))...
      ./phasequantizer*2*pi;
    FXkP_quant= FXkP_abs.*exp(i.*FXkP_angle);
    XkP  = ifft2(FXkP_quant); %new output
    %% calculate delta_g
    termA =  abs(image).*(XkP./abs(XkP)); termA(isnan(termA)) = 0;
    termB =  abs(image).*(XkP./abs(XkP)); termB(isnan(termA)) = 0;
    delta_g = termA -XkP+termA-termB;
    %% apply IO step
    new     = zeros(N11,N22);
    %% apply IO on zero paddings
    new(:,1:N2start-1)              = Xk(:,1:N2start-1) +...
        alpha/beta^(iter-1) * delta_g(:,1:N2start-1);
    new(:,N2end+1:N22)              = Xk(:,N2end+1:N22) +...
        alpha/beta^(iter-1) * delta_g(:,N2end+1:N22);
    new(1:N1start-1,N2start:N2end)  = Xk(1:N1start-1,N2start:N2end) +...
        alpha/beta^(iter-1) * delta_g(1:N1start-1,N2start:N2end);
    new(N1end+1:N11,N2start:N2end)  = Xk(N1end+1:N11,N2start:N2end) +...
        alpha/beta^(iter-1) * delta_g(N1end+1:N11,N2start:N2end);
    %% apply IO in support
    if imagetype == -2
        new(N1start:N1end,N2start:N2end) = ...
            Xk(N1start:N1end,N2start:N2end) + alpha/beta^(iter-1) * ...
            delta_g(N1start:N1end,N2start:N2end);
        Xk             = new;
    end
    %% record residual
    Xtem = zeros(N11,N22);
    Xtem(N1start:N1end,N2start:N2end) = ...
    projection(XkP(N1start:N1end,N2start:N2end),imagetype);
    intensityerr = norm(abs(fft2(Xtem))-Y_abs,'fro')^2;
    intensityerr = sqrt(intensityerr)/norm(Y_abs,'fro');
    iterres(iter)  = intensityerr;
    %% record the change of XkP
    iterchange(iter) = norm(TemMatrix-XkP,'fro')/norm(TemMatrix,'fro');
    %% display data
    fprintf('%6.0f residual = %6.4f p change = %6.6f p \n',iter,...
    100*intensityerr,100*iterchange(iter))
    Xrec   = XkP;
    figure(1)
    subplot(1,2,1)
    imagesc(abs(Xrec(start1:end1,start2:end2)));
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram reconstructed image')
    subplot(1,2,2)
    imagesc(abs(FXkP_quant))
    axis equal
    axis([1 end1-start1 1 end2-start2])
    colormap(gray)
    title('Hologram')
    drawnow;
end


%% plot the magnitudes of the recovered image

figure(2)
subplot(1,2,1)
imagesc(abs(Xrec));
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title('Hologram reconstructed image')
subplot(1,2,2)
imagesc(abs(FXkP_quant))
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title(['Hologram with ' num2str(magnitudequantizer)...
    ' discrete magnitude and ' num2str(phasequantizer) ' discrete phase'])
savefig(['h_full_' num2str(itermax) 'MIO_'... 
    num2str(magnitudequantizer) '_' num2str(phasequantizer) '.fig'])

iterres_256_10_MIO = iterres;
save(['h_' num2str(itermax) 'MIO_' num2str(magnitudequantizer)...
    '_' num2str(phasequantizer) '_residue.mat'],'iterres_256_10_MIO')

figure(3)
subplot(1,2,1)
imagesc(abs(image));
axis equal
axis([1 imn11 1 imn22])
colormap(gray)
title('original full image')
subplot(1,2,2)
imagesc(abs(image(start1:end1,start2:end2)));
axis equal
axis([1 end1-start1 1 end2-start2])
colormap(gray)
title('original partial image')
savefig('h_original.fig')


load('h_50IO_100_4_residue.mat')
load('h_50IO_256_10_residue.mat')
load('h_50OO_100_4_residue.mat')
load('h_50OO_256_10_residue.mat')
load('h_50MIO_100_4_residue.mat')
load('h_50MIO_256_10_residue.mat')

figure(4)
plot(1:length(iterres_100_4_IO),iterres_100_4_IO,'.-b',...
    1:length(iterres_256_10_IO),iterres_256_10_IO,'.-y',...
    1:length(iterres_100_4_OO),iterres_100_4_OO,'-og',...
    1:length(iterres_256_10_OO),iterres_256_10_OO,'-or',...
    1:length(iterres_100_4_MIO),iterres_100_4_MIO,'.-m',...
    1:length(iterres_256_10_MIO),iterres_256_10_MIO,'.-k')
axis([1 length(iterres_100_4_IO) 0 0.5])
legend('IO with 100/4 Quantized','IO with 256/10 Quantized',...
    'OO with 100/4 Quantized','OO with 256/10 Quantized',...
    'MIO with 100/4 Quantized','MIO with 256/10 Quantized');
title('relative residual at each iteration for Lena holography')
xlabel('iteration')
ylabel('relative residual')
savefig('h_residualcompare.fig')
