clc
clear
load('data.mat')
load('reflect.mat')

%% plot x0,y0,z0 color matching functions 
if 1
    figure(1)
    plot(400:10:700,x,'r',400:10:700,y,'b',400:10:700,z,'g')
    xlabel('wavelength (nm)')
    ylabel('color matching functions')
    hleg3 = legend('x_0','y_0','z_0');
    hold off
end

%% plot l0(?), m0(?), and s0(?) color matching functions corresponding to the long medium and short cones.
if 1
    A_inverse=[0.2430,0.8560,-0.0440;-0.3910,1.1650,0.0870;0.0100,-0.0080,0.5630];
    xyz(1,:)=x;
    xyz(2,:)=y;
    xyz(3,:)=z;
    Cone=A_inverse*xyz;
    
    figure(2)
    plot(400:10:700,Cone(1,:),'r',400:10:700,Cone(2,:),'b',400:10:700,Cone(3,:),'g')
    xlabel('wavelength (nm)')
    ylabel('cone responses')
    hleg3 = legend('l_0','m_0','s_0');
    hold off
end
%%  spectrum of the D65 and fluorescent illuminants
if 1
    figure(3)
    plot(400:10:700,illum1,'r',400:10:700,illum2,'b')
    xlabel('wavelength (nm)')
    ylabel('D_65 and fluorescent illuminations')
    hleg3 = legend('D_65','Fluorescent');
    hold off
end

%% chromaticity plot
if 1
    chromat=zeros(2,31);
    for i=1:31
        chromat(1,i)=Cone(1,i)/sum(Cone(:,i));
        chromat(2,i)=Cone(2,i)/sum(Cone(:,i));
    end
    CIE_1931=[0.73467,0.26533,0.0;0.27376,0.71741,0.00883;0.16658,0.00886,0.82456];
    CIE_709=[0.640,0.330,0.030;0.300,0.600,0.100;0.150,0.060,0.790];
    figure(4)
    plot(chromat(1,:),chromat(2,:),'r')
    hold on
    xlabel('chromaticity x')
    ylabel('chromaticity y')
    hleg3 = legend('chromaticities');
    h=patch(CIE_1931(:,1),CIE_1931(:,2),'w');
    text(CIE_1931(1,1),CIE_1931(1,2),'R_{1931}')
    text(CIE_1931(2,1),CIE_1931(2,2),'G_{1931}')
    text(CIE_1931(3,1),CIE_1931(3,2),'B_{1931}')
    set(h,'facealpha',.1);
    set(h,'edgealpha',.5);
    
    %fill(CIE_709(:,1),CIE_709(:,2),'w')
    h=patch(CIE_709(:,1),CIE_709(:,2),'w');
    text(CIE_709(1,1),CIE_709(1,2),'R_{709}')
    text(CIE_709(2,1),CIE_709(2,2),'G_{709}')
    text(CIE_709(3,1),CIE_709(3,2),'B_{709}')
    set(h,'facealpha',.1);
    set(h,'edgealpha',.5);
    
    D65=[0.3127,0.3290,0.3583];
    EE=[0.3333,0.3333,0.3333];
    plot(0.3127,0.3290,'x',0.3333,0.3333,'*')
    text(0.3127,0.3290,'D65')
    text(0.333,0.3333,'EE')
    hold off
end

%%  two images (same image) obtained from D65 and fluorescent light source
if 1
   CIE_709=[0.640,0.330,0.030;0.300,0.600,0.100;0.150,0.060,0.790];
   I=zeros(170,256,31);
   for i=1:31
       I(:,:,i)=R(:,:,i).*illum1(:,i);
   end
   XYZ=zeros(170,256,13,3);

   for i=1:31
       XYZ(:,:,i,1)=I(:,:,i)*x(1,i);
       XYZ(:,:,i,2)=I(:,:,i)*y(1,i);
       XYZ(:,:,i,3)=I(:,:,i)*z(1,i);
   end
   XYZ=squeeze(sum(XYZ,3));
   
   for i=1:31
       XYZ(:,:,i,1)=I(:,:,i)*x(1,i);
       XYZ(:,:,i,2)=I(:,:,i)*y(1,i);
       XYZ(:,:,i,3)=I(:,:,i)*z(1,i);
   end
   XYZ=squeeze(sum(XYZ,3));
   
   D65=[0.3127;0.3290;0.3583];
   k=inv(CIE_709)*D65;
   k_dia=[k(1),0,0;0,k(2),0;0,0,k(3)];
   M=CIE_709.'*k_dia
   M=inv(M);
   RGB=zeros(size(XYZ));
   RGB(:,:,1)=M(1,1)*XYZ(:,:,1)+M(1,2)*XYZ(:,:,2)+M(1,3)*XYZ(:,:,3);
   RGB(:,:,2)=M(2,1)*XYZ(:,:,1)+M(2,2)*XYZ(:,:,2)+M(2,3)*XYZ(:,:,3);
   RGB(:,:,3)=M(3,1)*XYZ(:,:,1)+M(3,2)*XYZ(:,:,2)+M(3,3)*XYZ(:,:,3);
   RGB=(RGB-min(min(min(RGB))))./(max(max(max(RGB)))-min(min(min(RGB))));
   RGB=RGB.^(1/2.2);
   image(RGB);
   imwrite(uint8(RGB*255),'illum1.tif')
end

if 1
   CIE_709=[0.640,0.330,0.030;0.300,0.600,0.100;0.150,0.060,0.790];
   I=zeros(170,256,31);
   for i=1:31
       I(:,:,i)=R(:,:,i).*illum2(:,i);
   end
   XYZ=zeros(170,256,13,3);

   for i=1:31
       XYZ(:,:,i,1)=I(:,:,i)*x(1,i);
       XYZ(:,:,i,2)=I(:,:,i)*y(1,i);
       XYZ(:,:,i,3)=I(:,:,i)*z(1,i);
   end
   XYZ=squeeze(sum(XYZ,3));
   
   for i=1:31
       XYZ(:,:,i,1)=I(:,:,i)*x(1,i);
       XYZ(:,:,i,2)=I(:,:,i)*y(1,i);
       XYZ(:,:,i,3)=I(:,:,i)*z(1,i);
   end
   XYZ=squeeze(sum(XYZ,3));
   
   D65=[0.3127;0.3290;0.3583];
   k=inv(CIE_709)*D65;
   k_dia=[k(1),0,0;0,k(2),0;0,0,k(3)];
   M=CIE_709.'*k_dia
   M=inv(M);
   RGB=zeros(size(XYZ));
   RGB(:,:,1)=M(1,1)*XYZ(:,:,1)+M(1,2)*XYZ(:,:,2)+M(1,3)*XYZ(:,:,3);
   RGB(:,:,2)=M(2,1)*XYZ(:,:,1)+M(2,2)*XYZ(:,:,2)+M(2,3)*XYZ(:,:,3);
   RGB(:,:,3)=M(3,1)*XYZ(:,:,1)+M(3,2)*XYZ(:,:,2)+M(3,3)*XYZ(:,:,3);
   RGB=(RGB-min(min(min(RGB))))./(max(max(max(RGB)))-min(min(min(RGB))));
   RGB=RGB.^(1/2.2);
   image(RGB);
   imwrite(uint8(RGB*255),'illum2.tif')
end

%% Compare Monitor chromaticity with 709 RGB
if 1
    [X,Y]=meshgrid(0:0.005:1,1:-0.005:0);
    Z=1-X-Y;
    M=CIE_709.';
    M=inv(M);
    RGB=zeros([size(X),3]);
    RGB(:,:,1)=M(1,1)*X+M(1,2)*Y+M(1,3)*Z;
    RGB(:,:,2)=M(2,1)*X+M(2,2)*Y+M(2,3)*Z;
    RGB(:,:,3)=M(3,1)*X+M(3,2)*Y+M(3,3)*Z;
    RGB(logical(RGB<0))=1;
    figure(5)
    image([0:.005:1],[0:.005:1],RGB);
    axis('xy');
    hold on
    plot(chromat(1,:),chromat(2,:),'r')
    hold off
    
    
end

