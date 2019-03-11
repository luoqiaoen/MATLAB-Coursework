%The Cornu Spiral
clear
dx=0.01;
t=0:dx: 10;
x=0;
y=zeros(401);
y=y(1,:);
x1 =0;
y1 =zeros(401 );
y1=y1(1,:);
%dx implies the increment distance between the intervals
%"t" represents alpha
%"x" represents C as a function of alpha
%"y" represents S as a function of alpha
for m = 1:401 %these "for" loops are used to evaluate the integrals
    for n = 1 :m
        x(n)=cos((pi.*t(n).^2)./2).*dx;
        y(m) = y(m)+x(n);
        x1(n)=sin((pi.*t(n).^2)./2).*dx;
        y1(m)=y1(m)+x1(n);
    end
end
Y=fliplr(y);
Y1=fliplr(y1);
%this flips the current original graph to the position
%in the third quadrant
T=-4:0.01:4;
BY(1:401)=-Y;
BY(401:801)=y;
BY1(1:401)=-Y1;
BY1(401:801) =y1 ;
plot3(BY,BY1,T)
view(0,90)
grid on
%this combines the two existing graphs into one spiral
%plots the original spiral (FIGURE l)
%rotates it for the viewer to view from a birds eye view
figure(2)%plots the second figure which is a side view of C and S in relation to alpha
plot(T,BY,'r') %this plots the C graph in relation to alpha
hold on
plot(T,BY1) %this plots the S graph in relation to alpha
grid on