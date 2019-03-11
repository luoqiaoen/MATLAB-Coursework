function y=num2bit(th_element,b,thetmin,thetmax)
d=(thetmax-thetmin)/(2^b-1);
rnd=round((th_element-thetmin)/d);
y=zeros(1,b);
for i=b:-1:1
   ratio=rnd/(2^(i-1));
   if ratio >= 1
      y(b-i+1)=1;
      rnd=rnd-2^(i-1);
   else
   end
end
