function y=bit2num(th_bit,b,thetmin,thetmax)
y=0;
for i=1:b
   y=y+th_bit(b+1-i)*2^(i-1);
end
y=thetmin+(thetmax-thetmin)*y/(2^b-1);
