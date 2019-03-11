function [detS, ri, M]=Ray_s(ro, do, f,z)
%This function is for output ray vector of
%a single lens system
To=[1, do;0,1];
Sf=[1,0;-(1/f) ,1];
Ti=[1,z;0,1];
S=Ti*Sf*To;
%Checking determinant for overall matrix
detS=det(S);
%"image" ray coordinate is ri
ri=S*ro;
%check magnification
M = S(1,1);
end
