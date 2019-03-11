function [phi, Ex, Ey] = FEM_Strip(d_x, d_y,s_length, s_vol)

dimension_x=d_x;
dimension_y=d_y;
strip_voltage=s_vol;
strip_index_x=1:s_length;
BC_voltage=0;
BC_value = 1e60;


strip_index_y= round(dimension_y/3-1/6);
epsilon_index= strip_index_y-1;
epsilon_value=[1 100];

element_no=(dimension_x-1).*(dimension_y-1).*2;
element_x=zeros(element_no,3);
element_y=zeros(element_no,3);
K_domain=zeros(dimension_x.*dimension_y,dimension_x.*dimension_y);
K_temp=zeros(dimension_x.*dimension_y,dimension_x.*dimension_y);
epsilon=ones(1,element_no);

for a=1:epsilon_index-1
    epsilon(1+(a-1).*(dimension_x-1)*2:(a-1).*(dimension_x-1)*2+(dimension_x-1)*2)=epsilon_value(1);
end

for b=epsilon_index:dimension_y-1
    epsilon((b-1).*(dimension_x-1)*2+1:(b-1).*(dimension_x-1)*2+(dimension_x-1)*2)=epsilon_value(2);
end    

alpha_e=epsilon;

a1=zeros(1,element_no);a2=zeros(1,element_no);a3=zeros(1,element_no);
b1=zeros(1,element_no);b2=zeros(1,element_no);b3=zeros(1,element_no);
c1=zeros(1,element_no);c2=zeros(1,element_no);c3=zeros(1,element_no);
area=zeros(1,element_no);

BC=zeros(1,dimension_x.*dimension_y);
    
for n=1:dimension_y-1
    for m=1:dimension_x-1 
        element_x((n-1).*2.*(dimension_x-1)+(m.*2-1),1:3)=[m  m  m+1];
        element_y((n-1).*2.*(dimension_x-1)+(m.*2-1),1:3)=[n  n+1  n+1];
        element_x((n-1).*2.*(dimension_x-1)+(m.*2),1:3)=[m  m+1  m+1];
        element_y((n-1).*2.*(dimension_x-1)+(m.*2),1:3)=[n  n+1  n];
    end    
end    

for o=1:element_no
    a1(o)=element_x(o,2).*element_y(o,3)-element_y(o,2).*element_x(o,3);
    a2(o)=element_x(o,3).*element_y(o,1)-element_y(o,3).*element_x(o,1);
    a3(o)=element_x(o,1).*element_y(o,2)-element_y(o,1).*element_x(o,2);
    b1(o)=element_y(o,2)-element_y(o,3);
    b2(o)=element_y(o,3)-element_y(o,1);
    b3(o)=element_y(o,1)-element_y(o,2);
    c1(o)=element_x(o,3)-element_x(o,2);
    c2(o)=element_x(o,1)-element_x(o,3);
    c3(o)=element_x(o,2)-element_x(o,1);    
end    

for p=1:element_no
    area(p)=-0.5.*(b1(p).*c2(p)-b2(p).*c1(p));
end    

for q=1:dimension_y-1
    for r=1:dimension_x-1
        alpha1=r+(q-1).*dimension_x;
        beta1=r+q.*dimension_x;
        gamma1=r+q.*dimension_x+1;
        alpha2=r+(q-1).*dimension_x;
        beta2=r+q.*dimension_x+1;
        gamma2=r+1+(q-1).*dimension_x;
        element_index1=r.*2-1+(q-1).*(dimension_x-1).*2;
        element_index2=r.*2+(q-1).*(dimension_x-1).*2;
        
        K_domain(alpha1,alpha1)=K_temp(alpha1,alpha1)+alpha_e(element_index1)./(4.*area(element_index1)).*(b1(element_index1).*b1(element_index1)+c1(element_index1).*c1(element_index1));
        K_domain(alpha1,beta1)=K_temp(alpha1,beta1)+alpha_e(element_index1)./(4.*area(element_index1)).*(b1(element_index1).*b2(element_index1)+c1(element_index1).*c2(element_index1));
        K_domain(alpha1,gamma1)=K_temp(alpha1,gamma1)+alpha_e(element_index1)./(4.*area(element_index1)).*(b1(element_index1).*b3(element_index1)+c1(element_index1).*c3(element_index1)); 
        
        K_domain(beta1,alpha1)=K_temp(beta1,alpha1)+alpha_e(element_index1)./(4.*area(element_index1)).*(b2(element_index1).*b1(element_index1)+c2(element_index1).*c1(element_index1));
        K_domain(beta1,beta1)=K_temp(beta1,beta1)+alpha_e(element_index1)./(4.*area(element_index1)).*(b2(element_index1).*b2(element_index1)+c2(element_index1).*c2(element_index1));
        K_domain(beta1,gamma1)=K_temp(beta1,gamma1)+alpha_e(element_index1)./(4.*area(element_index1)).*(b2(element_index1).*b3(element_index1)+c2(element_index1).*c3(element_index1));
        
        K_domain(gamma1,alpha1)=K_temp(gamma1,alpha1)+alpha_e(element_index1)./(4.*area(element_index1)).*(b3(element_index1).*b1(element_index1)+c3(element_index1).*c1(element_index1));
        K_domain(gamma1,beta1)=K_temp(gamma1,beta1)+alpha_e(element_index1)./(4.*area(element_index1)).*(b3(element_index1).*b2(element_index1)+c3(element_index1).*c2(element_index1));
        K_domain(gamma1,gamma1)=K_temp(gamma1,gamma1)+alpha_e(element_index1)./(4.*area(element_index1)).*(b3(element_index1).*b3(element_index1)+c3(element_index1).*c3(element_index1));
        
        K_temp=K_domain;
        
        K_domain(alpha2,alpha2)=K_temp(alpha2,alpha2)+alpha_e(element_index2)./(4.*area(element_index2)).*(b1(element_index2).*b1(element_index2)+c1(element_index2).*c1(element_index2));
        K_domain(alpha2,beta2)=K_temp(alpha2,beta2)+alpha_e(element_index2)./(4.*area(element_index2)).*(b1(element_index2).*b2(element_index2)+c1(element_index2).*c2(element_index2));
        K_domain(alpha2,gamma2)=K_temp(alpha2,gamma2)+alpha_e(element_index2)./(4.*area(element_index2)).*(b1(element_index2).*b3(element_index2)+c1(element_index2).*c3(element_index2)); 
        
        K_domain(beta2,alpha2)=K_temp(beta2,alpha2)+alpha_e(element_index2)./(4.*area(element_index2)).*(b2(element_index2).*b1(element_index2)+c2(element_index2).*c1(element_index2));
        K_domain(beta2,beta2)=K_temp(beta2,beta2)+alpha_e(element_index2)./(4.*area(element_index2)).*(b2(element_index2).*b2(element_index2)+c2(element_index2).*c2(element_index2));
        K_domain(beta2,gamma2)=K_temp(beta2,gamma2)+alpha_e(element_index2)./(4.*area(element_index2)).*(b2(element_index2).*b3(element_index2)+c2(element_index2).*c3(element_index2));
        
        K_domain(gamma2,alpha2)=K_temp(gamma2,alpha2)+alpha_e(element_index2)./(4.*area(element_index2)).*(b3(element_index2).*b1(element_index2)+c3(element_index2).*c1(element_index2));
        K_domain(gamma2,beta2)=K_temp(gamma2,beta2)+alpha_e(element_index2)./(4.*area(element_index2)).*(b3(element_index2).*b2(element_index2)+c3(element_index2).*c2(element_index2));
        K_domain(gamma2,gamma2)=K_temp(gamma2,gamma2)+alpha_e(element_index2)./(4.*area(element_index2)).*(b3(element_index2).*b3(element_index2)+c3(element_index2).*c3(element_index2));
        
        K_temp=K_domain;
    end    
end    

for N=1:dimension_x
    K_domain(N,N)=BC_value;
    K_domain((dimension_y-1).*(dimension_x)+N,(dimension_y-1).*(dimension_x)+N)=BC_value;
    BC(N)=BC_voltage.*BC_value;
    BC((dimension_y-1).*(dimension_x)+N)=BC_voltage.*BC_value;
end    

for M=1:dimension_y
    K_domain(dimension_x.*M,dimension_x.*M)=BC_value;
    BC(M)=BC_voltage.*BC_value;
end        

for i=strip_index_y
    for j=strip_index_x
        BC((i-1).*dimension_x+j)=strip_voltage.*BC_value;
        K_domain((i-1).*dimension_x+j,(i-1).*dimension_x+j)=BC_value;
    end    
end   

phi_temp= K_domain\BC';
phi=transpose(reshape(phi_temp,dimension_x,dimension_y));
[Ex,Ey]=gradient(phi);
Ex=-Ex;Ey=-Ey;
% figure(1)
% hold on
% for n=1:1:dimension_x-1
%     for m=1:1:dimension_y-1
%         plot([n n+1 n+1 n],[m m m+1 m]);
%     end
% end
% hold off
% axis image
% title('Mesh Scheme')

