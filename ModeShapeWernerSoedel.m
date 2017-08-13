%TYRE parameters
E       =  10.4e+9;            %Young's modulus
b       =  15.2e-2;            %Belt ring section width in m
h       =  3.1e-3;             %Belt ring section heigth in m
I       =  b*h^3/12;           %second moment of inertia
D       =  (E*h^3)/12;         %Bending stiffness of belt ring 
K       =  E*h;                %Membrane stiffness of belt ring
Dnorm   =   D*b;                %Normalized bending stiffness of belt ring
Knorm   =   K*b;                %Normalized Membrane stiffness of belt ring
kr      =   192.9e+6;          %Radial stiffness of sidewall in kg/m.sec^2
kt      =   648.7e+5;          %Tangential stiffness of sidewall in kg/m.sec^2
krnorm  =   kr*b;               %Normalized Radial stiffness of sidewall in kg/m.sec^2
ktnorm  =   kt*b;               %Normalized Tangenial stiffness of sidewall in kg/m.sec^2
rho     =   8.1e+3;             %Density in in kg/m^3
A       =   b*h;                %Belt ring section Area in m^2
a       =   0.3;                %Belt ring radius in m
rhof    =   rho*h;
hf      =   h;                 %homogenous thickness of the shell
m       =(rho*A) + ((1/3)*rhof*hf*b); 
n       =1:8;

  
%Precalculated terms
%K1=zeros(length(n),1);
K1=(((n.^2+1)/(a^2*m)) .* (((n.^2*E*I)./a^2)+(E*A))) + ((ktnorm+krnorm)/m);
%K2=zeros(length(n),1);
K2=(((n.^2.*(n.^2-1).^2)/(a^6*m^2))*E^2*I*A) + ((krnorm*ktnorm)+(krnorm*((n.^2.*E)/a^2).*((I/a^2)+A))+(ktnorm*(E/a^2)*(((n.^4*I)/a^2)+A)))/(m^2);

%Calculating natural frequencies of each mode

ValNatfreq1= (sqrt ((K1/2).*(1-(sqrt(1-(4*(K2./K1.^2))))))) *.159 ;
ValNatfreq2= (sqrt ((K1/2).*(1+(sqrt(1-(4*(K2./K1.^2))))))) *.159;

plot(n,ValNatfreq1,'-*r');
title('Mode frequencies - Soedel')
xlabel('Mode number');
ylabel('Natural frequency in Hz')

figure,plot(n,ValNatfreq2,'-og');

%natfreq =  natural frequencies of the modes;
%n =  number of modes
% %From kinetic energy considerations, assuming no surges in the homogeneous
%foundation itself, one-third of its mass per unit area has to be added to the h of the shell.
% for i=1:length(n)
%     fprintf('For mode %f the Radial natural frequency is %f Hz and the Tangential natural frequency is %f Hz\n\n',n(i),ValNatfreq1(i),ValNatfreq2(i));
% end
%Ratio of radial to tangential deformation amplitude
%ValNatfreq1=zeros(n,1);
%ValNatfreq2=zeros(n,1);
%ValNatfreq1=zeros(length(n),1);
%RTmbyRm1=((m*ValNatfreq1(n)^2)-krnorm-((n^4*E*I)/a^4)-((E*A)/a^2))/((n^3*E*I/a^4)+(n^2*E*A/a^2));
%RTmbyRm2=((m*ValNatfreq2(n)^2)-krnorm-((n^4*E*I)/a^4)-((E*A)/a^2))/((n^3*E*I/a^4)+(n^2*E*A/a^2));
%TTmbyRm1=((n^3*E*I/a^4)+(n^2*E*A/a^2))/((m*ValNatfreq1(n)^2)-ktnorm-((n^4*E*I)/a^4)-((n^2*E*A)/a^2));
%TTmbyRm2=((n^3*E*I/a^4)+(n^2*E*A/a^2))/((m*ValNatfreq2(n)^2)-ktnorm-((n^4*E*I)/a^4)-((n^2*E*A)/a^2));
%fprintf('The ratio of radial to tangential deformation amplitude for %f is %f and for Tangential mode is %f\n\n',RTmbyRm,TTmbyRm);
