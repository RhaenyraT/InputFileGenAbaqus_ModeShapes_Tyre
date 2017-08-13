%Tyre parameters

E=10.4e+9;            %Young's modulus
b=15.2e-2;            %Belt ring section width in m
h=3.1e-3;             %Belt ring section heigth in m
I=b*h^3/12;           %second moment of inertia
D=(E*h^3)/12;         %Bending stiffness of belt ring 
K=E*h;                %Membrane stiffness of belt ring
Dnorm=D*b;            %Normalized bending stiffness of belt ring
Knorm=K*b;            %Normalized Membrane stiffness of belt ring
kr=192.9e+6;          %Radial stiffness of sidewall in kg/m.sec^2
kt=648.7e+5;          %Tangential stiffness of sidewall in kg/m.sec^2
krnorm=kr*b;          %Normalized Radial stiffness of sidewall in kg/m.sec^2
ktnorm=kt*b;          %Normalized Tangenial stiffness of sidewall in kg/m.sec^2
rho=8.1e+9;           %Density in in kg/m^3
A=b*h;                %Belt ring section Area in m^2
R=0.3;                %Belt ring radius in m
rhof=rho*h;
hf=h;                 %homogenous thickness of the shell
%From kinetic energy considerations, assuming no surges in the homogeneous
%foundation itself, one-third of its mass per unit area has to be added to the h of the shell.
m=rho*A+((1/3)*rhof*hf*b); 
%natfreq =  natural frequencies of the modes;
%n =  number of modes
n=40;
%Natfreq=zeros(length(n),4);
%for i=1:length(n)

%CharacEq=[((D*n.^4/R^4)+kr-(K/R^4)-(rho*A*omega^2)) ((D*n.^3/R^4)+(K*n/R^2));((D*n.^3/R^4)+(K*n/R^2)) ((D*n.^2/R^4)+kt+(K*n.^2/R^2)-(rho*A*omega^2))];
%CharacEqDet=det(CharacEq);
%charac=((kr*kt)-(((D*K*n.^2)/R^6)*(1+2*n.^2-n.^4))+(((D*n.^2)/R^4)*(kr+(kt*n.^2)))+((K/R^2)*((kr*n.^2)-kt))-((2*K*n.^2)/R^2))+omeg^2*((-rho*A*(kr+kt))-(((rho*A*D*n.^2)/R^4)*(1+n.^2))+(((rho*A*K)/R^2)*(1-n.^2)))+omeg^4*(rho^2*A^2);



charac1=[rho^2*A^2 0 ((-rho*A*(kr+kt))-(((rho*A*Dnorm*n.^2)/R^4)*(1+n.^2))-(((rho*A*Knorm)/R^2)*(1+n.^2))) 0 ((kr*kt)+(((Dnorm*Knorm*n.^2)/R^6)*(1-2*n.^2+n.^4))+(((Dnorm*n.^2)/R^4)*(kr+(kt*n.^2)))+((Knorm/R^2)*((kr*n.^2)+kt)))];
Roots(charac1)








