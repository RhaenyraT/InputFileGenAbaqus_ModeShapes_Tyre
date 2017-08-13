
%   Tyre parameters
E=10.4e+9;            %     Young's modulus
b=15.2e-2;            %     Belt ring section width in m
d=3.1e-3;             %     Belt ring section heigth in m
I=(b*d^3)/12;
D=(E*d^3)/12;         %     Bending stiffness of belt ring 
K=E*d;                %     Membrane stiffness of belt ring
Dnorm=D*b;              %     Normalized bending stiffness of belt ring
Knorm=K*b;              %     Normalized Membrane stiffness of belt ring
kr=(192.9e+6)*b;          %     Radial stiffness of sidewall in kg/m.sec^2
kt=(648.7e+5)*b;          %     Tangential stiffness of sidewall in kg/m.sec^2
rho=8.1e+3;           %     Density in in kg/m^3
A=b*d;                %     Belt ring section Area in m^2
r=0.3;                %     Belt ring radius in m
P0=2.2e+5;
%zEGELAAR -undeformed tyre sidewall at 2.2e+5 inflation pressure
% r= 0.3                Belt ring radius in m
% rs=.054               radius of tyre sidewall arc
% t=.1                  thickness of sidewall 
% ls=.121               length of sidewall arc
% G=1.6e+6              shear modulus of sidewall
% w/2=62.3deg           half angle of sidewall
% L=.0867               height of tyre sidewall
%------------------------------------------------------------------------------------------------------------------------------------------------------------
%%
%   Preallocating
Nbelt=800;                       %     number of discretized belt Nodes
Nmodes=40;                       %     number of Modes
Del=zeros(Nmodes,1);         
B=zeros(Nmodes,1);
Bs=zeros(Nmodes,1);
Dels=zeros(Nmodes,1);
ValNatfreq=zeros(Nmodes,2);
Tm_div_Rm=zeros(Nmodes,2);
phid=zeros(Nmodes,1);
U=zeros(2*Nbelt,1);
theta = (2*pi/Nbelt)*(0:(Nbelt-1));    %     central angle of belt arch,

ModeStructArr(1:Nmodes) = struct('frq',0,'Tm_div_Rm',1,'theta',theta,'dr',zeros(size(theta)),'dt',zeros(size(theta)));

RadialModes = struct('Nmode',1:Nmodes,'ModeData',ModeStructArr);
TangentialModes = struct('Nmode',1:Nmodes,'ModeData',ModeStructArr);
% ------------------------------------------------------------------------------------------------------------------------------------------------------------
%%
%   NATURAL FREQUENCIES OF BELT RING
for n=1:Nmodes
        %    Precalculated terms
        Lambda4 = (A *rho)^2; 
        Lambda2 = -A*kr*rho -A*kt*rho -(A*E*I*(n^2)*rho/(r^4)) - (A*E*I*(n^4)*rho/(r^4)) - (A*E*A*rho/(r^2)) - (A*E*A*rho*(n^2)/(r^2)) + ((A*rho*P0*b)/r) - ((A*rho*P0*b*n^2)/r);
        Lambda0 = (kr*kt) + ((E*I*E*A*(n^2))/(r^6)) - ((2*E*I*E*A*(n^4))/(r^6)) + ((E*I*E*A*(n^6))/(r^6)) + ((E*I*kr*(n^2))/(r^4)) + ((E*I*kt*(n^4))/(r^4)) + ((E*A*kt)/(r^2)) + ((E*A*kr*(n^2))/(r^2)) + ((P0*b)*((-(E*I*n^2)/r^5)+((E*I*n^4)/r^5)-((E*A*n^2)/r^3)+((E*A*n^4)/r^3)-(kt/r)+((kr*n^2)/r)));
        
        Dels = Lambda2^2 - (4*Lambda4*Lambda0);
        
        
        %   Calculating natural frequencies of each mode in rad/s
        RadialModes.ModeData(n).frq = (sqrt((-Lambda2 - sqrt(Dels))/(2*Lambda4))) * .159; % radial
        TangentialModes.ModeData(n).frq = (sqrt((-Lambda2 + sqrt(Dels))/(2*Lambda4))) * .159; % tangential
        
        
        %   Ratio of radial to tangential deformation amplitude
        RadialModes.ModeData(n).Tm_div_Rm = -(E*I*(n^4)/(r^4) + E*A/(r^2) + kr - rho*A*(RadialModes.ModeData(n).frq^2))/(E*I*(n^3)/(r^4) + E*A*n/(r^2));      %Radial
        TangentialModes.ModeData(n).Tm_div_Rm = -(E*I*(n^4)/(r^4) + E*A/(r^2) + kr - rho*A*(TangentialModes.ModeData(n).frq^2))/(E*I*(n^3)/(r^4) + E*A*n/(r^2));      %Tangential
        
        %   Find out phi angle between the mode shapes
%         phid(n)=180/(2*n);

        % Calculate radial/tangential displacements
        RadialModes.ModeData(n).dr = cos(n*RadialModes.ModeData(n).theta);
        RadialModes.ModeData(n).dt = RadialModes.ModeData(n).Tm_div_Rm*sin(n*RadialModes.ModeData(n).theta);
        
        TangentialModes.ModeData(n).dr = cos(n*TangentialModes.ModeData(n).theta);
        TangentialModes.ModeData(n).dt = TangentialModes.ModeData(n).Tm_div_Rm*sin(n*TangentialModes.ModeData(n).theta);

%         fprintf('For mode %f the Radial natural frequency is %f Hz and the Tangential natural frequency is %f Hz\n',n,RadialModes.ModeData(n).frq,TangentialModes.ModeData(n).frq);
%         fprintf('The ratio of radial to tangential deformation amplitude for the above radial mode is %f and for Tangential mode is %f\n\n',RadialModes.ModeData(n).Tm_div_Rm, TangentialModes.ModeData(n).Tm_div_Rm);
end
%%
%Plot radial and tangential natural frequencies
figure,plot(1:Nmodes,[RadialModes.ModeData(1:Nmodes).frq],'-.r*')
xlabel('Mode number')
ylabel('Natural frequency in Hz')
% ylabel('Radial natural frequency in Hz')
hold ,plot(1:Nmodes,[TangentialModes.ModeData(1:Nmodes).frq],':bs');
xlabel('Mode number')
% ylabel('Tangential natural frequency in Hz')

%ZERO MODES 
omega1=sqrt(kt/(rho*A))*.159;
omega2=sqrt((E*A+(kr*r^2))/(r^2*rho*A))*.159;

%plotting mode shapes in 3d polar form
% for n=1:Nmodes
%     axis equal;
%     ux = (30+RadialModes.ModeData(n).dr) .* cos(RadialModes.ModeData(n).theta) - RadialModes.ModeData(n).dt .* sin(RadialModes.ModeData(n).theta);
%     uy = (30+RadialModes.ModeData(n).dr) .* sin(RadialModes.ModeData(n).theta) + RadialModes.ModeData(n).dt .* cos(RadialModes.ModeData(n).theta);
%     subplot(8,5,n)
%     plot(ux,uy,'linewidth',4);
%     hold on
%     plot(30*cos(theta),30*sin(theta),'g','linewidth',2);
%     axis square
% end
    
%modal participation factors


%Modal force 
