
%**************************************************************************************************************************************************
%  FULL TIRE INPUTS ^ SI UNITS ^ 
%**************************************************************************************************************************************************

%  TREAD stiffness, dampin
Kr_gen                  =  838.5E+6;      
Kt_gen                  =  120E+6;     
Cr_gen                  =  0;       
Ct_gen                  =  0;        
ltread                  =  15.5e-3;

%  SIDEWALL stiffness
kt_gen                  =  178115.135;      
%  BELT
sinput.Tire_Belt_Radius =  0.3504;
EI                      =  1.91;           
EA                      =  9.9e+6;       
rhoA                    =  4.68;
Poisson                 =  0.3;

%  WHEEL mass, inertia
mw                      =  1.9 ;
Ixx                     =  4.08037e-002;
Iyy                     =  6.86247e-002;
Izz                     =  4.08037e-002;
Ixy                     =  1.28938e-010;
Iyz                     =  1.94720e-011;
Izx                     =  -1.15992e-007;

% Tire deformation
Tire_Def                =  0.02   ;

%**************************************************************************************************************************************************
%  SIDEWALL INPUTS & EVALUATION *_* 
%**************************************************************************************************************************************************

sinput.NoofBeltElem     =  200;     
sinput.P                =  220000;
sinput.urneg            =  -(0:.001:.05)';
sinput.urpos            =  (0:.001:.003)';  
%INITIAL GUESS
sinput.w                =  1.74 ;       
sinput.lstring          =  .0869;  %0.0862 for 210e3  , .0869
sinput.Tire_Rim_Radius  =  0.2756; %.2763  for 210e3, 0.2756

% NOT ABOVE .004 m

[qmem_P0neg]   =   NL_sidewallneg(sinput);
Tr_neg         =  -2*qmem_P0neg*2*pi*sinput.Tire_Belt_Radius*sinput.P/sinput.NoofBeltElem;
SW_neg         =  [sinput.urneg Tr_neg];
SW_negSort     =   sortrows(SW_neg);

[qmem_P0pos]   =  NL_sidewallpos(sinput);
Tr_pos         =  -2*qmem_P0pos*2*pi*sinput.Tire_Belt_Radius*sinput.P/sinput.NoofBeltElem;
SW_pos         =  [sinput.urpos Tr_pos];
SW_posSort     =  sortrows(SW_pos);

%INPUT for .inp  [ur Tr]
SW             = [SW_negSort(1:end-1,:);SW_posSort(2:end,:)]; 
SW(:,1)        = SW(:,1); 

%Interpolated Value of Tr at '0'deformation.
Tr_0           = interp1(SW(:,1),SW(:,2),0,'linear');

%CHECK : +ve F for Tension, -ve for Compression
figure,plot(SW(:,1),SW(:,2),'-*r');   
%**************************************************************************************************************************************************

%  Area, Youngs modulus, density , Inflation 
beam_section_height =  0.00152191;
beam_section_width  =  .152; 
A                   =   beam_section_width*beam_section_height;
E                   =   EA / A;
rho                 =   rhoA / A;
b_prime             =  .135;
P_tire              =  sinput.P * b_prime; %sinput.P * b_prime;

%  NUMBER OF ELEMENTS , NODES
NoofRimElem      =  sinput.NoofBeltElem;
NoofConnElem     =  sinput.NoofBeltElem;
NodeBelt         =  sinput.NoofBeltElem;
NodeRim          =  NodeBelt;
NodeTread        =  NodeBelt;

% TIRE DISCRETIZATION ANGLE
theta  =  ((0:(360/NodeBelt):(360-(360/NodeBelt)))*(pi/180))';

% DISCRETIZED TIRE STIFFNESS 
Kr     =  Kr_gen*(360/NodeTread)*(pi/180)*sinput.Tire_Belt_Radius;
Kt     =  Kt_gen*(360/NodeTread)*(pi/180)*sinput.Tire_Belt_Radius;
Cr     =  Cr_gen*(360/NodeTread)*(pi/180)*sinput.Tire_Belt_Radius;
Ct     =  Ct_gen*(360/NodeTread)*(pi/180)*sinput.Tire_Belt_Radius;
kt     =  kt_gen*(360/NodeTread)*(pi/180)*((sinput.Tire_Belt_Radius));

% Road placed 1 mm below the tyre at the beginning 
road_reference   =  -(sinput.Tire_Belt_Radius+ltread+.001);

% NODAL POSITION of TREAD, RIM, SIDEWALL
Yb     =  (sinput.Tire_Belt_Radius*sin(theta))  ;                        %Yb=vertical global coordinate of a belt point
Xb     =  (sinput.Tire_Belt_Radius*cos(theta)) ;                         %Yb=vertical global coordinate of a belt point
Yt     =  (Yb  +  ((ltread)*sin(theta)))  ;                       %Yt=global vertical coordinate of the tread element's tip
Xt     =  (Xb  +  ((ltread)*cos(theta)))  ;
Xr     =  (sinput.Tire_Rim_Radius*cos(theta));
Yr     =  (sinput.Tire_Rim_Radius*sin(theta)); 

% NSETS NUMBERS
Nodes_Tread            = (1:NodeTread)';
% Nodes_Tread_Tspr       = (NodeTread+1:NodeTread*2)';
Nodes_Belt             = ((NodeTread+1):(NodeTread+NodeBelt))';
% Nodes_Belt_Tspr        = (((NodeTread*2+NodeBelt)+1):((NodeTread*2+NodeBelt)+NodeBelt))';
Nodes_Rim              = (((NodeTread*2)+1):((NodeTread*2+NodeBelt)))';
% Nodes_ghost            = (((NodeTread*2+NodeBelt)+NodeBelt)+NodeBelt+1):(((NodeTread*2+NodeBelt)+NodeBelt)+NodeBelt+NodeBelt);

Belt_Normal_x  = cos(theta);
Belt_Normal_y  = sin(theta);
Belt_Normal_z  = zeros(length(theta),1);

% CONNECTOR LOCAL ORIENTATIONS
X_x_r   =  cos(theta);
X_y_r   =  sin(theta);
X_z_r   =  zeros(NoofConnElem,1);
Y_x_r   =  sin(theta);
Y_y_r   =  -cos(theta);
Y_z_r   =  zeros(NoofConnElem,1);
X_x_t   =  -sin(theta);
X_y_t   =  cos(theta);
X_z_t   =  zeros(NoofConnElem,1);
Y_x_t   =  cos(theta);
Y_y_t   =  sin(theta);
Y_z_t   =  zeros(NoofConnElem,1);

% ELSETS NUMBERS
ELSET_Belt   =  (50001:(50000+sinput.NoofBeltElem))';
ELSET_Rim    =  (60001:(60000+NoofRimElem))';

% CONNECTOR ELSET numbers
a      =  (85001:(85000+NoofConnElem))';
b      =  ((85000+NoofConnElem+1):(85000+(NoofConnElem*2)))';
c      =  (((85000+(NoofConnElem*2))+1):(85000+(NoofConnElem*3)))';

%**************************************************************************************************************************************************
%  WRITING SUPPORTING INPUT FILES 
%**************************************************************************************************************************************************
Nodes_Beltfile = fopen('E:/2DTire LoadVsDef/200elem/Nodes_Belt.inp','w+');
fprintf(Nodes_Beltfile,'*node,nset=belt\n');
for s=1:NodeBelt
fprintf(Nodes_Beltfile,'%d,%f,%f\n',Nodes_Belt(s),Xb(s),Yb(s));
end
fclose(Nodes_Beltfile);

%****************************************************************************************************************************

Nodes_Treadfile = fopen('E:/2DTire LoadVsDef/200elem/Nodes_Tread.inp','w+');
fprintf(Nodes_Treadfile,'*node,nset=tread\n');
for s=1:NodeTread
fprintf(Nodes_Treadfile,'%d,%f,%f\n',Nodes_Tread(s),Xt(s),Yt(s));
end
fclose(Nodes_Treadfile);

%****************************************************************************************************************************

Nodes_Rimfile = fopen('E:/2DTire LoadVsDef/200elem/Nodes_Rim.inp','w+');
fprintf(Nodes_Rimfile,'*node,nset=rim\n');
for s=1:NodeBelt
fprintf(Nodes_Rimfile,'%d,%f,%f\n',Nodes_Rim(s),Xr(s),Yr(s));
end
fclose(Nodes_Rimfile);

%****************************************************************************************************************************

Elements_Beltfile = fopen('E:/2DTire LoadVsDef/200elem/Elements_Belt.inp','w+');
for s=1:NodeBelt-1
fprintf(Elements_Beltfile,'%d,%d,%d\n',ELSET_Belt(s),Nodes_Belt(s),Nodes_Belt(s+1));
end
fprintf(Elements_Beltfile,'%d,%d,%d\n',ELSET_Belt(NodeBelt),Nodes_Belt(NodeBelt),Nodes_Belt(1));
fclose(Elements_Beltfile);

%****************************************************************************************************************************

Elements_Rimfile = fopen('E:/2DTire LoadVsDef/200elem/Elements_Rim.inp','w+');
for s=1:NodeBelt-1
fprintf(Elements_Rimfile,'%d,%d,%d\n',ELSET_Rim(s),Nodes_Rim(s),Nodes_Rim(s+1));
end
fprintf(Elements_Rimfile,'%d,%d,%d\n',ELSET_Rim(NodeBelt),Nodes_Rim(NodeBelt),Nodes_Rim(1));
fclose(Elements_Rimfile);

%**************************************************************************************************************************************************
%  WRITING THE INPUT FILE   *_*   *_*    *_*    *_*    *_*    *_*    *_* 
%**************************************************************************************************************************************************

%B21 INFLATION & LOADING FILE GENERATION
fileID = fopen('E:/2DTire LoadVsDef/200elem/1_infload_200.inp','w+');

fprintf(fileID,'***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%NODE DEFINITION
fprintf(fileID,'*restart,WRITE,FREQUENCY=1\n');
fprintf(fileID,'**NODE DEFINITION\n');
fprintf(fileID,'*include, input=Nodes_Belt.inp\n');
fprintf(fileID,'*include, input=Nodes_Rim.inp\n');
fprintf(fileID,'*include, input=Nodes_Tread.inp\n');
fprintf(fileID,'*node,nset=rim_reference\n');
fprintf(fileID,'9999999,0,0\n');
fprintf(fileID,'*node,nset=road_reference\n');
fprintf(fileID,'4000000,0,%f\n',road_reference);

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%RIM DEFINITION
fprintf(fileID,'**RIM ELEMENT\n');
fprintf(fileID,'*element,type=R2D2,elset=Elements_Rim\n');
fprintf(fileID,'*include, input=Elements_Rim.inp\n');

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%BELT DEFINITION
fprintf(fileID,'**BELT ELEMENT\n');
fprintf(fileID,'*element,type=B21,elset=Elements_Belt\n');
fprintf(fileID,'*include, input=Elements_Belt.inp\n');

fprintf(fileID,'***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%  CONNECTOR BEHAVIOUR DEFINITION  *_* 

%TREAD CARTESIAN BEHAVIOUR
fprintf(fileID,'*connector behavior, name=Behaviour_TreadCartesian\n');
fprintf(fileID,'*CONNECTOR CONSTITUTIVE REFERENCE\n');
fprintf(fileID,'%f,0,0,0,0,0\n',ltread);
fprintf(fileID,'*connector elasticity, component=1\n');
fprintf(fileID,'%f\n',Kr);
fprintf(fileID,'*connector damping, component=1\n');
fprintf(fileID,'%f\n',Cr);
fprintf(fileID,'*connector stop, component=1\n');
fprintf(fileID,'0,100\n');
fprintf(fileID,'*connector elasticity, component=2\n');
fprintf(fileID,'%f\n',Kt);
fprintf(fileID,'*connector damping, component=2\n');
fprintf(fileID,'%f\n',Ct);
fprintf(fileID,'*connector stop, component=2\n');
fprintf(fileID,'-100,100\n');

fprintf(fileID,'***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%SIDEWALL CARTESIAN TANGENTIAL BEHAVIOUR
fprintf(fileID,'*connector behavior, name=Behaviour_SwallCartesian\n');
fprintf(fileID,'*CONNECTOR CONSTITUTIVE REFERENCE\n');
fprintf(fileID,'%f,0,0,0,0,0\n',sinput.Tire_Belt_Radius-sinput.Tire_Rim_Radius);
fprintf(fileID,'*connector elasticity, component=2\n');
fprintf(fileID,'%f\n',kt);
fprintf(fileID,'*connector stop, component=2\n');
fprintf(fileID,'-100,100\n');

fprintf(fileID,'***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%SIDEWALL RADIAL AXIAL BEHAVIOUR
fprintf(fileID,'*connector behavior, name=Behaviour_SwallRadial\n');
fprintf(fileID,'*CONNECTOR CONSTITUTIVE REFERENCE\n');
fprintf(fileID,'%f,0,0,0,0,0\n',sinput.Tire_Belt_Radius-sinput.Tire_Rim_Radius);
fprintf(fileID,'*connector elasticity, component=1,nonlinear,REGULARIZE=OFF\n');
for o=1:length(SW)
    fprintf(fileID,'%f,%f\n',SW(o,2),SW(o,1));
end
fprintf(fileID,'*connector stop, component=1\n');
fprintf(fileID,'0,%f\n',sinput.Tire_Belt_Radius-sinput.Tire_Rim_Radius+sinput.urpos(end));

%**************************************************************************************************************************************************
%  CONNECTOR ELEMENT DEFINITION *_* 
%**************************************************************************************************************************************************

for i=1:NodeBelt
    
 
%  CONNECTOR LOOP  *_* *_* *_* *_* *_* 

fprintf(fileID,'**TREAD CARTESIAN CONNECTOR ELEMENT No.%d\n',i);
fprintf(fileID,'*element,type=conn2d2, elset=Tread_Cart%d\n',a(i));
fprintf(fileID,'%d,%d,%d\n',a(i),Nodes_Belt(i),Nodes_Tread(i));
fprintf(fileID,'*orientation, name=ori%d, DEFINITION=COORDINATES, SYSTEM=RECTANGULAR\n',a(i));
fprintf(fileID,'%f,%f,%f,%f,%f,%f\n',X_x_r(i),X_y_r(i),X_z_r(i),Y_x_r(i),Y_y_r(i),Y_z_r(i));
fprintf(fileID,'*connector section, elset=Tread_Cart%d, behavior=Behaviour_TreadCartesian\n',a(i));
fprintf(fileID,'cartesian\n');
fprintf(fileID,'ori%d,\n',a(i));

fprintf(fileID,'***********************************************************************************************************************\n');

fprintf(fileID,'**SIDEWALL TANGENTIAL CONNECTOR ELEMENT No.%d\n',i);
fprintf(fileID,'*element,type=conn2d2, elset=SW_cart%d\n',b(i));
fprintf(fileID,'%d,%d,%d\n',b(i),Nodes_Rim(i),Nodes_Belt(i));
fprintf(fileID,'*orientation, name=ori%d, DEFINITION=COORDINATES, SYSTEM=RECTANGULAR\n',b(i));
fprintf(fileID,'%f,%f,%f,%f,%f,%f\n',X_x_r(i),X_y_r(i),X_z_r(i),Y_x_r(i),Y_y_r(i),Y_z_r(i));
fprintf(fileID,'*connector section, elset=SW_cart%d, behavior=Behaviour_SwallCartesian\n',b(i));
fprintf(fileID,'cartesian\n');
fprintf(fileID,'ori%d,\n',b(i));

fprintf(fileID,'***********************************************************************************************************************\n');

fprintf(fileID,'**SIDEWALL RADIAL PIN CONNECTOR ELEMENT No.%d\n',i);
fprintf(fileID,'*element,type=conn2d2, elset=SW_axial%d\n',c(i));
fprintf(fileID,'%d,%d,%d\n',c(i),Nodes_Rim(i),Nodes_Belt(i));
fprintf(fileID,'*orientation, name=ori%d, DEFINITION=COORDINATES, SYSTEM=RECTANGULAR\n',c(i));
fprintf(fileID,'%f,%f,%f,%f,%f,%f\n',X_x_r(i),X_y_r(i),X_z_r(i),Y_x_r(i),Y_y_r(i),Y_z_r(i));
fprintf(fileID,'*connector section, elset=SW_axial%d, behavior=Behaviour_SwallRadial\n',c(i));
fprintf(fileID,'axial\n');
fprintf(fileID,'ori%d,\n',c(i));

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

end

%**************************************************************************************************************************************************
% DEFINING ELSET FOR OUTPUT *_* 
%**************************************************************************************************************************************************

 fprintf(fileID,'*elset,elset=Treadcartesian_All\n'); 
for u=1:NodeBelt
 fprintf(fileID,'Tread_Cart%d\n',a(u));   
end

fprintf(fileID,'***********************************************************************************************************************\n');

fprintf(fileID,'*elset,elset=SWcartesian_All\n'); 
for u=1:NodeBelt
 fprintf(fileID,'SW_cart%d\n',b(u));
end

fprintf(fileID,'***********************************************************************************************************************\n');

fprintf(fileID,'*elset,elset=SWaxial_All\n'); 
for u=1:NodeBelt
 fprintf(fileID,'SW_axial%d\n',c(u));
end
 
fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%RIGID BODY DEFINITION
fprintf(fileID,'**RIGID BODY DEFINITION\n');
fprintf(fileID,'*RIGID BODY, ELSET=Elements_Rim, REF NODE=rim_reference\n');
fprintf(fileID,'*mass,elset=Elements_Rim\n');
fprintf(fileID,'%f\n',mw);
fprintf(fileID,'*rotary inertia,ELSET=Elements_Rim\n');
fprintf(fileID,'%f,%f,%f,%f,%f,%f\n',Ixx,Iyy,Izz,Ixy,Iyz,Izx);

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%MATERIAL DEFINITION
fprintf(fileID,'**MATERIAL DEFINITION\n');
fprintf(fileID,'*Material, name=steel\n');
fprintf(fileID,'*density\n');
fprintf(fileID,'%f\n',rho);%
fprintf(fileID,'*Elastic\n');
fprintf(fileID,'%f,%f\n',E,Poisson);
%BEAM SECTION DEFINITION
fprintf(fileID,'**BEAM SECTION DEFINITION\n');
fprintf(fileID,'*Beam Section, elset=Elements_Belt, material=steel, poisson =%f, section=RECT\n',Poisson);
fprintf(fileID,'%f, %f\n',beam_section_width,beam_section_height);
fprintf(fileID,'0.,0.,-1.\n');

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%ROAD GENERATION
fprintf(fileID,'**ROAD GENERATION\n');
fprintf(fileID,'*surface, type=segments, name=road\n');
fprintf(fileID,'START,-50.0,%f\n',road_reference);
fprintf(fileID,'LINE ,1.0,%f\n',road_reference);
fprintf(fileID,'*rigid body, analytical surface=road, ref node=4000000\n');

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%DEFINING TREAD TIP SURFACE 
fprintf(fileID,'**DEFINING SURFACE INTERACTION\n');
fprintf(fileID,'*SURFACE,NAME=treadsurf,TYPE=NODE\n');
fprintf(fileID,'tread\n');

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%DEFINING BELT SURFACE 
fprintf(fileID,'**DEFINING SURFACE INTERACTION\n');
fprintf(fileID,'*SURFACE,NAME=beltsurf,TYPE=NODE\n');
fprintf(fileID,'belt\n');

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%DEFINING RIM SURFACE 
fprintf(fileID,'**DEFINING SURFACE INTERACTION\n');
fprintf(fileID,'*SURFACE,NAME=rimsurf,TYPE=NODE\n');
fprintf(fileID,'rim\n');

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%DEFINING RADIAL BELT INFLATION SURFACE 
fprintf(fileID,'**DEFINING SURFACE INTERACTION\n');
fprintf(fileID,'*SURFACE,NAME=InfLoad,TYPE=ELEMENT\n');
fprintf(fileID,'Elements_Belt,SPOS\n');

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%CONTACT PAIR GENERATION
fprintf(fileID,'*CONTACT PAIR,INTERACTION=ROUGH1\n');
fprintf(fileID,'treadsurf,road\n');
fprintf(fileID,'*SURFACE INTERACTION,NAME=ROUGH1\n');
fprintf(fileID,'*FRICTION\n');
fprintf(fileID,'0.0\n');
fprintf(fileID,'*CONTACT PAIR,INTERACTION=ROUGH2\n');
fprintf(fileID,'beltsurf,road\n');
fprintf(fileID,'*SURFACE INTERACTION,NAME=ROUGH2\n');
fprintf(fileID,'*FRICTION\n');
fprintf(fileID,'0.0\n');
fprintf(fileID,'*CONTACT PAIR,INTERACTION=ROUGH3\n');
fprintf(fileID,'rimsurf,road\n');
fprintf(fileID,'*SURFACE INTERACTION,NAME=ROUGH3\n');
fprintf(fileID,'*FRICTION\n');
fprintf(fileID,'0.0\n');

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%INFLATION - STANDARD
fprintf(fileID,'*step, nlgeom=yes\n');
fprintf(fileID,'inflation of 240 kpa \n');
fprintf(fileID,'*STATIC\n');
fprintf(fileID,'1E-5,0.1\n');
fprintf(fileID,'*BOUNDARY\n');
fprintf(fileID,'4000000,1,6,0.\n');
fprintf(fileID,'9999999,1,6,0.\n');
fprintf(fileID,'*DSLOAD\n');
fprintf(fileID,'InfLoad,P,%f\n',P_tire);
%FIELD OUTPUT FOR STEP
%*********************************************************
fprintf(fileID,'*OUTPUT, FIELD,variable=preselect\n');
fprintf(fileID,'*ELEMENT OUTPUT,ELSET=Treadcartesian_All\n');
fprintf(fileID,'CTF,CU\n');
fprintf(fileID,'*ELEMENT OUTPUT,ELSET=SWcartesian_All\n');
fprintf(fileID,'CTF,CU\n');
fprintf(fileID,'*ELEMENT OUTPUT,ELSET=SWaxial_All\n');
fprintf(fileID,'CTF,CU\n');
fprintf(fileID,'*node output,nset=road_reference\n');
fprintf(fileID,'U,RF\n');
fprintf(fileID,'*node output,nset=rim_reference\n');
fprintf(fileID,'U,RF\n');
%HISTORY OUTPUT FOR STEP
%*********************************************************
fprintf(fileID,'*output, history,FREQ=1\n');
fprintf(fileID,'*ELEMENT OUTPUT,ELSET=Treadcartesian_All\n');
fprintf(fileID,'CTF,CU\n');
fprintf(fileID,'*ELEMENT OUTPUT,ELSET=SWcartesian_All\n');
fprintf(fileID,'CTF,CU\n');
fprintf(fileID,'*ELEMENT OUTPUT,ELSET=SWaxial_All\n');
fprintf(fileID,'CTF,CU\n');
fprintf(fileID,'*node output,nset=road_reference\n');
fprintf(fileID,'U,RF\n');
fprintf(fileID,'*node output,nset=rim_reference\n');
fprintf(fileID,'U,RF\n');
fprintf(fileID,'*MONITOR,NODE=road_reference,DOF=2\n');
%*********************************************************
fprintf(fileID,'*END STEP\n');

fprintf(fileID,'**\n**\n***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

% LOADING OF TYRE - DISPLACEMENT CONTROLLED
fprintf(fileID,'*amplitude,name=ramp1\n');
fprintf(fileID,'0,0,2,1\n');
fprintf(fileID,'*STEP,NLGEOM=yes,INC=10000\n');
fprintf(fileID,'LOADING OF TYRE - DISPLACEMENT CONTROLLED\n');
fprintf(fileID,'*STATIC\n');
fprintf(fileID,',2.0\n');
fprintf(fileID,'*BOUNDARY,OP=MOD,amplitude=ramp1\n');
fprintf(fileID,'4000000,2,2,%f\n',Tire_Def+.001);
fprintf(fileID,'*END STEP\n');
fprintf(fileID,'***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

%ROAD PULL
% fprintf(fileID,'*amplitude,name=ramp2\n');
% fprintf(fileID,'0,0,2,1\n');
% fprintf(fileID,'*STEP,NLGEOM=yes,INC=200\n');
% fprintf(fileID,'ROAD PULL\n');
% fprintf(fileID,'*STATIC,DIRECT USER CONTROL\n');
% fprintf(fileID,'1E-2,2.0\n');
% fprintf(fileID,'*CHANGE FRICTION,INTERACTION=ROUGH1\n');
% fprintf(fileID,'*FRICTION\n');
% fprintf(fileID,'0.9\n');
% fprintf(fileID,'*BOUNDARY,OP=mod\n');
% fprintf(fileID,'9999999,1,4,0.\n');
% fprintf(fileID,'9999999,6,6,0.\n');
% fprintf(fileID,'4000000,2,2,%f\n',Tire_Def);
% fprintf(fileID,'*BOUNDARY,amplitude=ramp2\n');
% fprintf(fileID,'4000000,1,1,.02\n');
% fprintf(fileID,'*END STEP\n');
% fprintf(fileID,'***********************************************************************************************************************\n**********************************************************************************************************************\n**\n**\n');

% FINDING FREE ROLLING VELOCITY

% fprintf(fileID,'*amplitude,name=ramp2\n');
% fprintf(fileID,'0.0,0.0,1.0,1.0\n');
% fprintf(fileID,'*step, nlgeom=yes,inc=20000\n');
% fprintf(fileID,'Finding FRV,braking\n');
% fprintf(fileID,'*STATIC,direct user control\n');
% fprintf(fileID,'1E-5,2\n');
% fprintf(fileID,'*CHANGE FRICTION,INTERACTION=ROUGH1\n');
% fprintf(fileID,'*FRICTION\n');
% fprintf(fileID,'0.8\n');
% fprintf(fileID,'*BOUNDARY\n');
% fprintf(fileID,'9999999,2,3,0.\n');
% fprintf(fileID,'4000000,3,6,0\n');
% fprintf(fileID,'4000000,2,2,%f\n',Tire_Def);
% fprintf(fileID,'*BOUNDARY,TYPE=VELOCITY,amplitude=ramp2\n');
% fprintf(fileID,'9999999,6,6,31\n');%29.9772   29.6744     29.3716    %29.0688    %28.7660    %28.1604   %27.252   24.224
% fprintf(fileID,'4000000,1,1,11.11\n'); %30.5828   30.8856    31.1884    31.4912      31.794      32.3996      33.0052 34.882  33.308
% fprintf(fileID,'*END STEP\n');
% 
% fprintf(fileID,'*STEP,nlgeom=yes, unsymm=yes,inc=20000\n');
% fprintf(fileID,'Finding FRV,accelerating\n');
% fprintf(fileID,'*DYNAMIC,direct user control\n');
% fprintf(fileID,'1e-5,2\n');
% fprintf(fileID,'*BOUNDARY,OP=MOD\n');
% fprintf(fileID,'9999999,2,3,0.\n');
% fprintf(fileID,'4000000,3,6,0\n');
% fprintf(fileID,'4000000,2,2,%f\n',Tire_Def);
% fprintf(fileID,'*BOUNDARY,TYPE=VELOCITY,amplitude=ramp2\n');
% fprintf(fileID,'9999999,6,6,31.5\n');%29.9772   29.6744     29.3716    %29.0688    %28.7660    %28.1604   %27.252   24.224
% fprintf(fileID,'4000000,1,1,11.11\n'); %30.5828   30.8856    31.1884    31.4912      31.794      32.3996      33.0052 34.882  33.308
% fprintf(fileID,'*END STEP\n');


% 
% fprintf(fileID,'*STEP,INC=300,nlgeom=yes, unsymm=yes\n');
% fprintf(fileID,'test_ACC (40 kmph)\n');
% fprintf(fileID,'*STATIC\n');
% fprintf(fileID,'0.5, 1.0\n');
% fprintf(fileID,'*BOUNDARY,TYPE=VELOCITY\n');
% fprintf(fileID,'9999999,6,6,32\n');
% fprintf(fileID,'4000000,1,1,11.1\n');
% fprintf(fileID,'*END STEP\n');

fclose(fileID);

