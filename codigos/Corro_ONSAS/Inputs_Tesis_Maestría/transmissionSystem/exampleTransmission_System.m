%% Transimission System input
%
%%
clear all, close all
dirOnsas             = [ pwd '\..\..\..\..\ONSAS\ONSAS' ]                  ;
dirTowerGeometry     = [pwd '/../..\auxTowerGeometry']                  ;
dirCableGeometry     = [pwd '/../..\auxCableGeometry_TSystem']          ;
dirFuerzaViento      = [pwd '/../..\Fuerza_Viento']                     ; 
addpath (dirOnsas,dirCableGeometry,dirTowerGeometry,dirFuerzaViento)    ;
problemName = 'Transmission_Towers'                                     ;

%% Pricnipal Varaible Params ================= %=
NelemC     = 200  ;                              %=
finalTime  = 1   ;                             %= 
timeIncr   = 0.5 ;                            %=
% ===============================================
%% Towers Geometry

%====================== Tower 1 ======================

%Nodes Matrix:
Nodes_T1           = load('MatNodos.txt');
Num_nodes_T1       = size (Nodes_T1,1);

%Conec Matrix:
%Conec elem
Conec_elem_T1      = load('Conec_elem.txt');
Cambio_areas_T1    = load('Vec_sec.txt');
Conec_elem_T1(:,4) = Cambio_areas_T1;
Num_elem_T1        = size(Conec_elem_T1,1);

%Conec nod
Conec_Nod_T1       = load('Conec_Nod.txt');


%====================== Tower 2 ======================
%% Distance
DistVaneX_1 = 0 ;
DistVaneY_1 = 267.09 ;
DistVaneZ_1 = 0 ;
%====================== Tower 2 ======================
%% Nodes Matrix:
Nodes_T2           = [ Nodes_T1(:,1)+DistVaneX_1 , Nodes_T1(:,2)+DistVaneY_1 ,...
                                                     Nodes_T1(:,3)+DistVaneZ_1] ;
                                                 
Num_nodes_T2       = Num_nodes_T1;

%Conec Matrix:
%Conec elem
Conec_elem_T2      = load('Conec_elem.txt');
Conec_elem_T2(:,6) = Conec_elem_T2(:,6)+ Num_nodes_T1 ;
Conec_elem_T2(:,7) = Conec_elem_T2(:,7)+ Num_nodes_T1 ;
Num_elem_T2        = Num_elem_T1 ;

%Conec nod
Conec_Nod_T2       = Conec_Nod_T1 ;
Conec_Nod_T2(:,6) = Conec_Nod_T2(:,6) + Conec_Nod_T1(end,6) ;

%====================== Tower 3 ======================
%% Distance
DistVaneX_2 = 0 ;
DistVaneY_2 = 267.09 ;
DistVaneZ_2 = 0 ;
%====================== Tower 3 ======================
%% Nodes Matrix:
Nodes_T3           = [ Nodes_T2(:,1)+DistVaneX_2 , Nodes_T2(:,2)+DistVaneY_2 ,...
                                                     Nodes_T2(:,3)+DistVaneZ_2] ;
                                                 
Num_nodes_T3       = Num_nodes_T1;

%Conec Matrix:
%Conec elem
Conec_elem_T3      = load('Conec_elem.txt');
Conec_elem_T3(:,6) = Conec_elem_T3(:,6)+ Num_nodes_T1*2 ;
Conec_elem_T3(:,7) = Conec_elem_T3(:,7)+ Num_nodes_T1*2 ;
Num_elem_T3        = Num_elem_T1 ;

%Conec nod
Conec_Nod_T3       = Conec_Nod_T1 ;
Conec_Nod_T3(:,6) = Conec_Nod_T3(:,6) + Conec_Nod_T1(end,6)*2 ;

%====================== Joint Towers  ======================

Nodes_T      = [Nodes_T1  ;Nodes_T2 ;Nodes_T3 ]        ;
Num_nodes_t  = Num_nodes_T1+Num_nodes_T2 + Num_nodes_T3    ;
Conec_Nod_T  = [Conec_Nod_T1;Conec_Nod_T2;Conec_Nod_T3]    ;
Conec_elem_T = [Conec_elem_T1;Conec_elem_T2;Conec_elem_T3]  ;
Num_elem     = Num_elem_T1 +Num_elem_T2 + Num_elem_T3       ;


%Nodos de max x (de cadenas)
%Torre 1 arriba: 157 158    medio:111 112         abajo 77 78
%Torre 2 arriba: 349 350    medio:303 304         abajo 269 270

% nodesT1 = [157 158 111 112 77  78] ;
% nodesT2 = [349 350 303 304 269 270];
nodesT1 = [157 ] ;
nodesT2 = [349 ];
nodesT3 = nodesT2 + 192            ;

% ====================== Cable and Chain  ======================
%% Geometrical properties
%conductor
dc  = 0.0281*1 ;
lc1      = norm([DistVaneX_1 DistVaneY_1 DistVaneZ_1 ]);
lc2      = norm([DistVaneX_1 DistVaneY_1 DistVaneZ_1 ]);



%chain
da      = 0.1 ; la  = 3;
NelemA  = 1   ;


%mat and sec params position cable
matPositionC = 2   ;
secPostionC  = 8   ;

%mat and sec params position chain
matPositionA = 3   ;
secPostionA  = 9   ;

%% Insert cables
%initialize conectivity and nodes matrix
Conec_Nod  = Conec_Nod_T ;
Conec_Elem = Conec_elem_T;
Nodes      = Nodes_T     ;
nonHomogeneousInitialCondU0 = [];
%between first and second tower

indexVane= 1 ;
for i=1:size(nodesT1,2)
Node1 = nodesT1(i);
Node2 = nodesT2(i);
[Conec_Nod,Conec_Elem,Nodes,nonHomogeneousInitialCondU0] = cableGeometry(...
                                  NelemC,NelemA,la,Nodes,Node1,Node2,...
            matPositionC,secPostionC,matPositionA,secPostionA,...
                                                    Conec_Nod,Conec_Elem,nonHomogeneousInitialCondU0,indexVane) ;
end
%betwwen second and third tower
%Add nonHomogeneousInitialCondU0
indexVane=2 ;
for i=1:size(nodesT2,2)
Node1 = nodesT2(i);
Node2 = nodesT3(i);
[Conec_Nod,Conec_Elem,Nodes,nonHomogeneousInitialCondU0] = cableGeometry(...
                                  NelemC,NelemA,la,Nodes,Node1,Node2,...
            matPositionC,secPostionC,matPositionA,secPostionA,...
                                                    Conec_Nod,Conec_Elem,nonHomogeneousInitialCondU0,indexVane) ;
end


% nonHomogeneousInitialCondU0 = [];
%Convert 2 cell
Conec_Nod_cell  = mat2cell(Conec_Nod,ones(1,size(Conec_Nod,1)));
Conec_elem_cell = mat2cell(Conec_Elem,ones(1,size(Conec_Elem,1)));

Conec           = [Conec_Nod_cell;Conec_elem_cell];



% ======================================================================
%% --- MELCS parameters ---
materialsParams = cell(1,1) ; % M
elementsParams  = cell(1,1) ; % E
loadsParams     = cell(1,1) ; % L
crossSecsParams = cell(1,1) ; % C
springsParams   = cell(1,1) ; % S

%Tower params
Et = 300e9 ;  nut = 0.3  ;  rhot = 7850 ;
%Cable params
%--------------------
% DRAKE ACSR 7/26 2017 Foti
EA   = 29700e3    ;
GJ   = 159        ;
EI   = 2100       ;
massPerMeter = 1.8;
A    = pi*dc^2/4  ;
Ec   = EA/A       ;
J    = pi*dc^4/64 ;
Gc   = GJ/J       ;
rhoc = massPerMeter/A;
nuc  = Ec/(2*Gc)-1;

%Isolator chain params
Ea = 72e9  ;  nua  = 0.3 ;  rhoa = 2555 ;
						% 3 green 2 ing (tipo de deformacion unitaria)
materialsParams = {[ rhot 3 Et nut];[rhoc 1 Ec nuc ] ;[rhoa 3 Ea nua ]} ; 

elementsParams = { 1; [ 2 0 ];3 } ; 

%Constant Loads

A1 = 5100e-6;
d1 = sqrt (4*A1/pi);
A2 = 1790e-6;
d2 = sqrt (4*A2/pi);
A3 = 1790e-6;
d3 = sqrt (4*A3/pi);
A4 = 2270e-6;
d4 = sqrt (4*A4/pi);
A5 = 1790e-6;
d5 = sqrt (4*A5/pi);
A6 = 1790e-6;
d6 = sqrt (4*A6/pi);
A7 = 1790e-6;
d7 = sqrt (4*A7/pi);


crossSecsParams = { [3 d1];[3 d2];[3 d3];[3 d4];[3 d5];[3 d6];[3 d7];[3 dc];[3 da] } ;

springsParams = { [ inf  0  inf  0  inf   0 ];[ 0  inf  inf  0  0  inf ]} ;


% ======================================================================
% % Method Params
% 
% %Static
% controlDofs = [ 183 5 -1 ] ;
% 
% % analysis parameters
% stopTolIts       = 100       ;
% stopTolDeltau    = 1.0e-5   ;
% stopTolForces    = 1.0e-5   ;
% 
% %Loads
% targetLoadFactrNR  = 0      ;
% nLoadSteps         = 1     ;
% 						
% numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
%                              targetLoadFactrNR nLoadSteps ] ; 


%% Dynamic
% times


% analysis parameters
stopTolIts       = 30       ;
stopTolDeltau    = 1.0e-5   ;
stopTolForces    = 1.0e-5   ;

alphaHHT = -0.05 ;
numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;

% ======================================================================
%% Booleans
reportBoolean = 1 ;
storeBoolean = 1 ;
stabilityAnalysisBoolean = 0 ;
% ================================================================= 
%% Loads and supports:
%Loads---------------------
booleanSelfWeightZ = 1 ;
loadsParams   = {[ 1 1   1 0 0 0 0 0 ]} ;
%compute Fviento
DensidadAire = 1.2 ;
C_d = 1.5          ; 
Plot_ft = 0        ;
factAmpVel= 2    ;
[F_t,v_sim]=Fuerza_Viento (timeIncr,finalTime, DensidadAire, C_d, dc,lc1,...
                                            NelemC,Plot_ft,factAmpVel);
vmedia = mean(v_sim) ;

%Dampint
factorAmplDamping = 0.6;

 %nodalDispDamping = DensidadAire*C_d*dc*lc1/NelemC*vmedia/2  *  factorAmplDamping ;

Tinicio             = 100 ;
Tfinal              = finalTime;
loadFactorsFunc     = @(t) F_t(round(t/timeIncr)+1)*(t>Tinicio)*(t<Tfinal);

%  loadFactorsFunc = @(t) 6*t


% userLoadsFilename = 'myTransmissionSystemForce'

% ================================================================= 
%% Ploteo
plotParamsVector = [0 ];
plotParamsVector = [ 3 200 ];
printFlag = 0 ;

% Plot_Tower
ONSAS

save (datosTransmissionSystem)
												



