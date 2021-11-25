%% Transimission Towers input
%
%%
clear all, close all
dirOnsas = [ pwd '\..\..\..\..\ONSAS\ONSAS\src' ] ;
dirTowerGeometry = [pwd '/../..\auxTowerGeometry'] ;
addpath (dirOnsas,dirTowerGeometry)  ;
problemName = 'Transmission_Towers'  ;
% ======================================================================
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
DistVaneY_1 = 300 ;
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
DistVaneY_2 = 300 ;
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

Nodes     = [Nodes_T1  ;Nodes_T2 ;Nodes_T3 ]             ;
Num_nodes_t  = Num_nodes_T1+Num_nodes_T2 + Num_nodes_T3     ;
Conec_Nod_T  = [Conec_Nod_T1;Conec_Nod_T2;Conec_Nod_T3]     ;
Conec_elem_T = [Conec_elem_T1;Conec_elem_T2;Conec_elem_T3]  ;
Num_elem     = Num_elem_T1 +Num_elem_T2 + Num_elem_T3       ;

%Convert 2 cell
Conec_Nod_cell  = mat2cell(Conec_Nod_T,[ones(1,size(Conec_Nod_T,1))]);
Conec_elem_cell = mat2cell(Conec_elem_T,[ones(1,size(Conec_elem_T,1))]);
Conec           = [Conec_Nod_cell;Conec_elem_cell];%%

cd ..
% =====================================================================
%% --- MELCS parameters ---
materialsParams = cell(1,1) ; % M
elementsParams  = cell(1,1) ; % E
loadsParams     = cell(1,1) ; % L
crossSecsParams = cell(1,1) ; % C
springsParams   = cell(1,1) ; % S

E = 210e9 ;  nu = 0.3 ;  rho = 7850 ;
						% 3 green 2 ing (tipo de deformacion unitaria)
materialsParams = {[ rho 3 E nu ]} ;

elementsParams = { 1; [ 2 1 ]} ; 

loadsParams  =     {[ 1 1  1  0  0  0  0  0 ]};
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


crossSecsParams = { [3 d1];[3 d2];[3 d3];[3 d4];[3 d5];[3 d6];[3 d7] } ;

springsParams = { [ inf  0  inf  0  inf   0 ]} ;


% ======================================================================
%% Method Params

%Static

controlDofs = [ 183 5 -1 ] ;

% analysis parameters
stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-8 ;
stopTolForces    = 1.0e-8  ;

%Loads
booleanSelfWeightZ = 1 ;
targetLoadFactrNR  = 1e3    ;
nLoadSteps         = 10    ;
						
loadFactorsFunc = @(t) 0;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                             targetLoadFactrNR nLoadSteps ] ; 

% ======================================================================
%% Booleans
reportBoolean = 1 ;
storeBoolean = 1 ;
stabilityAnalysisBoolean = 0 ;

% ================================================================= 
%% Ploteo
plotParamsVector = [0 ];
plotParamsVector = [ 3 300 ]; 
printFlag = 0 ;

ONSAS



