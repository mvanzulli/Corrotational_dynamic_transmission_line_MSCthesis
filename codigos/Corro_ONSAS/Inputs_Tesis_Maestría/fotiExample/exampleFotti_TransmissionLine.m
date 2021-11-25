% ------------------------------------------------------------------------------
% example TransmissionLine
%
% ------------------------------------------------------------------------------

clear all, close all
dirOnsas = [ pwd '\..\..\..\..\ONSAS\ONSAS' ]          ;
dirCableGeometry = [pwd '/../..\auxCableGeometry']  ;
addpath (dirOnsas,dirCableGeometry)                 ;
problemName = 'Fotti_TransmissionLine'              ;



%% Geometric properties cable:
%-------------------- 
dc  = 0.0281*1 ; lc  = 267.09;
% Ac  = pi*dc^2/4 ;
% Iyc = pi*dc^4/64;
% Izc = Iyc        ; 
% Jc  = 2*Iyc     ;

%% Geometric properties isolator chain:
%-------------------- 
da  = 0.1 ; la  = 3;
% Aa  = pi*da^2/4  ;
% Iya = pi*da^4/64 ;
% Iza = Iya        ; 
% Ja  = 2*Iya      ;

%% Material properties cable :
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

%% Material properties isolator chain :
%--------------------
Ea   = 72e9    ;
nua  = 0.3      ; 
rhoa = 2555	    ;
Ga	= Ea/(2*(1+nua));	

%% Elem CrosSec and elem Params Cells:
%--------------------
%1 cable 2 isolator chain
materialsParams = {[rhoc 1 Ec nuc ] ;[rhoa 3 Ea nua ]} ;

crossSecsParams = {[3 dc] ; [3 da ]} ;	
% crossSecsParams = {[3 dc] ; [3 da]} ;	
elementsParams  = { 1;[2 1]; 3} ;

%% Nodes and Conec matrix:
%-----------------------

%Nelem Cable
NelemC  =150;
NnodesC = NelemC+1;

%Nelem Aisladora
NelemA  = 1;
NnodesA = NelemA+1;

%Nelem Tottales
Nelem = NelemC+NelemA*2	;
Nnodes = Nelem+1;


 Nodes1 = [ 	 zeros(NelemA,1)    zeros(NelemA,1) 		         ((NelemA):-1:1)'*la/NelemA ;
 			zeros(NnodesC,1)  (0:(NelemC))'*lc/NelemC            zeros(NnodesC,1) 	   ;	 
 			zeros(NelemA,1)     lc*ones(NelemA,1) 		          (1:(NelemA))'*la/NelemA ] ;
        
% Nodes1 = [ 	 zeros(NnodesC,1)  (0:(NelemC))'*lc/NelemC            zeros(NnodesC,1) 	  ] ;
													   
% 				   %Material             %Element                   %Loads                 %CrossSection             %Springs 					 
% auxConecElem  =[ [(ones(NelemA,1)*2 )    (ones(NelemA,1)*2)      (zeros(NelemA,1))         (ones(NelemA,1)*2)      (zeros(NelemA,1)) ...  %ElemNodes...
%                 (1:(NelemA))'                               (2:(NnodesA))'                              zeros(NelemA,2)];  
%                  [(ones(NelemC,1)*1 )    (ones(NelemC,1)*3)      (zeros(NelemC,1))         (ones(NelemC,1)*1)      (zeros(NelemC,1)) ...  %ElemNodes...
%                 (NelemA+1:(NelemA+NelemC))'                 (NelemA+2:NelemA+1+NelemC)'                  zeros(NelemC,2)];
%                  [(ones(NelemA,1)*2 )    (ones(NelemA,1)*2)      (zeros(NelemA,1))         (ones(NelemA,1)*2)      (zeros(NelemA,1)) ...  %ElemNodes...
%                 (NelemA+NelemC+1:NelemA+NelemC+NelemA)'	 	(NelemA+NelemC+2:NelemA+NelemC+NelemA+1)']   zeros(NelemA,2)  ] ;  
%   
% auxConecElem  =[ [(ones(NelemC,1)*1 )    (ones(NelemC,1)*3)      (zeros(NelemC,1))         (ones(NelemC,1)*1)      (zeros(NelemC,1)) ...  %ElemNodes...
%                 (NelemA+1:(NelemA+NelemC))'                 (NelemA+2:NelemA+1+NelemC)'                  zeros(NelemC,2)] ] ;  
%   
 auxConecNodes =[ (zeros(Nnodes,1)*1 )    (ones(Nnodes,1))         (zeros(Nnodes,1))         (zeros(Nnodes,1))        (zeros(Nnodes,1))    (1:Nnodes)'] ;
% 
% % Fixed nodes
% auxConecNodes(1,5)              = 1 ;
% auxConecNodes(end,5)            = 1 ;
% auxConecNodes(NnodesA,5)        = 2 ;
% auxConecNodes(Nnodes-NnodesA+1,5) = 2 ;


 %%Loades  Geometry using cableGeometry script
%mat and sec params position cable
matPositionC = 1   ;
secPostionC  = 1   ;

%mat and sec params position chain
matPositionA = 2   ;
secPostionA  = 2   ;

Node1 = 1       ;
Node2 = 2  ;
Nodesaux = [Nodes1(1,:);Nodes1(end,:)];
nodesT1 = [Node1;Node2];

% cd 
% cd ..\auxCableGeometry

indexVane = 2;
Conec_Nod = [ [0     1     1     0     1     1];...
              [0     1     1     0     1     2]];
Conec_Elem=[];
nonHomogeneousInitialCondU0 = [] ;
[Conec_Nod,Conec_Elem,Nodes,nonHomogeneousInitialCondU0] = cableGeometry(...
                                      NelemC,NelemA,la,Nodesaux,Node1,Node2,...
                       matPositionC,secPostionC,matPositionA,secPostionA,...
               Conec_Nod,Conec_Elem,nonHomogeneousInitialCondU0,indexVane) ;
% cd ..

% nonHomogeneousInitialCondU0 = []
%Build Conec Cell
Conec = cell (Nelem+Nnodes,1);
for i = 1:Nnodes
    Conec{i,1}=Conec_Nod(i,:) ;
end
for i =  1:Nelem
    Conec{Nnodes+i,1} =Conec_Elem(i,:);
 end
% Nodes = Nodes1 ;
% Conec = cell (Nelem+Nnodes,1);
% for i = 1:Nnodes
%     Conec{i,1}=auxConecNodes(i,:) ;
% end
% for i =  1:Nelem
%     Conec{Nnodes+i,1} = auxConecElem(i,:);
% end
%% Suports---------------------
springsParams    = {[ inf  0  inf  0  inf  0 ];[ inf  inf  inf  inf  inf  inf ]} ;
% springsParams    = {[ inf  0  inf  0  inf  0 ];[ inf  0  inf  0  inf  0 ]} ;

%% Parameters and tolerances:
%-----------------------
% times
finalTime  = 350    ;
timeIncr   =  1  ;

nLoadSteps = finalTime/timeIncr ;

% tolerances
stopTolDeltau = 1e-5        ; 
stopTolForces = 1e-7        ;
stopTolIts    = 30          ;

%method params
%Newmark
%~ DeltaNW    =  0.5               ;
%~ AlphaNW    =  0.25              ;
%~ numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

%HHT
alphaHHT = -0.1 ;
numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;
    
%% Loads and supports:
%Loads---------------------
booleanSelfWeightZ = 1 ;
loadsParams   = {[ 1 1   1 0 0 0 0 0 ]} ;
%compute Fviento
DensidadAire = 1.2 ;
C_d = 1.5          ; 
vmedia = 10 ;

% Dampint
 factorAmplDamping = 1;
% 
%   nodalDispDamping = DensidadAire*C_d*dc*lc/NelemC*vmedia/2  *  factorAmplDamping ;
Tinicio = 50;
Tfinal  = finalTime;
factorAmplForce = 1;


W = @(t)  20*(t-Tinicio)/(finalTime-Tinicio)*(t>Tinicio) ;
loadFactorsFunc =  @(t) 0.5 * DensidadAire *W(t)^2 *C_d * dc*lc/(NelemC);


%% Booleans control and plot:
%-----------------------
%controlDof
controlDofs = [ Nelem/2 5 1 ] ;
%Plots
plotParamsVector = [ 3 100 ];
printFlag = 0 ;
%Booleans
flagOutputMatrices = 0;
storeBoolean = 1  ;
reportBoolean = 1 ;

ONSAS

save (datosFotiExample)