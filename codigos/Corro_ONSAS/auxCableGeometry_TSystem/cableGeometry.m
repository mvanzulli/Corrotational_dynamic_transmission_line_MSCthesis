 %Create matrix cable

function[Conec_Nod,Conec_Elem,Nodes,nonHomogeneousInitialCondU0] = cableGeometry(...
             NelemC,NelemA,la,Nodes_T,Node1,Node2,...
            matPositionC,secPostionC,matPositionA,secPostionA,...
            Conec_Nod_T,Conec_Elem_T,nonHomogeneousInitialCondU0,indexVane) 
%Cable length
 lc = norm(Nodes_T(Node1,:)-Nodes_T(Node2,:));
                                                                             
%Nnodes
NnodesC = NelemC+1;
NnodesA = NelemA+1;

%Nelem Tottales
Nelem = NelemC+NelemA*2	;
Nnodes = Nelem+1;

% El nodo inicial y el nodo final ya estan definidos
Nodes_Cable     = [ Nodes_T(Node1,1)*ones(NelemA-1,1)                           Nodes_T(Node1,2)*ones(NelemA-1,1)                         (Nodes_T(Node1,3)-(((NelemA-1):-1:1)'*la/NelemA)) ;
                    (linspace(Nodes_T(Node1,1),Nodes_T(Node2,1),NelemC+1))'  (linspace(Nodes_T(Node1,2),Nodes_T(Node2,2),NelemC+1))'       ((linspace(Nodes_T(Node1,3),Nodes_T(Node2,3),NelemC+1)-la)');	 
                    Nodes_T(Node2,1)*ones(NelemA-1,1)                   Nodes_T(Node2,2)*ones(NelemA-1,1)                              ( Nodes_T(Node2,3)-(((NelemA-1):-1:1)'*la/NelemA))] ;
Nnodes_Cable    = size(Nodes_Cable,1);

%Initial condition
% MaxZ     = 8;
% 
% B        = -MaxZ ;
% A        = -B*4/lc^2;
% Catenaria = @(x) A*(x-lc/2).^2+B;
% uZcables      = Catenaria(Nodes_Cable(:,2)-Nodes_Cable(1,2));
% Nodes_Cable(:,3)=Nodes_Cable(:,3)+uZcables ;
%     
%% Add to Nodes matrix
if indexVane==2
  Nodes_Cable(1,:) = [] ;  
end
Nodes           =   [ Nodes_T ; Nodes_Cable ]   ;
NumNodesTowers  =   size(Nodes_T,1)             ;


                        %Material                           %Element                   %Loads                                %CrossSection                %Springs 					 
auxConecElem_Cable  =[ [(ones(NelemA-1,1)*matPositionA )    (ones(NelemA-1,1)*2)      (zeros(NelemA-1,1))         (ones(NelemA-1,1)*secPostionA)      (zeros(NelemA-1,1)) ...  %ElemNodes...
                      (NumNodesTowers+1:(NumNodesTowers+NelemA-1))'          (NumNodesTowers+2:(NumNodesTowers+NelemA))'  ];  
                       [(ones(NelemC,1)*matPositionC )    (ones(NelemC,1)*3)      (zeros(NelemC,1))         (ones(NelemC,1)*secPostionC)      (zeros(NelemC,1)) ...  %ElemNodes...
                       (NumNodesTowers+NelemA-1:(NumNodesTowers+NelemA+NelemC-1)-1)'                 (NumNodesTowers+NelemA+1-1:NumNodesTowers+NelemA+NelemC-1)'                 ];
                        [(ones(NelemA-1,1)*matPositionA )    (ones(NelemA-1,1)*2)      (zeros(NelemA-1,1))         (ones(NelemA-1,1)*secPostionA)      (zeros(NelemA-1,1)) ...  %ElemNodes...
                       (NumNodesTowers+NelemA+NelemC+1:NumNodesTowers+NelemA+NelemC)'	 	(NumNodesTowers+NelemA+NelemC+2:NumNodesTowers+NelemA+NelemC+1)'   ] ] ;  

%add tower nodes to conect                   
auxConecElem_Cable = [ [matPositionA                              2                       0                     secPostionA             0 ...
                            Node1                           NumNodesTowers+1        ];
                                auxConecElem_Cable;
                       [matPositionA                              2                       0                     secPostionA             0    ...
                            size(Nodes,1)                       Node2               ] ];

if indexVane==1                  
        auxConecNodes     =[ (zeros(Nnodes_Cable,1) )    (ones(Nnodes_Cable,1))         (ones(Nnodes_Cable,1))         (zeros(Nnodes_Cable,1))        (zeros(Nnodes_Cable,1))    (size(Nodes_T,1)+1:size(Nodes,1))'] ;
        auxConecElem_Cable (2,:)=[];
elseif indexVane==2
        auxConecNodes      =[ (zeros(Nnodes_Cable-1,1))    (ones(Nnodes_Cable-1,1))         (zeros(Nnodes_Cable-1,1))         (zeros(Nnodes_Cable-1,1))        (zeros(Nnodes_Cable-1,1))    (size(Nodes_T,1)+1:size(Nodes,1))'] ;
        cadenaAMofif       = auxConecElem_Cable (1,:);
        cadenaAMofif(end)  = cadenaAMofif(end)-1;
        cadenaAMofif(end-1)= cadenaAMofif(end)-1;
        cadenaAMofif (1:end-2)=[2 3 0 8 0]  ;
        auxConecElem_Cable (1,:) = cadenaAMofif       ;

% auxConecNodes      =[ (zeros(Nnodes_Cable,1)*1 )    (ones(Nnodes_Cable,1))         (zeros(Nnodes_Cable,1))         (zeros(Nnodes_Cable,1))        (zeros(Nnodes_Cable,1))    (size(Nodes_T,1)+1:size(Nodes,1))'] ;  
end

%Add springs
auxConecNodes(1,5)   = 2;
auxConecNodes(end,5) = 2;
%Add conecElem conecNodes 
Conec_Nod  = [Conec_Nod_T;auxConecNodes];
Conec_Elem = [Conec_Elem_T;auxConecElem_Cable];



%% Initial condition
MaxZ     = 5.2; %foti
% MaxZ     =  1; %Largo 100%
Coef     = 10e-3 ;
B        = -MaxZ ;
A        = -B*4/lc^2;
% Catenaria = @(x) A*(x-lc/2).^2+B;
 Catenaria= @(x) cosh(Coef*(x-lc/2))-cosh(-Coef*lc/2);
 Catenaria=@(x)-Catenaria(x)*B/-Catenaria(lc/2);
Derivada = @(x)B*Coef*sinh(Coef*(x-lc/2));

% Derivada  = @(x) 2*A*(x-lc/2)   ;
uZcables      = Catenaria(Nodes_Cable(:,2)-Nodes_Cable(1,2));
thetaYcables  =     Derivada (Nodes_Cable(:,2)-Nodes_Cable(1,2)); 
nodesLocationCable = size(Nodes,1)-size(Nodes_Cable,1);
nodesCable         = nodesLocationCable+1:size(Nodes,1);
% nonHomogeneousInitialCondU0 = [nonHomogeneousInitialCondU0;...
%                                nodesCable' 5*ones(size(nodesCable,2),1) uZcables;
%                                nodesCable' 1*ones(size(nodesCable,2),1) thetaYcables];
     nonHomogeneousInitialCondU0 = [nonHomogeneousInitialCondU0;...
                               nodesCable' 5*ones(size(nodesCable,2),1) uZcables];
nodesCableVane2 = [];                     
if indexVane ==2
nodesCableVane2  = [nodesCableVane2;nodesCable] ;
node1CablesVane2 = nodesCableVane2(1)    ;       
nodEndCableVane2 = nodesCableVane2(end) ;        
end

end