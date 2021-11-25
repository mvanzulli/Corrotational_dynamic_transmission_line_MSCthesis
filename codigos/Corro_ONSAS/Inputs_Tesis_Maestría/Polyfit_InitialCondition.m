% ------------------------------------------------------------------------------
% Polifyt Initia Condition TransmissionLine
%
% ------------------------------------------------------------------------------
%Initial cable Node 
Node1 = NelemA+2              ;
dofZ1  = (NelemA+1)*6-1        ;
%Final cable Node
Node2 = NelemA+NelemC+1       ;
dofZ2  = (NelemA+NelemC+1)*6-1 ;


dispZConverged = matUs(dofZ1:6:dofZ2,end)   ;
VecX           = (0 :lc/NelemC :lc)'     ;
polinomioDisp      = polyfit(VecX,dispZConverged,5);


Z0             = polyval(polinomioDisp,VecX) ;

Error = Z0 - dispZConverged ;

figure (1)
hold on
plot (VecX,dispZConverged,'r')
plot (VecX,Z0,'k')
