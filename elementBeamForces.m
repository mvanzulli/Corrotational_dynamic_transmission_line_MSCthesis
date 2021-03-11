function  [ fs, ks, stress, rotData ]= elementBeamForces( ...
  elemCoords, elemCrossSecParams, elemConstitutiveParams, solutionMethod, Ue, Udote, Udotdote, elemrho ) ;

elemCoords = elemCoords(:)       ;
xs         = elemCoords(1:2:end) ;

booleanCSTangs = 0 ;

% --- material constit params ---
rho = elemrho ;
E   = elemConstitutiveParams(2) ;
nu  = elemConstitutiveParams(3) ;
G   = E/(2*(1+nu)) ;
% -------------------------------

% --- cross section ---
if elemCrossSecParams(1) == 1
    Area = elemCrossSecParams( 2 ) ;
    J    = elemCrossSecParams( 3 ) ;
    Iyy  = elemCrossSecParams( 4 ) ;
    Izz  = elemCrossSecParams( 5 ) ;
    %
    if length( elemCrossSecParams ) > 5
        Jrho =  diag( elemCrossSecParams( 6:8 ) ) ;
    else
        Jrho = rho * diag( [ J Iyy Izz ] ) ;
    end
elseif elemCrossSecParams(1) == 2
    Area = elemCrossSecParams(2)*elemCrossSecParams(3)      ;
    Iyy  = elemCrossSecParams(2)*elemCrossSecParams(3)^3/12 ;
    Izz  = elemCrossSecParams(3)*elemCrossSecParams(2)^3/12 ;
    if elemCrossSecParams(2)==elemCrossSecParams(3)
        J    = 1/3*0.40147*elemCrossSecParams(2)^4 ;
    else
        error('rectangular section type not implemented yet, please create an issue')
    end
    Jrho = rho * diag( [ J Iyy Izz ] ) ;
elseif elemCrossSecParams(1) == 3
    diameter = elemCrossSecParams(2) ;
    Area = pi*diameter^2/4           ;
    Iyy  = pi*diameter^4/64          ;
    Izz  = Iyy                       ;
    J    = Iyy + Izz ;
    Jrho = rho * diag( [ J Iyy Izz ] ) ;
else
  error(' section type not implemented yet, please create an issue')
end
% -------------------------------

%--- Auxiliar matrices ---
I3 = eye(3)     ;
O3 = zeros(3)   ;
O1 = zeros(1,3) ;

permutIndxs = [1:2:5 2:2:6 ([1:2:5]+6) ([2:2:6]+6) ] ;

dg       = Ue      ( permutIndxs ) ;
if solutionMethod > 2
  ddotg    = Udote   ( permutIndxs ) ;
  ddotdotg = Udotdote( permutIndxs ) ;
end

% --- Global thetas ---
tg1 = dg (4:6);
tg2 = dg (10:12);

%% --- Rotation matrices ---
Rg1 = expon( tg1 ) ;
Rg2 = expon( tg2 ) ;

x21 = xs(4:6) - xs(1:3) ;
d21 = dg(7:9) - dg(1:3) ;

lo = sqrt( ( x21       )' * ( x21       ) ) ; %
l  = sqrt( ( x21 + d21 )' * ( x21 + d21 ) ) ; %

% --- rotation matrix to reference configuration ---
Ro = beamRefConfRotMat( x21 ) ;

% --- rigid rotation ---

% --- deformed x axis ---
e1 = ( x21 + d21 ) / l   ;

q1 = Rg1 * Ro * [0 1 0]' ;
q2 = Rg2 * Ro * [0 1 0]' ;
q  = ( q1 + q2 ) / 2     ;

%  ---  deformed z local axis  ---
e3 = cross (e1, q) ;
e3 = e3 / norm(e3) ;

%  ---  deformed y local axis  ---
e2 = cross (e3, e1);

%  --- rotation matrix  ---
Rr = [ e1 e2 e3 ] ;

% --- local displacements ---

%  ---  axial displacement  ---
u  = l - lo;

%  ---  local rotations  ---
Re1 = Rr' * Rg1 * Ro;
Re2 = Rr' * Rg2 * Ro;

tl1 = logar( Re1 ) ;
tl2 = logar( Re2 ) ;

locDisp = [ u tl1' tl2' ] ;

%% --- local force vector and tangent stiffness matrix ---
[fl, kl, strain, stress] = beamLocalStaticForces (u, tl1, tl2, lo, E, G, Area, Iyy, Izz, J ) ;

%% transformation to the new local coordinates
q  = Rr' *  q ;
q1 = Rr' * q1 ;

nu = q(1)/q(2);
nu11 = q1(1)/q(2);
nu12 = q1(2)/q(2);
nu21 = 2*nu-nu11;
nu22 = 2-nu12;

De1 = invTs( tl1 ) ;
De2 = invTs( tl2 ) ;

% --- matrix for transformation between global and relative ---
% --- rotations/moments ---
H  = [  1   O1   O1 ; ...
       O1' De1   O3 ; ...
       O1'  O3  De2 ] ;

fe = H' * fl ;
Dh1 = dinvTs( tl1, fl(2:4) ) * De1 ;
Dh2 = dinvTs( tl2, fl(5:7) ) * De2 ;

Kh = [ 0   O1   O1
      O1' Dh1   O3
      O1'  O3  Dh2 ] ;

ke = H' * kl * H + Kh ;

% --- transformation to the global coordinates ---
r = [ -e1' O1  e1' O1 ]' ;

B = [ r'
   -nu/l*e3' (1-nu12/2)*e1'+nu11/2*e2'  nu/l*e3' 1/2*(-nu22*e1'+nu21*e2')
   -e3'/l e2' e3'/l 0 0 0
    e2'/l e3' -e2'/l 0 0 0
   -nu/l*e3' 1/2*(-nu12*e1'+nu11*e2')  nu/l*e3' (1-nu22/2)*e1'+nu21/2*e2'
   -e3'/l 0 0 0 e3'/l e2'
    e2'/l 0 0 0 -e2'/l e3'];

fg = B' * fe ;

A  = (I3-e1*e1')/l;

Dr=[A  O3 -A  O3
    O3 O3  O3 O3
   -A  O3  A  O3
    O3 O3  O3 O3];

G=[0   0    nu/l  nu12/2  -nu11/2  0  0  0    -nu/l  nu22/2  -nu21/2  0
   0   0    1/l     0        0     0  0  0    -1/l     0        0     0
   0  -1/l  0       0        0     0  0  1/l   0       0        0     0]';

II=[O3 I3 O3 O3
    O3 O3 O3 I3];

P = II - [G'; G'] ;

F = P' * fe(2:7);

sF=[skew(F(1:3))
    skew(F(4:6))
    skew(F(7:9))
    skew(F(10:12))];

EE=[Rr O3 O3 O3
    O3 Rr O3 O3
    O3 O3 Rr O3
    O3 O3 O3 Rr];

nab=[0
    (nu*(fe(2)+fe(5))+fe(3)+fe(6))/l
    (fe(4)+fe(7))/l];

Kg = B' * ke * B + Dr * fe(1) - EE*sF*G'*EE' + EE*G*nab*r' ;


%% --- transformation to the new global coordinates ---

Dg1 = Ts( tg1 ) ;
Dg2 = Ts( tg2 ) ;

q=[fg(1:3)
   Dg1'*fg(4:6)
   fg(7:9)
   Dg2'*fg(10:12)];

Dk1=dTs(tg1,fg(4:6));
Dk2=dTs(tg2,fg(10:12));

H=[I3 O3  O3 O3
   O3 Dg1 O3 O3
   O3 O3  I3 O3
   O3 O3  O3 Dg2];

Kt = H' * Kg * H ;

Kt( 4:6 , 4:6 ) = Kt( 4:6 , 4:6 ) + Dk1 ;
Kt(10:12,10:12) = Kt(10:12,10:12) + Dk2 ;

Kt = (Kt+Kt')/2;

Finte   = zeros(size(q)) ;
dofscomb = [ 1:2:5 2:2:6 7:2:11 8:2:12 ] ;

Finte( dofscomb ) = q ;
KTe = zeros( size(Kt));

KTe( dofscomb, dofscomb ) = Kt ;

fs = {Finte} ;
ks = {KTe};

rotData = {locDisp, Rr} ;

if solutionMethod > 2
%% ---- interpolation functions ---
  % --- linear ---
  N1 = @(x) 1 -x/lo              ;
  N2 = @(x) x/lo				 ;

  % --- cubic ---
  N3 = @(x) x*(1-x/lo)^2		 ;
  N4 = @(x) -(1-x/lo)*(x^2)/lo   ;
  N5 = @(x) (1-3*x/lo)*(1-x/lo)  ;
  N6 = @(x) (3*x/lo-2)*(x/lo)	 ;

  N7 = @(x) N3(x)+N4(x)          ;
  N8 = @(x) N5(x)+N6(x)-1		 ;
  % -------------------------------------

  P1    = @(x) [ 0      0       0  0      0      0 ; ...
                 0      0   N3(x)  0      0  N4(x) ; ...
                 0 -N3(x)	      0  0 -N4(x)      0 ] ;

  ul    = @(x) P1(x) * [ tl1; tl2 ] ; % Eq. 38

  P2    = @(x) [ N1(x)      0      0  N2(x)      0     0 ; ...
                     0  N5(x)      0      0  N6(x)     0 ; ...
                     0      0  N5(x)      0      0 N6(x) ] ;

  N     = @(x) [ N1(x)*I3   O3   N2(x)*I3    O3 ];

  H1    = @(x) N(x) + P1( x ) * P - 1*skew( ul(x) ) * G' ;

  wdoter= G' * EE' * ddotg ;% Eq. 65

  A1    = [   O1          O1   O1       O1 ;
               0 -1    0  O1 	 0  1  0	O1 ;
               0  0 	-1  O1	 0  0  1  O1 ] ;

  udotl = @(x)  P1(x) * P * EE' * ddotg ;

  H1dot = @(x)  N7(x)/(l^2)*A1*(r' * ddotg) - skew( udotl(x) ) * G' ;

  ET = [skew(wdoter)      O3         		O3 			 O3			;
        O3		skew(wdoter)  		O3   		 O3			;
        O3			O3 			skew(wdoter)     O3			;
        O3			O3  			O3       skew(wdoter)   ];

  C1 = @(x)  skew(wdoter)*H1(x) + H1dot(x) -H1(x)*ET;

  udot    = @(x) Rr*H1(x)*EE'*ddotg;
  udotdot = @(x) Rr*H1(x)*EE'*ddotdotg+Rr*C1(x)*EE'*ddotg;

  %% --- Matrix to compute wdot y wdtotdot ---

  H2 = @(x) P2(x)*P+G';

  wdot  = @(x) Rr*H2(x)*EE'*ddotg;


  A2    = [   O1    O1    O1     O1;
            0 0  1  O1  0 0 -1   O1;
            0 -1 0  O1  0 1  0   O1 ] ;

  H2dot    = @(x)	N8(x)/l^2*A2*(r'*ddotg) ;

  C2       = @(x) skew(wdoter)*H2(x) + H2dot(x) - H2(x)*ET ;

  wdotdot  = @(x) Rr*H2(x)*EE'*ddotdotg  + Rr*C2(x)*EE'*ddotg ;

  %% --- Tensor dyadc of Intertia ---
  %compute Rg(x)
  thethaRoof  = @(x) P2(x)*[tl1;tl2] ;
  Rex         = @(x) expon(thethaRoof(x)) ;
  Rgx  	      = @(x) Rr*Rex(x)*Ro';

  Irho		  = @(x) Rgx(x)*Ro*(Jrho)*(Rgx(x)*Ro)';
  Irhoe       = @(x) Rr'*Irho(x)*Rr;

  % --- Compute interial force by quadrature ---
  xIntPoints = [ -sqrt(3/5)     0  sqrt(3/5)  ] ;
  wIntPoints = [        5/9	  8/9        5/9  ] ;

  IntegrandoForce  = @(x)  H1(x)'*Rr'*Area*rho*udotdot(x) ...
                         + H2(x)'*Rr'*( ...
                           Irho(x)*wdotdot(x)...
                           + skew(wdot(x)) * Irho(x) * wdot(x) ...
                         ) ;

 IntegrandoMassMatrix  = @(x) 1*H1(x)'*Area*rho*H1(x)+1*H2(x)'*Irhoe(x)*H2(x);

  %% --- Compute C3 and C4 ---

  h1 = @(x) H1(x) * ddotg ;
  h2 = @(x) H2(x) * ddotg ;

  rElem = [ [-1 0 0]   O1  [1 0 0] O1];

  F1    = [skew(udot(0))' skew(wdot(0))' skew(udot(lo))' skew(wdot(lo))']';

  C3  = @(x) -skew(h1(x))*G'  + (N7(x)/l^2)*A1*(ddotg*rElem)...
                +skew(wdoter)*P1(x)*P + H1(x)*F1*G';

  C4  = @(x) -skew(h2(x))*G' + (N8(x)/l^2)*A2*ddotg*rElem + H2(x)*F1*G';

  %% --- Compute Gyroscopic Matrix---
  IntegrandoGyroMatrix  = @(x)  H2(x)' * ( ( skew(wdoter) * Irhoe(x) ) - skew( Irhoe(x) * wdoter) ) * H2(x) ...
                              + H1(x)' * Area*rho*(C1(x) + C3(x))  + H2(x)'*Irhoe(x)*(C2(x)+C4(x)) ;

  sumForce = zeros (12, 1 ) ;
  sumGyro  = zeros (12    ) ;
  sumMass  = zeros (12    ) ;


  for ind = 1 : length( xIntPoints )
    sumForce = sumForce ...
      + lo/2 * wIntPoints( ind ) * IntegrandoForce     ( lo/2 * (xIntPoints( ind ) + 1) ) ;
    %
    sumGyro = sumGyro ...
      + lo/2 * wIntPoints( ind ) * IntegrandoGyroMatrix( lo/2 * (xIntPoints( ind ) + 1) ) ;
    %
    sumMass = sumMass ...
      + lo/2 * wIntPoints( ind ) * IntegrandoMassMatrix( lo/2 * (xIntPoints( ind ) + 1) ) ;
  end

  Fine       = EE * sumForce      ;
  GyroMatrix = EE * sumGyro * EE' ;
  MassMatrix = EE * sumMass * EE' ;

  %% --- Add Bt Matrix ---

  Bt=[I3   O3       O3      O3
      O3 inv(Dg1)'    O3      O3
      O3     O3      I3      O3
      O3     O3      O3      inv(Dg2)' ];
  MassMatrix =MassMatrix*Bt ;
  GyroMatrix =GyroMatrix*Bt ;
  %% --- Switch base ---
  Fine       = Cambio_Base(Fine);
  GyroMatrix = Cambio_Base(GyroMatrix);
  MassMatrix = Cambio_Base(MassMatrix);

  fs{3} = Fine ;
  ks{2} = GyroMatrix ;
  ks{3} = MassMatrix ;

end
