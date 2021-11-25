
function [F_t,v_sim] = Fuerza_Viento (timeIncr,finalTime, DensidadAire, C_d, dc,lc,NelemC,Plot_ft,factAmpVel)
%% Velocidad de viento

viento = load('Viento_tormenta.txt');
V_0 = find(viento(:,2) == 0);
viento(V_0,:) = [] ;



Tiempo = viento(:,1);
Tiempo = Tiempo(:) - Tiempo(1);
Velocidad = viento(:,2);
Velocidad(1) = 0;
Adicional = 0;
Velocidad = Velocidad + Adicional; % 5 m/s mas para contemplar la componente fluctuante 

Velocidad_2 = viento(:,3);
Velocidad(1) = 0;

% Tiempo y velocidades para la simulacion

dt = timeIncr;

t_sim = (Tiempo(1) : dt : Tiempo(end))' ;

v_sim = zeros(length(t_sim),1);
v_sim(1) = Velocidad(1);
v_sim(end) = Velocidad(end);

for i = 2 : length(t_sim)-1
    
    V_v = find( Tiempo >  t_sim(i) );
    V_1 = Velocidad(V_v(1)-1);
    V_2 = Velocidad(V_v(1));
    t_1 = Tiempo(V_v(1)-1);
    t_2 = Tiempo(V_v(1));
    
    v_sim(i) = (t_sim(i) - t_1)/(t_2-t_1)*(V_2-V_1) + V_1;
end

v_sim=v_sim*factAmpVel;
q_t = 0.5 * DensidadAire * v_sim.^2 ; % Carga (KPa)



F_t = q_t * C_d * dc*lc/(NelemC);


if Plot_ft == 1
    figure
    plot(t_sim,F_t)
    ylabel('Fuerza total sobre el elemento [N]')
    xlabel('Tiempo [s]')
end


v_sim = v_sim-Adicional;
v_sim(1:3) = 0;
if Plot_ft == 1
    figure
    plot(t_sim,v_sim)
    ylabel('Velocidad media del viento [m/s]')
    xlabel('Tiempo [s]')
end

end


