%%%%%%%%% REPRESENTAR EL ENSANCHAMIENTO EN FUNCION DE LA LONGITUD %%%%%%%%%
%%%%%%%%% RECORRIDA, DEL VALOR DEL CHIP Y DEL SIGNO DE B_0_2 %%%%%%%%%%%%%%

clear, clc;

% Valores de tau_0 y dispersion
tau_0 = 20e-12;
beta_0_2 = 20e-24*1e-3; %ps^2/km 

% Valores tomados de chirp y de longitud 
Longitud = 0:(1e3):(100e3);
C=[-2 -1 0 1 2];

colores = ['g' 'k' 'r' 'b' 'y' ];

figure(1)
grid on
title('Ensanchamiento del pulso segun L (km) y C (Chirp) para b_0" = + 20 ps^2/km')
xlabel('Longitud: en metros')
ylabel('Ensanchamiento del pulso a la salida: tau_0_salida/tau_0_entrada')

for Chirp = C
    
    relacion_ensanch = []; %vector que recogera el ensanchamiento en funcion de L 
    
    for L=Longitud
         l_salida = sqrt(tau_0.^2 + (beta_0_2.*L).^2./(tau_0.^2)*(1+Chirp.^2)-2*Chirp*beta_0_2.*L);
         relacion_ensanch(end+1) = abs(l_salida/tau_0);
    end

    hold on,plot(Longitud, relacion_ensanch, colores(find(C==Chirp)), 'LineWidth',3)
end

legend('C= -2','C= -1','C= 0','C= 1','C= 2') % Cambiar segun los valores de chirp