%%%% CODIGO CLASE 1/12 - PASO BAJO EQUIVALENTE %%%%

%% 0.- Inicializacion variables
clear, clc
% Constates
c = 2.99793e8;
f_0 = 193e12;
C = 2; % Parametro Chirp

% Parto ahora del dominio frecuencial para llegar al temporal en todo el
% bloque de tratamiento digital de senal

N_frec = 2.^15; % El paso bajo es mas rapido - Exige menos puntos
f_i = linspace(f_0-4e12, f_0+4e12, N_frec+1).'; % 4THz para cada lado desde f0
f_i_PB = f_i-f_0;  % Cambio de variable para formalismo paso bajo eq.
f_max = f_i(length(f_i),1)-f_i(1,1); 
f_muestreo = f_i(2,1)-f_i(1,1);

t_muestreo = 1/f_max;
t_i = t_muestreo.*(0:1:N_frec).';

%% 1.- Definicion senal entrada dominio tiempo

% Parametros de entrada
A_0 = 1;
tau_0 = 20e-12;
tau_d = 2000e-12; % Factor de desplazamiento si fuera necesario

% Con chirp
e_t_0 = A_0.*exp(-(t_i-tau_d).^2./(2*tau_0.^2)).*cos(-C.*((t_i-tau_d).^2)./(2*tau_0.^2));
%e_t_0 = A_0.*exp(-(1+1j*C).*(t_i-tau_d).^2.*inv(2*tau_0.^2));
%e_t_0 = A_0.*exp(-0.5.*((t_i-tau_d)./tau_0).^2);

% En el dominio frecuencial
E_f_0 = fftshift(fft(e_t_0));

%% Representacion temporal
figure, plot(t_i, abs(e_t_0)); 
title('|E_0_t|')
% Ahora tenemos una amplitud compleja que tiene modulo y fase, sin embargo
% a la entrada solo tengo modulo

%% Representacion frecuencial
figure(2), plot(f_i,abs(E_f_0))

%% 2.- Medio dispersivo: Parametros fibra optica

alfa_0 = 0.2; % dB/km
alfa_0_un = 0.1*log(10)*alfa_0/1000; % m-1, necesito unidades naturales

n = 1.452;       % indice de grupo
L = 50e3;        % Longitud en metros
v_grupo = c/n;
t_grupo = L/v_grupo;
beta_0_1 = 0;
beta_0_2 = -20e-24*1e-3; %ps^2/km

% H(L,-w) = H*(L,w) (Conjugado)
% Desprecio beta0 (afecta a todo igual y se quita)
% Desprecio beta0' porque me moveria la ventana. 

H_f_L = exp(-0.5*alfa_0_un*L).*exp(-1i*beta_0_1*L*2*pi*(f_i-f_0)).*exp(-1i*0.5.*beta_0_2*L*(2*pi*(f_i-f_0)).^2);

%% Representacion funcion transferencia
%figure, plot(f_i, abs(H_f_L)) % Modulo de la funcion de transferencia
%figure, plot(f_i, unwrap(angle(H_f_L))) % Modulo de la funcion de transferencia
%set(gca,'Xlim',[f_0-1e12 f_0+1e12])
% La dependecia de la fase de la funcion de transferencia con respecto a la
% fase es cuadratica. Como el tiempo de grupo es la derivada de esta fase
% respecto a w, el tiempo de grupo es lineal respecto a omega.

%% Senal a la salida en el dominio frecuencia
E_f_L = E_f_0 .* H_f_L;

%% Calculo senal salida dominio tiempo
e_t_L = ifft(fftshift(E_f_L));

%% Representacion
figure, plot(t_i,abs(e_t_0)), hold on
plot(t_i,abs(e_t_L),'r')
%set(gca,'Ylim',[0 2.24e-3])

%% Hallamos fase instantanea
fase = unwrap(angle(e_t_L));
f_inst = inv(2*pi).*diff(fase)./diff(t_i);
% Como se distribuyen las componentes espectrales
figure, plot(t_i(1:length(t_i)-1,1),f_inst) 

%% 
subplot(211)
plot(t_i,abs(e_t_0)), hold on
plot(t_i,abs(e_t_L),'r')
legend('|E(0,t)|','|E(L,t)|')
%set(gca,'Xlim',[0.2e-9 0.8e-9])    % Varia segun la representacion
%set(gca,'Ylim',[0 0.05])           % Varia segun la representacion
subplot(212)
plot(t_i(1:length(t_i)-1,1),f_inst)
legend('Frecuencia instantanea')
%set(gca,'Xlim',[0.2e-9 0.8e-9])    % Varia segun la representacion
set(gca,'Ylim',[-2e11 2e11])