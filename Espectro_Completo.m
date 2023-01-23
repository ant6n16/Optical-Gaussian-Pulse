%%%% CODIGO CLASE 1/12 - ESPECTRO COMPLETO %%%%

%% 0.- Inicializacion variables
clear, clc
% Constates
c = 2.99793e8;
f_0 = 193e12;
C = 2; % Parametro Chirp

% Variables temporales y frecuenciales
t_muestreo = 1/f_0/11;  % 11 puntos por cada periodo de la portadora
f_max = 1/(2.*t_muestreo);

% Numero de puntos tomados en la ventana de frecuencias
N_freq = 2^24;
f_muestreo = f_max/N_freq;
f_i = f_muestreo.*(0:1:N_freq).';
t_i = t_muestreo.*(0:1:N_freq).';

%% 1.- Definicion senal entrada dominio tiempo

% Parametros de entrada
A_0 = 1;
tau_0 = 20e-12;
tau_d = 4000e-12; % Factor de desplazamiento si fuera necesario

% Con chirp
e_t_0 = A_0.*exp(-(t_i-tau_d).^2./(2*tau_0.^2)).*cos(2.*pi.*f_0.*(t_i-tau_d)-C.*((t_i-tau_d).^2)./(2*tau_0.^2));
%e_t_0 = A_0.*exp(-(t_i-tau_d).^2.*inv(2*tau_0.^2)).*cos(2.*pi.*f_0.*(t_i-tau_d));

% En el dominio frecuencial
E_f_0 = fft(e_t_0,2*N_freq); % Dos veces N_freq para meter las frec. negativas
E_f_0 = E_f_0(1:N_freq+1,1); % Deshago para no tener el doble de puntos

%% Representacion temporal
figure(2), plot(t_i, e_t_0);
set(gca,'Xlim',[3.5e-9 4.5e-9]) % Variar en funcion de la representacion

%% Representacion frecuencial
figure, plot(f_i,abs(E_f_0))

%% 2.- Medio dispersivo: Parametros fibra optica

alfa_0 = 0.2; % dB/km
alfa_0_un = 0.1*log(10)*alfa_0/1000; % m-1, necesito unidades naturales

n = 1.452;  % indice de grupo
L = 50e3;   % Longitud, en metros
v_grupo = c/n;
t_grupo = L/v_grupo;
beta_0_1 = 0;
beta_0_2 = -20e-24*1e-3; %ps^2/km

% H(L,-w) = H*(L,w) (Conjugado)
% Desprecio beta0 (afecta a todo igual y se quita)
% Desprecio beta0' porque me moveria la ventan de representacion.

H_f_L = exp(-0.5*alfa_0_un*L).*exp(-1i*beta_0_1*L*2*pi*(f_i-f_0)).*exp(-1i*0.5.*beta_0_2*L*(2*pi*(f_i-f_0)).^2);

%% Representacion funcion transferencia
figure, plot(f_i, abs(H_f_L)) % Modulo de la funcion de transferencia
figure, plot(f_i, unwrap(angle(H_f_L))) % Fase de la funcion de transferencia
set(gca,'Xlim',[f_0-1e12 f_0+1e12])
% La dependecia de la fase de la funcion de transferencia con respecto a la
% fase es cuadratica. Como el tiempo de grupo es la derivada de esta fase
% respecto a w, el tiempo de grupo es lineal respecto a omega.

%% Senal a la salida en el dominio frecuencia
E_f_L = E_f_0 .* H_f_L;

%% Calculo senal salida dominio tiempo
e_t_L = real(ifft([ E_f_L; conj(E_f_L(N_freq:-1:2,1)) ] ));
% Tengo que desdoblar mi ventana de la frecuencia para tener las
% frecuencias negativas, por eso anado la parte conjugada

e_t_L = e_t_L(1:N_freq+1,1); % Me quedo con la ventana que me interesa

%% Representacion
figure, plot(t_i,e_t_0), hold on
plot(t_i,e_t_L,'r')
set(gca,'Xlim',[3.5e-9 4.5e-9]) % Varia en funcion de la representacion
set(gca,'Ylim',[-0.2 0.2])      % Varia en funcion de la representacion