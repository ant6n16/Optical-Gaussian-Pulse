%%%% CODIGO TREN PULSOS - ESPECTRO COMPLETO %%%%

%% 0.- Inicializacion variables
clear, clc
% Constates
c = 2.99793e8;
f_0 = 193e12;

% Variables temporales y frecuenciales
t_muestreo = 1/f_0/11;  % 11 puntos por cada periodo de la portadora
f_max = 1/(2.*t_muestreo);

N_freq = 2^24;
f_muestreo = f_max/N_freq;
f_i = f_muestreo.*(0:1:N_freq).';
t_i = t_muestreo.*(0:1:N_freq).';

%% 1.- Definicion senal entrada dominio tiempo

% Parametros de entrada
A_0 = 1;
tau_0 = 20e-12;

% Secuencia de pulsos de entrada
secuencia = [1 0 1 1 1 0 0 1];
k = length(secuencia); % Numero de pulsos del tren de pulsos
tau_d = (4000e-12-200e-12*(k/2-1)):200e-12:(4000e-12+200e-12*(k/2));
% Tener cuidado si meto muchos pulsos necesito ventana mayor para que
% quepan, tendria que poner mas puntos en la estructura de datos

% Crea los diferentes pulsos espaciados
pulsos = zeros(size(t_i,1),k);

for ii=1:k

        pulsos(:,ii) = secuencia(ii).*A_0.*exp(-(t_i-tau_d(ii)).^2.*inv(2*tau_0.^2)).*cos(2.*pi.*f_0.*(t_i-tau_d(ii)));
        
end

% Conforma un unico pulso con todos los pulsos definidos para reaprovechar
% el codigo
e_t_0 = zeros(length(pulsos),1);
for ii = 2:2:k
    
    e_t_0 = e_t_0 + pulsos(:,ii-1) + pulsos(:,ii);
    
end

E_f_0 = fft(e_t_0,2*N_freq); % Dos veces N_freq para meter las frec. negativas
E_f_0 = E_f_0(1:N_freq+1,1); % Deshago para no tener el doble de puntos

%% Representacion temporal
figure(2), plot(t_i, e_t_0);
set(gca,'Xlim',[3e-9 5e-9]) % Cambia segun la representacion

%% Representacion frecuencial
figure, plot(f_i,abs(E_f_0))

%% 2.- Medio dispersivo: Parametros fibra optica

alfa_0 = 0.2; % dB/km
alfa_0_un = 0.1*log(10)*alfa_0/1000; % m^(-1) - necesito unidades naturales

n = 1.452;             % Indice de grupo
L = 50e3;              % Longitud: en metros
v_grupo = c/n;
t_grupo = L/v_grupo;
beta_0_1 = 0;
beta_0_2 = -20e-24*1e-3; %ps^2/km

% H(L,-w) = H*(L,w) (Conjugado)
% Desprecio beta0 (afecta a todo igual y se quita)
% Desprecio beta0' porque me moveria la ventana. Si lo necesito despues
% puedo tener en cuenta el tg (minuto 06:30) - Por si acaso

H_f_L = exp(-0.5*alfa_0_un*L).*exp(-1i*beta_0_1*L*2*pi*(f_i-f_0)).*exp(-1i*0.5.*beta_0_2*L*(2*pi*(f_i-f_0)).^2);

%% Representacion funcion transferencia
figure, plot(f_i, abs(H_f_L)) % Modulo de la funcion de transferencia
figure, plot(f_i, unwrap(angle(H_f_L))) % Modulo de la funcion de transferencia
set(gca,'Xlim',[f_0-1e12 f_0+1e12])
% La dependecia de la fase de la funcion de transferencia con respecto a la
% fase es cuadratica. Como el tiempo de grupo es la derivada de esta fase
% respecto a w, el tiempo de grupo es lineal respecto a omega.

%% Senal a la salida en el dominio frecuencia
E_f_L = E_f_0 .* H_f_L;

%% Calculo senal salida dominio tiempo
e_t_L = real(ifft([ E_f_L; conj(E_f_L(N_freq:-1:2,1)) ]));
e_t_L = e_t_L(1:N_freq+1,1);

%% Representacion
figure, plot(t_i,e_t_0,'b'), hold on
plot(t_i,e_t_L,'r')
set(gca,'Xlim',[3e-9 5e-9])    % Cambia segun la representacion
%set(gca,'Ylim',[-0.2 0.2])    % Cambia segun la representacion