%%%% CODIGO TREN PULSOS - PASO BAJO EQUIVALENTE %%%%

%% 0.- Inicializacion variables
clear, clc
% Constates
c = 2.99793e8;
f_0 = 193e12;

% Parto ahora del dominio frecuencial para llegar al temporal en todo el
% bloque de tratamiento digital de senal

N_freq = 2.^15; % El paso bajo es mas rapido - Exige menos puntos
f_i = linspace(f_0-4e12, f_0+4e12, N_freq+1).'; % 4THz para cada lado desde f0
f_i_PB = f_i-f_0;  % Cambio de variable para formalismo paso bajo eq.
f_max = f_i(length(f_i),1)-f_i(1,1); 
f_muestreo = f_i(2,1)-f_i(1,1);

t_muestreo = 1./(f_max);
t_i = t_muestreo.*(0:1:N_freq).';

%% 1.- Definicion senal entrada dominio tiempo

% Parametros de entrada
A_0 = 1;
tau_0 = 20e-12;

% Secuencia de pulsos de entrada:
secuencia = [1 0 1 1 1 0 0 1];
k = length(secuencia); % Numero de pulsos del tren de pulsos
tau_d = (2000e-12 - 200e-12*(k/2-1)) : 200e-12 : (2000e-12+200e-12*(k/2));
% Tener cuidado si meto muchos pulsos necesito ventana mayor para que
% quepan, tendria que poner mas puntos en la estructura de datos

% Crea los diferentes pulsos espaciados
pulsos = zeros(size(t_i,1),k);

for ii=1:k

        pulsos(:,ii) = secuencia(ii).*A_0.*exp(-0.5.*((t_i-tau_d(ii))./tau_0).^2);
        
end

% Conforma un unico pulso con todos los pulsos definidos para reaprovechar
% el codigo
e_t_0 = zeros(length(pulsos),1);
for ii = 2:2:k
    
    e_t_0 = e_t_0 + pulsos(:,ii-1) + pulsos(:,ii);
    
end

E_f_0 = fftshift(fft(e_t_0));

%% Representacion temporal
figure, plot(t_i, abs(e_t_0)); 
title('|E_0_t|')
set(gca,'Xlim',[1e-9 3e-9])
% Ahora tenemos una amplitud compleja que tiene modulo y fase, sin embargo
% a la entrada solo tengo modulo, la fase es cero

%% Representacion frecuencial
figure(2), plot(f_i,abs(E_f_0))

%% 2.- Medio dispersivo: Parametros fibra optica

alfa_0 = 0.2; % dB/km
alfa_0_un = 0.1*log(10)*alfa_0/1000; % m^-1, queremos unidades naturales

n = 1.452;                  % Indice de grupo
L = 50e3;                   % Longitud: En metros
v_grupo = c/n;
t_grupo = L/v_grupo;
beta_0_1 = 0;
beta_0_2 = 20e-24*1e-3;    %ps^2/km

% H(L,-w) = H*(L,w) (Conjugado)
% Desprecio beta0 (afecta a todo igual y se quita)
% Desprecio beta0' porque me moveria la ventana. Si lo necesito despues
% puedo tener en cuenta el tg (minuto 06:30) - Por si acaso

H_f_L = exp(-0.5*alfa_0_un*L).*exp(-1i*beta_0_1*L*2*pi*(f_i-f_0)).*exp(-1i*0.5.*beta_0_2*L*(2*pi*(f_i-f_0)).^2);

%% Representacion funcion transferencia
%figure, plot(f_i, abs(H_f_L)) % Modulo de la funcion de transferencia
%figure, plot(f_i, unwrap(angle(H_f_L))) % Modulo de la funcion de transferencia
%set(gca,'Xlim',[f_0-1e12 f_0+1e12])
% La dependecia de la fase de la funcion de transferencia con respecto a la
% fase es cuadratica. Como el tiempo de grupo es la derivada de esta fase
% respecto a w, el tiempo de grupo es lineal respecto a omega.

%% Senal a la salida en el dominio frecuencia

%E_f_L = zeros(size(E_f_0));
%e_t_L = zeros(size(E_f_0));

%for ii=1:k

%      E_f_L(:,ii) = E_f_0(:,ii).*H_f_L;    
%      e_t_L(:,ii) = ifft(fftshift(E_f_L(:,ii)));
      
%end

E_f_L = E_f_0.*H_f_L;
e_t_L = ifft(fftshift(E_f_L));

%% Representacion
figure, plot(t_i,abs(e_t_0)), hold on
plot(t_i,abs(e_t_L),'r')
%set(gca,'Ylim',[0 2.24e-3])

%%
fase = unwrap(angle(e_t_L));
f_inst = inv(2*pi).*diff(fase)./diff(t_i);
% Como se distribuyen las componentes espectrales
figure, plot(t_i(1:length(t_i)-1,1),f_inst) 

%% 
subplot(211)
plot(t_i,abs(e_t_0)), hold on
plot(t_i,abs(e_t_L),'r')
legend('|E(0,t)|','|E(L,t)|')
set(gca,'Xlim',[0.9e-9 3.3e-9])     % Cambia segun la representacion

subplot(212)
plot(t_i(1:length(t_i)-1,1),f_inst)
legend('Frecuencia instantanea')
set(gca,'Xlim',[0.9e-9 3.3e-9])     % Cambia segun la representacion
set(gca,'Ylim',[-1e11 1e11])        