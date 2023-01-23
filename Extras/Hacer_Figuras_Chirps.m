e_t_c0 = e_t_0; % 0 2, -20
e_L_c0 = e_t_L;

e_t_c1 = e_t_0; % 2, 20
e_L_c1 = e_t_L;

e_t_c2 = e_t_0; % -2, -20
e_L_c2 = e_t_L;

e_t_c3 = e_t_0; % -2, 20
e_L_c3 = e_t_L;

e_t_c4 = e_t_0; % -4
e_L_c4 = e_t_L;
%%

subplot(421)
plot(t_i,e_t_c0)
set(gca,'Xlim',[3.5e-9 4.5e-9])
title('Chirp = 2  , beta_0" = -20')
subplot(422)
plot(t_i,e_L_c0)
set(gca,'Xlim',[3.5e-9 4.5e-9])


subplot(423)
plot(t_i,e_t_c1)
set(gca,'Xlim',[3.5e-9 4.5e-9])
title('Chirp = 2  , beta_0" = 20')
subplot(424)
plot(t_i,e_L_c1)
set(gca,'Xlim',[3.5e-9 4.5e-9])

subplot(425)
plot(t_i,e_t_c2)
title('Chirp = -2  , beta_0" = -20')
set(gca,'Xlim',[3.5e-9 4.5e-9])
subplot(426)
plot(t_i,e_L_c2)
set(gca,'Xlim',[3.5e-9 4.5e-9])

subplot(427)
plot(t_i,e_t_c3)
title('Chirp = -2  , beta_0" = 20')
set(gca,'Xlim',[3.5e-9 4.5e-9])
subplot(428)
plot(t_i,e_L_c3)
set(gca,'Xlim',[3.5e-9 4.5e-9])
%%
subplot(529)
plot(t_i,e_t_c4)
title('Chirp = -4')
set(gca,'Xlim',[3.5e-9 4.5e-9])

subplot(5,2,10)
plot(t_i,e_L_c4)
set(gca,'Xlim',[3.5e-9 4.5e-9])




