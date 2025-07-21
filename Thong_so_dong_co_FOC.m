%% Thong So Dong Co.
VDC = 100;
Rs = 3.2;
Ld = 0.008;
Lq = 0.008;
P = 8; %(poles)
zp = 4;
J = 1.4*10^(-5); 
B = 4.5*10^(-6); %( kg.m^2)
landa_f = 0.0471405 ; %(Wb)
Kt = 3*zp*landa_f/2;
%% Tieu chuan chat luong qua do.
wn_s = 200;
Zeta_s = 1;
wn_i = 2000;
Zeta_i = 1;
%% Thong so PI mien lien tuc.
Kp_i = 2*Zeta_i*wn_i*Lq - Rs;
Ki_i = wn_i^2*Lq;
Kp_s = (2*Zeta_s*wn_s*J - B)/Kt;
Ki_s = (J*wn_s^2)/Kt;
%% Chu ky lay mau he roi rac.
Ts = 1*10^(-3);
Ts_idq = 1*10^(-4);
T_pulse = 5*10^(-5);
%% He so doi tuong roi rac.
tao_i = exp(-(Rs*Ts_idq)/Lq);
tao_s = exp(-(B*Ts)/J); 
bi = (1-tao_i)/Rs;
ai = tao_i;  % Giz = bi/z-ai;
bw = (1-tao_s)*(Kt/B);
aw = tao_s;  % Gwz = bw/z-aw;
%% Thong so PI roi rac.
ri = exp(-(Ts_idq*Zeta_i*wn_i));
Phi_i = Ts_idq*wn_i*sqrt(1-Zeta_i^2); % PTDT: z^2 - 2*cos(Phi_i)*ri*z + ri^2 = 0.
Ki_iz = (-2*cos(Phi_i)*ri +ri^2 + 1)/(bi*Ts_idq);
Kp_iz = -(ri^2 - ai)/bi + Ki_iz*Ts_idq/2;
%---------------------------------------------------------------------------------
Phi_s = Ts*wn_s*sqrt(1-Zeta_s^2);     % PTDT: z^2 - 2*cos(Phi_s)*rs*z + rs^2 = 0.
rs = exp(-(Ts*Zeta_s*wn_s));
Ki_sz = (-2*cos(Phi_s)*rs +rs^2 + 1)/(bw*Ts);
Kp_sz = -(rs^2 - aw)/bw + Ki_sz*Ts/2;
%% Loc dau vao.
tao_f = exp(-(Ki_sz*Ts)/Kp_sz);
A = 1 - tao_f;
C = tao_f;
%% Thong so SMO
Tp_SMO = 5*10^(-5);
Psi = (1 - exp(-(Rs*Tp_SMO)/Lq))/Rs;
Phi = exp(-(Rs*Tp_SMO)/Lq);
wc = 500;
Wc = exp(-500*Tp_SMO);
Omega = exp(-wc*Tp_SMO);
Ksa = 100;
a = 1.2;