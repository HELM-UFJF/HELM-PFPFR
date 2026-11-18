function mpc = case3_tutorial_c_perdas

%% MATPOWER Case Format : Version 2
mpc.version = '2';


%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;


%% bus data
% OBS: 1) Bs       ==> Nao precisara dividir por 2 dentro do codigo
% OBS: 2) Reg      ==> Estatismo, em p.u., na base da maquina
% OBS: 3) base_maq ==> Potencia base da maquina, em MVA
mpc.bus = [
%  bus_i  type	    Pd      Qd     Gs Bs/2 area	     Vm    Va baseKV   zone
%  Vmax    Vmin    PG     QG     Qmim      Qmax     Reg    base_maq
	1       3        0	    0      0	0	1       1.00    0	230     1       1.05    0.95     230     0     -20000     15000  0.05  100;
	2       2       80      0      0	0	1       1.00	0	230     1       1.05	0.95      50     0     -10000     10000  0.05  100;
	3       1      200    124      0    0   1       1.00  	0  	230     1       1.05    0.95       0     0     -99999     99999  0       0;
  ];

%% Adicao de 20% em cada barra de carga
% mpc.bus(:,3) = mpc.bus(:,3)*1.2;
%% Reducao de 20% em cada barra de carga
% mpc.bus(:,3) = mpc.bus(:,3)*0.8;

%% branch data
% OBS: 1) SHL (b) ==> deve-se /2 dentro do codigo e ja esta em [pu]
%      2) r, x    ==> [pu]
%      3) TAP     ==> onde nao tem transformador, o TAP deve ser 1
mpc.branch = [
%  fbus	tbus r      x       b	       rateA	rateB	rateC  TAP	  angle	status	  angmin	angmax
    1	2	0.01008	0.05040	0.05125*2	130     130     130     1       0       1      -360     360;
	1	3	0.00744	0.03720	0.03875*2	130     130     130     1       0       1      -360     360;
	2	3	0.01272	0.06360	0.06375*2	130     130     130     1       0       1      -360     360;
];