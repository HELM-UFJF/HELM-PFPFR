%                                                                       JMJ
%
%       "Aquilo que Deus quer para nos e muito mais belo do que tudo
%              que poderiamos pedir com a nossa imaginacao."
%
%                                 Chiara Corbella Petrillo


%% ________________________________________________________________________
%%
%%                        1. INFORMACOES INICIAIS
%% ________________________________________________________________________

%% a) Comentarios Iniciais
% Este programa executa o tutorial de 3 barras do artigo (Secao 3.5):
%   - coeficientes das series de V3, Q2 e lambda e convergencia da soma
%     parcial e da estimativa de Wynn (Tabela 3 e Figura 4);
%   - caso com limite de Q na barra 2 (Q2max = 0.60 p.u.).
% Dados: case3_tutorial_c_perdas. As duas unidades tem R = 5 % na base de
% 100 MVA; o despacho da referencia, 233.49 MW, da 60 Hz na carga nominal,
% e todas as cargas sao aumentadas em 20 %.


%% b) Limpeza de Figuras e Variaveis
clc
clear all
close all
format short
warning off


%% c) Pastas das Funcoes e dos Dados
addpath('funcoes', 'dados');

if ~exist('resultados', 'dir')
    mkdir('resultados');
end


%% d) Cabecalho
disp('....................................................................')
disp('.                                                                  .')
disp('.          PF-PFR PELO HELM: TUTORIAL DE 3 BARRAS                  .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                       2. DADOS E CENARIO
%% ________________________________________________________________________

%% a) Sistema e Despacho da Referencia
mpc = case3_tutorial_c_perdas();
D   = FPRP_Dados_MPC(mpc);
D.Pg(1) = 2.33489;


%% b) Carga Aumentada em 20 %
D.Pd = 1.2*D.Pd;
D.Qd = 1.2*D.Qd;



%% ________________________________________________________________________
%%
%%                3. COEFICIENTES E CONVERGENCIA DAS SERIES
%% ________________________________________________________________________

%% a) Series do HELM (sem limites)
% Rigidez das unidades e injecoes em lambda = 0
beta = zeros(D.N, 1);
beta(D.ger) = 1./D.Reg;
Pnet = D.Pg - D.Pd;
Qnet = -D.Qd;

% Tolerancia muito pequena para obter os coeficientes ate a ordem 16
[V1, lam1, hi] = FPRP_HELM_Nucleo(D.Ysr, D.Ysh, D.tipo, D.Vesp, Pnet, Qnet, beta, D.ref, 1e-16, 20, 30);


%% b) Referencia pelo NR
[Vn, lamNR] = FPRP_Resolver(D, 'NR', 1e-14, 30);

fprintf('\nlambda_NR = %.7f (f = %.4f Hz), |V3| = %.4f p.u.\n', lamNR, 60*(1 - lamNR), abs(Vn(3)));


%% c) Tabela de Coeficientes e Erros
fprintf('\n n   Re{V3[n]}     Im{V3[n]}     Q2[n]         lambda[n]     |soma-NR|  |Wynn-NR|\n');

% Estimativas de Wynn de lambda, uma antidiagonal nova por ordem
wl    = [];
nMax  = min(16, size(hi.lam, 2) - 1);
dados = zeros(nMax+1, 4);

for n = 0:nMax
    soma = sum(hi.lam(1:n+1));
    [wl, est] = Wynn_Passo(wl, soma);
    es = abs(soma - lamNR);
    ew = abs(real(est) - lamNR);
    dados(n+1,:) = [n, abs(hi.lam(n+1)), es, ew];

    if n <= 8
        fprintf('%2d %13.4e %13.4e %13.4e %13.4e %10.2e %10.2e\n', n, real(hi.V(3,n+1)), imag(hi.V(3,n+1)), hi.Qh(2,n+1), hi.lam(n+1), es, ew);
    end
end


%% d) Dados da Figura de Convergencia
% Valores nulos sao trocados por 1e-17 (escala logaritmica)
dados(:,2:4) = max(dados(:,2:4), 1e-17);
fid = fopen(fullfile('resultados', 'dados_convergencia.dat'), 'w');
fprintf(fid, 'n coef soma wynn\n');
fprintf(fid, '%d %.6e %.6e %.6e\n', dados.');
fclose(fid);



%% ________________________________________________________________________
%%
%%                     4. LIMITE DE Q NA BARRA 2
%% ________________________________________________________________________

%% a) Sem Limite
[Vh, lh, ~, ih] = FPRP_Resolver(D, 'HELM', 1e-12, 100);
fprintf('\nSem limite de Q: Q2 = %.3f p.u., f = %.4f Hz\n', ih.Qbus(2), 60*(1 - lh));


%% b) Com Q2max = 0.60 p.u.
D.Qmax(2) = 0.60;
[Vh, lh] = FPRP_Resolver(D, 'HELM', 1e-12, 100);
[Vn, ln] = FPRP_Resolver(D, 'NR', 1e-13, 30);

fprintf('Com Q2max = 0.60 p.u.: |V2| = %.4f p.u., f = %.4f Hz, max|V_HELM - V_NR| = %.1e p.u.\n', abs(Vh(2)), 60*(1 - lh), max(abs(Vh - Vn)));
