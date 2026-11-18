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
% Este programa varre o deficit de geracao (Secao 6.4 e Figura 9 do
% artigo): a carga ativa e reativa de todas as barras e multiplicada por
% c, com a geracao programada mantida, e o PF-PFR e resolvido pelo HELM e
% pelo NR, com limites de P e de Q. Para cada ponto: deficit (MW),
% frequencia pelos dois metodos, unidades saturadas, barras em limite de
% Q, max|V_HELM - V_NR| e tensao minima.


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
disp('.            PF-PFR PELO HELM: VARREDURA DO DEFICIT                .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                       2. DADOS E PARAMETROS
%% ________________________________________________________________________

%% a) Parametros
% tol  : tolerancia do HELM (Wynn) e do NR
% nmax : ordem maxima do HELM
% nit  : iteracoes maximas do NR por passagem
tol  = 1e-10;
nmax = 200;
nit  = 30;


%% b) Sistemas e Fatores de Carga
sistemas = {'118', 'pico', 'leve'};
fatores  = {1:0.01:1.30, 1:0.005:1.08, 1:0.01:1.10};



%% ________________________________________________________________________
%%
%%                            3. VARREDURA
%% ________________________________________________________________________

for is = 1:numel(sistemas)

    %% a) Ponto de Operacao
    [mpc, opc] = FPRP_Caso(sistemas{is});
    Pcarga = sum(mpc.bus(:,3));
    cs  = fatores{is};
    res = NaN(numel(cs), 8);

    fprintf('\n%s\n', sistemas{is});
    fprintf('%7s %8s %9s %9s %5s %5s %8s %6s %3s\n', 'c', 'MW', 'f HELM', 'f NR', 'nsat', 'nq', 'max|dV|', 'Vmin', 'ok');


    %% b) Pontos da Varredura
    for i = 1:numel(cs)

        %% b1) Carga Multiplicada por c
        mc = mpc;
        mc.bus(:,3) = cs(i)*mpc.bus(:,3);
        mc.bus(:,4) = cs(i)*mpc.bus(:,4);
        D = FPRP_Dados_MP(mc, opc);


        %% b2) HELM e NR
        [Vh, lh, Pgh, ih] = FPRP_Resolver(D, 'HELM', tol, nmax);
        [Vn, ln, Pgn, in] = FPRP_Resolver(D, 'NR', tol, nit);


        %% b3) Resultados do Ponto
        % ok: os dois metodos convergiram sem perda da referencia
        ok = ih.conv && in.conv && ~ih.perda && ~in.perda;
        res(i,:) = [(cs(i) - 1)*Pcarga, 60*(1 - lh), 60*(1 - ln), ih.nsat, ih.nq, max(abs(Vh - Vn)), min(abs(Vh)), ok];

        fprintf('%7.3f %8.0f %9.4f %9.4f %5d %5d %8.1e %6.3f %3d\n', cs(i), res(i,1), res(i,2), res(i,3), res(i,4), res(i,5), res(i,6), res(i,7), ok);
    end


    %% c) Gravacao (Figura 9)
    fid = fopen(fullfile('resultados', sprintf('dados_deficit_%s.dat', sistemas{is})), 'w');
    fprintf(fid, 'deficit f fNR nsat nq dV Vmin ok\n');
    fprintf(fid, '%.3f %.6f %.6f %d %d %.3e %.4f %d\n', res.');
    fclose(fid);
end
