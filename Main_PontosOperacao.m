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
% Este programa resume os pontos de operacao do estudo (Tabela 4 e
% Figura 6 do artigo): IEEE 118 e ACTIVSg2000 no pico e na hora leve,
% despachados por OPF AC com requisito de reserva (geracao/Gera_Base_FPRP).
% Para cada ponto: carga, unidades candidatas, unidades em Pmax, unidades
% com estatismo, rigidez das unidades e da carga (MW/Hz) e a folga de cada
% candidata, (Pmax - Pg)/Pmax, em ordem decrescente.


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
disp('.              PF-PFR PELO HELM: PONTOS DE OPERACAO                .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                       2. RESUMO DOS PONTOS
%% ________________________________________________________________________

%% a) Sistemas
sistemas = {'118', 'pico', 'leve'};
f0 = 60;

fprintf('\n%-20s %8s %6s %6s %6s %10s %10s\n', 'Ponto', 'Carga', 'Cand.', 'Pmax', 'Regul.', 'Unid.MW/Hz', 'Carga MW/Hz');


%% b) Laco sobre os Pontos de Operacao
for is = 1:numel(sistemas)

    %% b1) Caso, Candidatas e Unidades com Estatismo
    [mpc, opc, nome] = FPRP_Caso(sistemas{is});
    D = FPRP_Dados_MP(mpc, opc);

    switch sistemas{is}
        case '118'
            cand = mpc.gen(:,9) > 0 & mpc.gen(:,8) == 1;
        otherwise
            cand = ismember(mpc.genfuel(:), {'ng', 'coal'}) & mpc.gen(:,8) == 1 & mpc.gen(:,9) > mpc.gen(:,10);
    end
    g = find(cand);


    %% b2) Indicadores
    carga = sum(mpc.bus(:,3));
    emMax = sum(mpc.gen(g,2) >= mpc.gen(g,9) - 1e-3);
    kU    = sum(1./D.Reg)*mpc.baseMVA/f0;
    kL    = sum(D.KL)*mpc.baseMVA/f0;

    fprintf('%-20s %8.0f %6d %6d %6d %10.0f %10.0f\n', nome, carga, numel(g), emMax, numel(D.ger), kU, kL);


    %% b3) Folga de Cada Candidata (%), em Ordem Decrescente
    folga = 100*(mpc.gen(g,9) - mpc.gen(g,2))./mpc.gen(g,9);
    folga = sort(max(folga, 0), 'descend');
    x     = 100*((1:numel(folga)).' - 0.5)/numel(folga);

    fid = fopen(fullfile('resultados', sprintf('dados_folga_%s.dat', sistemas{is})), 'w');
    fprintf(fid, 'x folga\n');
    fprintf(fid, '%.3f %.4f\n', [x folga].');
    fclose(fid);
end
