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
% Este programa reproduz o caso de Wang et al. (Energy Proceedings, 2020):
% IEEE 39 barras, saida do gerador da barra 30, regulacao primaria sem
% limites (Secao 6.1 do artigo).
%   - O ganho de cada gerador (G31 a G39) vem da Tabela 1 de Wang et al.:
%     K_G = (P1 - P0)/|df| (p.u./Hz) e 1/R = 60*K_G;
%   - A sensibilidade da carga, K_L = kap*Pd, e calibrada pelo metodo da
%     secante para reproduzir df = -0.1143 Hz;
%   - HELM e NR sao comparados com a Tabela 1 de Wang et al.
% Como os ganhos e a sensibilidade vem da mesma fonte, a comparacao e uma
% verificacao de consistencia, e nao uma validacao independente.


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
disp('.        PF-PFR PELO HELM: IEEE 39 BARRAS (WANG ET AL.)            .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                    2. DADOS DO ARTIGO E DO CASO
%% ________________________________________________________________________

%% a) Tabela 1 de Wang et al. (p.u., base 100 MVA)
% Geracao inicial e final dos geradores G31 a G39 e desvio de frequencia
P0   = [6.7787 6.5000 6.3200 5.0800 6.5000 5.6000 5.4000 8.3000 10.0000];
P1   = [7.1360 6.7072 6.5685 5.2736 6.6963 5.8210 5.5612 8.5472 10.4192];
df_w = -0.1143;


%% b) IEEE 39 com o Gerador da Barra 30 Fora, sem Limites de P
[D, mpc] = Wang39_Dados(true);
D.Pmax(:) =  Inf;
D.Pmin(:) = -Inf;


%% c) Ganhos dos Geradores com Estatismo (barras 31 a 39)
[~, idx] = ismember(D.ger, 31:39);
KG    = (P1 - P0)/abs(df_w);
D.Reg = 1./(60*KG(idx).');



%% ________________________________________________________________________
%%
%%                  3. CALIBRACAO DA SENSIBILIDADE DA CARGA
%% ________________________________________________________________________

%% a) Pontos Iniciais da Secante
% Desvio de frequencia (Hz) pelo NR para kap = 1 e kap = 3
k0 = 1;
k1 = 3;
D.KL = k0*D.Pd;
[~, l0] = FPRP_Resolver(D, 'NR', 1e-12, 30);
D.KL = k1*D.Pd;
[~, l1] = FPRP_Resolver(D, 'NR', 1e-12, 30);
f0 = -60*l0 - df_w;
f1 = -60*l1 - df_w;


%% b) Metodo da Secante: df(kap) = df_w
for it = 1:20
    k2 = k1 - f1*(k1 - k0)/(f1 - f0);
    D.KL = k2*D.Pd;
    [~, l2] = FPRP_Resolver(D, 'NR', 1e-12, 30);
    f2 = -60*l2 - df_w;
    k0 = k1;
    f0 = f1;
    k1 = k2;
    f1 = f2;

    if abs(f2) < 1e-12
        break
    end
end
kap = k1;

fprintf('\nkap calibrado = %.5f (K_L total = %.4f p.u./Hz)\n', kap, kap*sum(D.Pd)/60);



%% ________________________________________________________________________
%%
%%                   4. HELM, NR E TABELA DE WANG ET AL.
%% ________________________________________________________________________

%% a) Solucoes com o kap Calibrado
D.KL = kap*D.Pd;
[Vh, lh, Pgh, ih] = FPRP_Resolver(D, 'HELM', 1e-10, 150);
[Vn, ln, Pgn]     = FPRP_Resolver(D, 'NR', 1e-12, 30);

fprintf('HELM: df = %.6f Hz (ordem %d) | NR: df = %.6f Hz | max|V_HELM - V_NR| = %.1e p.u.\n', -60*lh, ih.n, -60*ln, max(abs(Vh - Vn)));


%% b) Geracao de Cada Gerador
fprintf('\n  G     P0 (p.u.)  P1 Wang   P1 HELM   P1 NR    HELM - Wang\n');

for i = 1:9
    k = find(D.ger == 30 + i);
    fprintf('  G%d  %8.4f  %8.4f  %8.4f  %8.4f   %9.1e\n', 30 + i, P0(i), P1(i), Pgh(k), Pgn(k), Pgh(k) - P1(i));
end

fprintf('\nmax |HELM - Wang| = %.2e p.u.\n', max(abs(Pgh(idx > 0) - P1(idx(idx > 0)).')));


%% c) Arquivo de Dados
% Colunas: gerador, P0, P1 de Wang et al., P1 HELM e P1 NR (p.u.)
fid = fopen(fullfile('resultados', 'dados_wang39.dat'), 'w');
fprintf(fid, 'G P0 P1Wang P1HELM P1NR\n');

for i = 1:9
    k = find(D.ger == 30 + i);
    fprintf(fid, '%d %.6f %.6f %.6f %.6f\n', 30 + i, P0(i), P1(i), Pgh(k), Pgn(k));
end

fclose(fid);
