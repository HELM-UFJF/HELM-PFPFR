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
% Este programa mede os tempos de execucao do HELM (Wynn) e do NR
% (partida plana) nos casos do estudo do PF-PFR (Tabela 7 e Figura 10 do
% artigo).
%   - Cada execucao inclui o laco completo de limites de P e de Q
%     (FPRP_Resolver) e parte do caso ja convertido: a montagem de D nao e
%     cronometrada;
%   - Uma execucao de aquecimento e descartada e reporta-se a mediana de
%     nrep execucoes;
%   - Rodar no MATLAB, com os demais programas fechados.
% A tabela final e gravada em resultados/tempos.dat.


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
disp('.          PF-PFR PELO HELM: TEMPOS DE EXECUCAO (HELM x NR)        .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                     2. COMPUTADOR E PARAMETROS
%% ________________________________________________________________________

%% a) Especificacoes do Computador
fprintf('\nData           : %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('Arquitetura    : %s\n', computer);
fprintf('MATLAB/Octave  : %s\n', version);

try
    fprintf('Threads        : %d\n', maxNumCompThreads);
catch
end


%% b) Parametros da Medicao
% nrep : execucoes cronometradas por caso
% tol  : tolerancia do HELM (Wynn) e do NR
% nmax : ordem maxima do HELM
% nit  : iteracoes maximas do NR por passagem
nrep = 20;
tol  = 1e-10;
nmax = 200;
nit  = 30;



%% ________________________________________________________________________
%%
%%                              3. CASOS
%% ________________________________________________________________________

%% a) Definicao: sistema, barras das unidades perdidas e descricao
casos = {
    '118',  [],          '118, base';
    '118',  100,         '118, perda da barra 100';
    '118',  [100 89],    '118, perda das barras 100 e 89';
    'pico', [],          '2000 pico, base';
    'pico', 7099,        '2000 pico, perda de 1 nuclear';
    'pico', [7098 7099], '2000 pico, perda de 2 nucleares';
    'leve', [7098 7099], '2000 leve, perda de 2 nucleares';
    };
nc  = size(casos, 1);
res = zeros(nc, 9);


%% b) Laco sobre os Casos
for ic = 1:nc

    %% b1) Ponto de Operacao e Evento
    [mpc, opc] = FPRP_Caso(casos{ic,1});
    mpc.gen(ismember(mpc.gen(:,1), casos{ic,2}), 8) = 0;
    D = FPRP_Dados_MP(mpc, opc);


    %% b2) Aquecimento (descartado)
    [Vh, lh, Pgh, ih] = FPRP_Resolver(D, 'HELM', tol, nmax);
    [Vn, ln, Pgn, in] = FPRP_Resolver(D, 'NR', tol, nit);


    %% b3) Execucoes Cronometradas
    th = zeros(nrep, 1);
    tn = zeros(nrep, 1);

    for ir = 1:nrep
        tic
        [Vh, lh, Pgh, ih] = FPRP_Resolver(D, 'HELM', tol, nmax);
        th(ir) = toc;

        tic
        [Vn, ln, Pgn, in] = FPRP_Resolver(D, 'NR', tol, nit);
        tn(ir) = toc;
    end


    %% b4) Resultados do Caso
    % [tHELM tNR f passos ordens nsat iteracoesNR max|dV| convergiu]
    res(ic,:) = [median(th), median(tn), 60*(1 - lh), ih.passos, ih.ntot, ih.nsat, in.ntot, max(abs(Vh - Vn)), ih.conv && in.conv];

    fprintf('  caso %d de %d concluido (%s)\n', ic, nc, casos{ic,3});
end



%% ________________________________________________________________________
%%
%%                            4. RESULTADOS
%% ________________________________________________________________________

%% a) Tabela na Tela
fprintf('\nTempos (s), mediana de %d execucoes\n\n', nrep);
fprintf('%-34s %9s %9s %6s %9s %6s %6s %5s %6s %9s\n', 'Caso', 'HELM', 'NR', 'H/NR', 'f (Hz)', 'passos', 'ordens', 'nsat', 'it.NR', 'max|dV|');

for ic = 1:nc
    fprintf('%-34s %9.4f %9.4f %6.2f %9.4f %6d %6d %5d %6d %9.1e\n', casos{ic,3}, res(ic,1), res(ic,2), res(ic,1)/res(ic,2), res(ic,3), res(ic,4), res(ic,5), res(ic,6), res(ic,7), res(ic,8));
end

if any(res(:,9) == 0)
    fprintf('\nAtencao: algum caso nao convergiu (ver res(:,9)).\n');
end


%% b) Arquivo de Dados
fid = fopen(fullfile('resultados', 'tempos.dat'), 'w');
fprintf(fid, 'caso tHELM tNR razao f passos ordens nsat itNR\n');

for ic = 1:nc
    fprintf(fid, '%d %.5f %.5f %.3f %.4f %d %d %d %d\n', ic, res(ic,1), res(ic,2), res(ic,1)/res(ic,2), res(ic,3), res(ic,4), res(ic,5), res(ic,6), res(ic,7));
end

fclose(fid);
