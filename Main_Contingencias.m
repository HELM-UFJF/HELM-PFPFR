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
% Este programa simula as perdas de geracao do estudo (Secoes 6.2 e 6.3 do
% artigo) no IEEE 118 e no ACTIVSg2000 no pico e na hora leve:
%   - fluxo de potencia classico (referencia unica, MATPOWER) contra o
%     PF-PFR (Tabela 5 e Figura 7);
%   - PF-PFR com limites de P e Q pelo HELM e pelo NR, conferido pelo
%     MATPOWER, e o mesmo caso sem limites de Q (Tabela 6);
%   - perfil de tensoes do IEEE 118 na perda dupla (Figura 8).
% Requer o MATPOWER no path (runpf).


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

if exist('runpf', 'file') ~= 2
    error('Main_Contingencias: o MATPOWER (runpf) nao esta no path.');
end


%% d) Cabecalho
disp('....................................................................')
disp('.                                                                  .')
disp('.            PF-PFR PELO HELM: PERDAS DE GERACAO                   .')
disp('.                                                                  .')
disp('....................................................................')



%% ________________________________________________________________________
%%
%%                       2. DADOS E PARAMETROS
%% ________________________________________________________________________

%% a) Parametros
tol  = 1e-10;
nmax = 200;
nit  = 30;
kap  = 1;
opMP = mpoption('verbose', 0, 'out.all', 0, 'exp.use_legacy_core', 1);


%% b) Casos: sistema, barras das unidades perdidas, descricao
casos = {
    '118',  100,         '118, perda da barra 100';
    '118',  89,          '118, perda da barra 89';
    '118',  [100 89],    '118, perda das duas';
    'pico', 7099,        '2000 pico, 1 nuclear';
    'pico', [7098 7099], '2000 pico, 2 nucleares';
    'leve', 7099,        '2000 leve, 1 nuclear';
    'leve', [7098 7099], '2000 leve, 2 nucleares';
    };
nc = size(casos, 1);



%% ________________________________________________________________________
%%
%%                3. FLUXO CLASSICO E PF-PFR (TABELA 5)
%% ________________________________________________________________________

fprintf('\n%-26s %7s | %7s %7s %4s %6s | %8s %7s %7s %4s\n', 'Caso', 'Prmax', 'Pr', 'Smax', 'n>1', 'c/ Q', 'f (Hz)', 'Pr', 'Smax', 'n>1');

for ic = 1:nc

    %% a) Caso Base e Evento
    [mpc, opc] = FPRP_Caso(casos{ic,1});
    D0 = FPRP_Dados_MP(mpc, opc);
    Vb = FPRP_Resolver(D0, 'HELM', tol, nmax);
    Sb = FPRP_Carregamento(mpc, D0, Vb);
    mpc.gen(ismember(mpc.gen(:,1), casos{ic,2}), 8) = 0;
    D  = FPRP_Dados_MP(mpc, opc);


    %% b) Fluxo Classico (a referencia assume todo o deficit)
    mc = mpc;
    tp = D.tipo;
    tp(D.tipo == 0) = 3;
    mc.bus(:,2) = tp;
    gr = find(mc.gen(:,1) == D.num(D.ref) & mc.gen(:,8) == 1);

    r0 = runpf(mc, opMP);
    Sc = FPRP_Carregamento(mpc, D, r0.bus(:,8).*exp(1j*pi/180*r0.bus(:,9)));
    r1 = runpf(mc, mpoption(opMP, 'pf.enforce_q_lims', 1));


    %% c) PF-PFR
    [Vh, lh, Pgh] = FPRP_Resolver(D, 'HELM', tol, nmax);
    Sp = FPRP_Carregamento(mpc, D, Vh);
    ir = D.ger == D.ref;


    %% d) Linha da Tabela 5
    convQ = {'diverge', 'converge'};
    fprintf('%-26s %7.0f | %7.0f %7.1f %4d %6s | %8.3f %7.0f %7.1f %4d\n', casos{ic,3}, sum(mc.gen(gr,9)), ...
        sum(r0.gen(gr,2)), max(Sc), sum(Sc > 100), convQ{r1.success + 1}, ...
        60*(1 - lh), sum(Pgh(ir))*mpc.baseMVA, max(Sp), sum(Sp > 100));


    %% e) Dados da Figura 7 (20 ramos mais carregados)
    if strcmp(casos{ic,3}, '2000 pico, 2 nucleares') || strcmp(casos{ic,3}, '118, perda da barra 100')
        n  = 20;
        sc = sort(Sc, 'descend');
        sp = sort(Sp, 'descend');
        sb = sort(Sb, 'descend');
        nomeArq = 'dados_carreg_118.dat';

        if strcmp(casos{ic,1}, 'pico')
            nomeArq = 'dados_carreg_pico.dat';
        end

        fid = fopen(fullfile('resultados', nomeArq), 'w');
        fprintf(fid, 'k classico pfpr base\n');
        fprintf(fid, '%d %.3f %.3f %.3f\n', [(1:n).' sc(1:n) sp(1:n) sb(1:n)].');
        fclose(fid);
    end
end



%% ________________________________________________________________________
%%
%%               4. PF-PFR COM E SEM LIMITES DE Q (TABELA 6)
%% ________________________________________________________________________

fprintf('\n%-26s %6s %8s %8s %4s %4s %6s %8s %8s %8s\n', 'Caso', 'MW', 'f (Hz)', 'f s/Q', 'nsat', 'nq', 'Vmin', 'dV_NR', 'dV_MP', 'dPref');

for ic = 1:nc

    %% a) Caso e Evento
    [mpc, opc] = FPRP_Caso(casos{ic,1});
    perdida = ismember(mpc.gen(:,1), casos{ic,2}) & mpc.gen(:,8) == 1;
    Pperda  = sum(mpc.gen(perdida,2));
    mpc.gen(perdida, 8) = 0;
    D = FPRP_Dados_MP(mpc, opc);


    %% b) HELM e NR com Limites de P e Q
    [Vh, lh, Pgh, ih] = FPRP_Resolver(D, 'HELM', tol, nmax);
    [Vn, ln]          = FPRP_Resolver(D, 'NR', tol, nit);


    %% c) Conferencia pelo MATPOWER
    c = FPRP_Confere_MP(mpc, D, Vh, lh, Pgh, kap, ih.qstate);


    %% d) Mesmo Caso sem Limites de Q
    Ds = D;
    Ds.Qmin(:) = -Inf;
    Ds.Qmax(:) =  Inf;
    [Vs, ls] = FPRP_Resolver(Ds, 'HELM', tol, nmax);


    %% e) Linha da Tabela 6
    fprintf('%-26s %6.0f %8.4f %8.4f %4d %4d %6.3f %8.1e %8.1e %8.1e\n', casos{ic,3}, Pperda, 60*(1 - lh), 60*(1 - ls), ...
        ih.nsat, ih.nq, min(abs(Vh)), max(abs(Vh - Vn)), c.dV, abs(c.dPref));


    %% f) Dados da Figura 8 (IEEE 118, perda das duas unidades)
    if strcmp(casos{ic,3}, '118, perda das duas')
        lim = ih.qstate ~= 0;

        fid = fopen(fullfile('resultados', 'dados_tensao_118.dat'), 'w');
        fprintf(fid, 'barra comQ semQ lim\n');
        fprintf(fid, '%d %.5f %.5f %d\n', [D.num abs(Vh) abs(Vs) lim].');
        fclose(fid);

        fid = fopen(fullfile('resultados', 'dados_tensao_118_lim.dat'), 'w');
        fprintf(fid, 'barra comQ\n');
        fprintf(fid, '%d %.5f\n', [D.num(lim) abs(Vh(lim))].');
        fclose(fid);
    end
end
