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
% Este programa gera os pontos de operacao do estudo do PF-PFR (fluxo de
% potencia com regulacao primaria) nos sistemas IEEE 118 e ACTIVSg2000.
% Requer o MATPOWER (runopf e toggle_reserves).
%
% Despacho com reserva de regulacao primaria, como num mercado:
%   1) Requisito de reserva do sistema (MW), atendido pelas unidades com
%      estatismo: no ACTIVSg2000 (Texas sintetico), o criterio de
%      contingencia de recurso do ERCOT, de 2750 MW, que corresponde a
%      perda das duas maiores unidades; no IEEE 118, a maior unidade (N-1);
%   2) A reserva de cada unidade e limitada ao que o seu estatismo entrega
%      antes do primeiro estagio do alivio de carga por subfrequencia:
%      Rmax_g = (df_ERAC/f0)/R * Pmax_g, com R = 5 % e f_ERAC = 59,3 Hz,
%      isto e, 23,3 % do Pmax. Reserva alem disso nao atua a tempo;
%   3) O OPF AC escolhe quais unidades seguram a reserva (custo de reserva
%      pequeno e igual para todas, so para desempate). O resultado e uma
%      distribuicao heterogenea: unidades baratas carregadas e unidades
%      caras com folga;
%   4) IEEE 118 (PGLib, Pmin = 0): unidades despachadas abaixo de 1 MW sao
%      desligadas e o OPF e refeito (compromisso de unidades simplificado);
%   5) ACTIVSg2000: cargas, eolicas e bancos de capacitores das horas 5368
%      (pico) e 8504 (70 % do pico) iguais aos do estudo do CIA (arquivos
%      case_ACTIVSg2000_h5368/h8504), com margem termica de 5 % nas linhas.
%
% Os arquivos gerados guardam o ponto de operacao com os limites originais.


%% b) Limpeza de Figuras e Variaveis
clc
clear all
close all
format short
warning off


%% c) Pastas das Funcoes e dos Dados
% Rodar de dentro da pasta geracao; os casos sao gravados em ../dados.
% Requer o MATPOWER no path (runopf e toggle_reserves).
addpath(fullfile('..', 'funcoes'), fullfile('..', 'dados'));


%% d) Cabecalho
disp('....................................................................')
disp('.                                                                  .')
disp('.       PF-PFR PELO HELM: GERACAO DOS PONTOS DE OPERACAO (OPF)     .')
disp('.                                                                  .')
disp('....................................................................')


%% ________________________________________________________________________
%%
%%                       2. DADOS E PARAMETROS
%% ________________________________________________________________________

%% a) Parametros da reserva
% Estatismo das unidades
R     = 0.05;
% Primeiro estagio do alivio de carga por subfrequencia (Hz)
fERAC = 59.3;
f0    = 60;
% Fracao maxima do Pmax que pode ser reserva primaria
fMax  = ((f0 - fERAC)/f0)/R;
% Custo da reserva ($/MW), apenas para desempate
cRes  = 0.1;
% Opcoes do OPF
opcoes = mpoption('verbose', 0, 'out.all', 0);


%% ________________________________________________________________________
%%
%%                          3. IEEE 118 BARRAS
%% ________________________________________________________________________

%% a) Dados e unidades reais (os 35 compensadores sincronos tem Pmax = 0)
m   = pglib_opf_case118_ieee;
ng  = size(m.gen, 1);
uni = m.gen(:,9) > 0;

%% b) OPF com requisito de reserva e desligamento das unidades sem despacho
for it = 1:10
    opf = m;
    opf.reserves.zones = double(uni.' & opf.gen(:,8).' == 1);
    % Requisito N-1: maior unidade ligada
    opf.reserves.req  = max(opf.gen(uni & opf.gen(:,8) == 1, 9));
    opf.reserves.cost = cRes*ones(ng, 1);
    opf.reserves.qty  = fMax*opf.gen(:,9);
    opf = toggle_reserves(opf, 'on');
    res = runopf(opf, opcoes);
    if ~res.success
        error('Gera_Base_FPRP: o OPF do IEEE 118 nao convergiu.');
    end
    % Unidades reais despachadas abaixo de 1 MW sao desligadas
    desl = uni & m.gen(:,8) == 1 & res.gen(:,2) < 1;
    if ~any(desl)
        break
    end
    m.gen(desl, 8) = 0;
end

%% c) Ponto de operacao e gravacao
base = m;
base.bus(:,8)  = res.bus(:,8);
base.bus(:,9)  = res.bus(:,9);
base.gen(:,2)  = res.gen(:,2);
base.gen(:,3)  = res.gen(:,3);
[~, ig] = ismember(base.gen(:,1), base.bus(:,1));
base.gen(:,6)  = res.bus(ig, 8);
g  = find(uni & m.gen(:,8) == 1);
fo = (res.gen(g,9) - res.gen(g,2))./res.gen(g,9);
fprintf('IEEE 118: custo %.0f $/h | reserva %.0f MW | unidades %d | folga min %.2f, mediana %.2f, max %.2f | em Pmax: %d\n', ...
    res.f, sum(res.reserves.R(g)), numel(g), min(fo), median(fo), max(fo), sum(fo < 1e-3));
Grava_Caso('case118_fprp', base, false, ...
    {'IEEE 118 barras (PGLib-OPF v23.07) no ponto base do estudo do PF-PFR.', ...
     sprintf('OPF AC com requisito de reserva igual a maior unidade (%.0f MW),', opf.reserves.req), ...
     sprintf('reserva por unidade limitada a %.1f %% do Pmax; unidades sem despacho', 100*fMax), ...
     'desligadas (gerado por Gera_Base_FPRP.m). Limites originais.', ...
     'Fonte dos dados: IEEE PES Power Grid Library, licenca CC BY 4.0.'});


%% ________________________________________________________________________
%%
%%                          4. ACTIVSg2000
%% ________________________________________________________________________

%% a) Caso original (custos) e unidades com estatismo (gas e carvao)
m0   = case_ACTIVSg2000;
ng   = size(m0.gen, 1);
on   = m0.gen(:,8) == 1;
conv = on & ismember(m0.genfuel(:), {'ng', 'coal'}) & m0.gen(:,9) > m0.gen(:,10);
eol  = strcmp(m0.genfuel(:), 'wind');
parque = [2127 4153];
% Margem termica de 5 % (mesma do estudo do CIA)
rL   = 0.05;
% Requisito de reserva do ERCOT (MW)
req  = 2750;

for h = [5368 8504]

    %% b) Carga, compensacao e eolica da hora (iguais as do estudo do CIA)
    hb  = feval(sprintf('case_ACTIVSg2000_h%04d', h));
    opf = m0;
    opf.bus(:,3) = hb.bus(:,3);
    opf.bus(:,4) = hb.bus(:,4);
    opf.bus(:,6) = hb.bus(:,6);
    gp = find(ismember(m0.gen(:,1), parque) & eol & on);
    opf.gen(gp,2)  = hb.gen(gp,2);
    opf.gen(gp,9)  = hb.gen(gp,2);
    opf.gen(gp,10) = hb.gen(gp,2);
    opf.branch(:,6:8) = (1 - rL)*m0.branch(:,6:8);

    %% c) Requisito de reserva do sistema
    opf.reserves.zones = double(conv.');
    opf.reserves.req   = req;
    opf.reserves.cost  = cRes*ones(ng, 1);
    opf.reserves.qty   = fMax*opf.gen(:,9);
    opf = toggle_reserves(opf, 'on');

    %% d) OPF
    res = runopf(opf, opcoes);
    if ~res.success
        error('Gera_Base_FPRP: o OPF da hora %d nao convergiu.', h);
    end

    %% e) Ponto de operacao com os limites originais e gravacao
    base = m0;
    base.bus(:,3)  = opf.bus(:,3);
    base.bus(:,4)  = opf.bus(:,4);
    base.bus(:,6)  = opf.bus(:,6);
    base.bus(:,8)  = res.bus(:,8);
    base.bus(:,9)  = res.bus(:,9);
    base.gen(:,2)  = res.gen(:,2);
    base.gen(:,3)  = res.gen(:,3);
    [~, ig] = ismember(base.gen(:,1), base.bus(:,1));
    base.gen(:,6)  = res.bus(ig, 8);
    g  = find(conv);
    fo = (res.gen(g,9) - res.gen(g,2))./res.gen(g,9);
    fprintf('ACTIVSg2000 h%04d: custo %.0f $/h | reserva %.0f MW | unidades %d | folga min %.2f, mediana %.2f, max %.2f | em Pmax: %d\n', ...
        h, res.f, sum(res.reserves.R(g)), numel(g), min(fo), median(fo), max(fo), sum(fo < 1e-3));
    Grava_Caso(sprintf('case_ACTIVSg2000_fprp_h%04d', h), base, true, ...
        {sprintf('ACTIVSg2000 na hora %d da serie de 2016 (estudo do PF-PFR).', h), ...
         'Carga, eolicas e capacitores da hora iguais aos do estudo do CIA;', ...
         sprintf('OPF AC com requisito de reserva de %.0f MW (criterio do ERCOT) nas', req), ...
         sprintf('unidades a gas e carvao, reserva por unidade de ate %.1f %% do Pmax,', 100*fMax), ...
         sprintf('margem termica de %.0f %% (gerado por Gera_Base_FPRP.m). Limites originais.', 100*rL), ...
         'Fonte: Texas A&M ACTIVSg2000 (MATPOWER), licenca CC BY 4.0.'});
end
