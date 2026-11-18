function [V, lam, Pg, info] = FPRP_Resolver(D, metodo, tol, nmax)
%FPRP_RESOLVER  Fluxo de potencia com regulacao primaria (PF-PFR) e limites.
%
%   [V, lam, Pg, info] = FPRP_Resolver(D, metodo, tol, nmax)
%
%   Laco externo unico, comum ao HELM e ao NR. A cada passagem, o PF-PFR e
%   resolvido para o estado atual das unidades e das barras PV, e as
%   verificacoes seguem esta prioridade:
%     a) Limites de P: as unidades cujo valor de estatismo
%        Pg_sch + lambda/R sai de [Pmin, Pmax] sao fixadas no limite e
%        deixam o estatismo. Com D.satur = 'todos' (padrao de
%        FPRP_Dados_MP), todas as violadoras de uma vez; com 'um', so a de
%        maior violacao;
%     b) Liberacao ('todos'): a unidade saturada cujo valor de estatismo
%        voltou para dentro dos limites volta a regular;
%     c) Limites de Q: as barras PV fora de [Qmin, Qmax] passam a PQ no
%        limite violado;
%     d) Retorno a PV: a barra em Qmax com |V| > Vesp, ou em Qmin com
%        |V| < Vesp, volta a controlar a tensao.
%   A cada mudanca de estado, o HELM recomeca da solucao germinal e o NR
%   da partida plana. Liberacoes e retornos sao limitados a 3 por unidade
%   e por barra (protecao contra ciclos).
%
%   D      : dados (FPRP_Dados_MP ou FPRP_Dados_MPC), em p.u.
%   metodo : 'HELM' ou 'NR'
%   tol    : tolerancia (criterio de parada do HELM ou residuo do NR)
%   nmax   : ordem maxima do HELM ou numero maximo de iteracoes do NR
%
%   V      : tensoes complexas
%   lam    : desvio de frequencia, lambda = -(f - f0)/f0
%   Pg     : geracao final das unidades com estatismo (ordem de D.ger)
%   info   : conv, n, ntot, passos, nsat, ativo, perda, qstate, nq,
%            Qbus e seq (eventos [tipo barra]: 1 Pmax, 2 Pmin, 3 Qmax,
%            4 Qmin, 5 retorno a PV, 6 liberacao)



    %% ____________________________________________________________________
    %%
    %%                           1. PREPARACAO
    %% ____________________________________________________________________

    %% a) Dados da Rede e das Cargas
    N    = D.N;
    Y    = sparse(D.Ysr) + spdiags(D.Ysh(:), 0, N, N);
    Pd   = D.Pd(:);
    Qd   = D.Qd(:);
    Pesp = D.Pg(:);
    usaHELM = strcmpi(metodo, 'HELM');


    %% b) Unidades com Estatismo
    % Potencia programada de cada unidade: D.Pug (varias unidades por
    % barra) ou, na falta dele, a geracao especificada da barra
    ger = D.ger(:);
    ng  = numel(ger);
    Rg  = D.Reg(:);
    Pmx = D.Pmax(:);
    Pmn = D.Pmin(:);

    if isfield(D, 'Pug')
        Pug = D.Pug(:);
    else
        Pug = Pesp(ger);
    end

    % Carga sensivel a frequencia
    if isfield(D, 'KL')
        KL = D.KL(:);
    else
        KL = zeros(N, 1);
    end


    %% c) Limites de Q (sem limite: -Inf e +Inf)
    if isfield(D, 'Qmin')
        Qmn = D.Qmin(:);
        Qmx = D.Qmax(:);
    else
        Qmn = -Inf*ones(N, 1);
        Qmx =  Inf*ones(N, 1);
    end


    %% d) Estrategia de Saturacao
    todos = isfield(D, 'satur') && strcmpi(D.satur, 'todos');


    %% e) Estado Inicial
    % Todas as unidades no estatismo e todas as barras PV regulando
    ativo  = true(ng, 1);
    Plim   = zeros(ng, 1);
    nlib   = zeros(ng, 1);
    qest   = zeros(N, 1);
    nvolta = zeros(N, 1);
    seq    = zeros(0, 2);
    tolP   = 1e-9;
    tolQ   = 1e-9;
    tolV   = 1e-9;
    passo  = 0;
    ntot   = 0;
    perda  = false;



    %% ____________________________________________________________________
    %%
    %%                       2. LACO DE ESTADOS
    %% ____________________________________________________________________

    while true
        passo = passo + 1;


        %% a) Injecoes e Rigidez do Estado Atual
        % A unidade saturada troca o seu Pg programado pelo limite
        Pnet = Pesp - Pd;
        dPs  = Plim(~ativo) - Pug(~ativo);
        Pnet = Pnet + accumarray(ger(~ativo), dPs, [N 1]);

        % Rigidez de cada barra: unidades no estatismo e carga
        beta = accumarray(ger(ativo), 1./Rg(ativo), [N 1]);
        beta = beta + KL;


        %% b) Tipos de Barra e Q Especificado
        tipo = D.tipo;
        tipo(qest ~= 0) = 1;
        Qnet = -Qd;
        Qnet(qest == 1)  = Qmx(qest == 1)  - Qd(qest == 1);
        Qnet(qest == -1) = Qmn(qest == -1) - Qd(qest == -1);


        %% c) PF-PFR do Estado Atual
        if usaHELM
            [V, lam, inn] = FPRP_HELM_Nucleo(D.Ysr, D.Ysh, tipo, D.Vesp, Pnet, Qnet, beta, D.ref, tol, nmax);
            ntot = ntot + inn.nrodadas;
        else
            V0 = [];
            if isfield(D, 'V0')
                V0 = D.V0;
            end
            [V, lam, inn] = FPRP_NR_Nucleo(Y, tipo, D.Vesp, D.ref, Pnet, Qnet, beta, tol, nmax, V0);
            ntot = ntot + inn.nit;
        end


        %% d) Geracoes Resultantes
        Pg = Pug;
        Pg(ativo)  = Pug(ativo) + lam./Rg(ativo);
        Pg(~ativo) = Plim(~ativo);

        % Q gerado em cada barra
        S    = V.*conj(Y*V);
        Qbus = imag(S) + Qd;
        Qbus(qest == 1)  = Qmx(qest == 1);
        Qbus(qest == -1) = Qmn(qest == -1);



        %% ________________________________________________________________
        %%
        %%                     3. VERIFICACAO DOS LIMITES
        %% ________________________________________________________________

        %% a) Limites de P (prioridade)
        viol = zeros(ng, 1);
        viol(ativo) = max(max(Pg(ativo) - Pmx(ativo), Pmn(ativo) - Pg(ativo)), 0);

        if max(viol) > tolP

            % 'todos': todas as violadoras; 'um': a de maior violacao
            if todos
                fixa = find(viol > tolP);
            else
                [~, fixa] = max(viol);
            end

            for g = fixa(:).'
                ativo(g) = false;

                if Pg(g) > Pmx(g)
                    Plim(g) = Pmx(g);
                    seq(end+1, :) = [1 ger(g)];
                else
                    Plim(g) = Pmn(g);
                    seq(end+1, :) = [2 ger(g)];
                end
            end

            % Sem unidades e sem carga sensivel: reserva esgotada
            if ~any(ativo) && sum(KL) == 0
                perda = true;
                break
            end

            continue
        end


        %% b) Liberacao de Unidades Saturadas ('todos')
        if todos
            Pdr  = Pug + lam./Rg;
            emMx = ~ativo & Plim == Pmx & Pdr < Pmx - tolP;
            emMn = ~ativo & Plim == Pmn & Pdr > Pmn + tolP;
            lib  = find((emMx | emMn) & nlib < 3);

            if ~isempty(lib)
                ativo(lib) = true;
                nlib(lib)  = nlib(lib) + 1;
                seq = [seq; 6*ones(numel(lib), 1), ger(lib)];
                continue
            end
        end


        %% c) Limites de Q nas Barras PV
        pv   = find(D.tipo == 2 & qest == 0);
        vmax = pv(Qbus(pv) > Qmx(pv) + tolQ);
        vmin = pv(Qbus(pv) < Qmn(pv) - tolQ);

        if ~isempty(vmax) || ~isempty(vmin)
            qest(vmax) = 1;
            qest(vmin) = -1;
            seq = [seq; 3*ones(numel(vmax), 1), vmax(:); 4*ones(numel(vmin), 1), vmin(:)];
            continue
        end


        %% d) Retorno ao Controle de Tensao
        lim   = find(qest ~= 0 & nvolta < 3);
        Vk    = abs(V(lim));
        volta = lim((qest(lim) == 1 & Vk > D.Vesp(lim) + tolV) | (qest(lim) == -1 & Vk < D.Vesp(lim) - tolV));

        if ~isempty(volta)
            qest(volta)   = 0;
            nvolta(volta) = nvolta(volta) + 1;
            seq = [seq; 5*ones(numel(volta), 1), volta(:)];
            continue
        end


        %% e) Estado Consistente
        break
    end



    %% ____________________________________________________________________
    %%
    %%                            4. SAIDAS
    %% ____________________________________________________________________

    %% a) Unidades Saturadas Sempre no Limite
    Pg(~ativo) = Plim(~ativo);


    %% b) Informacoes
    info.conv   = inn.conv;
    info.n      = inn.n;
    info.ntot   = ntot;
    info.passos = passo;
    info.nsat   = sum(~ativo);
    info.ativo  = ativo;
    info.perda  = perda;
    info.qstate = qest;
    info.nq     = sum(qest ~= 0);
    info.Qbus   = Qbus;
    info.seq    = seq;
end
