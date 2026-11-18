function D = FPRP_Dados_MP(mpc, opc)
%FPRP_DADOS_MP  Estrutura de dados do PF-PFR a partir de um caso do MATPOWER.
%
%   D = FPRP_Dados_MP(mpc)
%   D = FPRP_Dados_MP(mpc, opc)
%
%   Converte um caso no formato padrao do MATPOWER (bus, gen, branch e,
%   opcionalmente, genfuel), unidade a unidade:
%     - barras renumeradas de 1 a N (D.num guarda a numeracao original);
%     - matriz nodal na convencao do MATPOWER: estado do ramo, tap (0 = 1),
%       defasagem, susceptancia de linha e shunts de barra (Gs + jBs);
%     - unidades ligadas (gen(:,8) = 1), varias por barra; a barra PV sem
%       unidade ligada passa a PQ;
%     - unidades com estatismo: opc.estatismo ou, por padrao, as ligadas
%       com Pmax > Pmin; estatismo opc.R na base da unidade, com potencia
%       nominal S_g = Pmax (Reg = R*baseMVA/Pmax, em p.u. do sistema);
%     - limites de Q por barra = soma dos limites das unidades ligadas; a
%       referencia nao e limitada;
%     - carga sensivel a frequencia K_L = kap*Pd.
%
%   opc.estatismo : vetor logico (linhas de mpc.gen) das unidades com
%                   estatismo
%   opc.R         : estatismo na base da unidade (padrao 0.05)
%   opc.kap       : sensibilidade da carga a frequencia (padrao 1)
%   opc.ref       : numero ORIGINAL da barra de referencia (padrao: a do
%                   caso)
%   opc.satur     : estrategia de saturacao, 'todos' (padrao) ou 'um'
%
%   D contem os campos usados por FPRP_Resolver e ainda Pug (P programado
%   de cada unidade com estatismo), gidx (linha de mpc.gen de cada unidade
%   com estatismo), num (numero original das barras) e baseMVA.



    %% ____________________________________________________________________
    %%
    %%                      1. OPCOES E NUMERACAO
    %% ____________________________________________________________________

    %% a) Valores Padrao
    if nargin < 2
        opc = struct();
    end

    if ~isfield(opc, 'R')
        opc.R = 0.05;
    end

    if ~isfield(opc, 'kap')
        opc.kap = 1;
    end

    if ~isfield(opc, 'satur')
        opc.satur = 'todos';
    end

    bus     = mpc.bus;
    gen     = mpc.gen;
    branch  = mpc.branch;
    baseMVA = mpc.baseMVA;
    N       = size(bus, 1);


    %% b) Numeracao Interna (1 a N)
    num  = bus(:,1);
    mapa = sparse(num, 1, (1:N).', max(num), 1);
    gb   = full(mapa(gen(:,1)));
    on   = gen(:,8) == 1;



    %% ____________________________________________________________________
    %%
    %%                            2. BARRAS
    %% ____________________________________________________________________

    %% a) Tipos de Barra
    % Barras PV (ou a antiga referencia) sem unidade ligada passam a PQ
    temGer = accumarray(gb(on), 1, [N 1]) > 0;
    tipo   = ones(N, 1);
    tipo((bus(:,2) == 2 | bus(:,2) == 3) & temGer) = 2;

    if isfield(opc, 'ref')
        ref = full(mapa(opc.ref));
    else
        ref = find(bus(:,2) == 3);
    end

    if ~temGer(ref)
        error('FPRP_Dados_MP: a barra de referencia nao tem unidade ligada.');
    end
    tipo(ref) = 0;


    %% b) Cargas
    D.N    = N;
    D.ref  = ref;
    D.tipo = tipo;
    D.Pd   = bus(:,3)/baseMVA;
    D.Qd   = bus(:,4)/baseMVA;
    D.KL   = opc.kap*D.Pd;


    %% c) Tensao Especificada (das unidades ligadas na barra)
    Vesp = bus(:,8);
    gOn  = find(on);
    Vesp(gb(gOn)) = gen(gOn,6);
    D.Vesp = Vesp;


    %% d) Geracao por Barra e Limites de Q
    D.Pg   = accumarray(gb(on), gen(on,2), [N 1])/baseMVA;
    D.Qmax = accumarray(gb(on), gen(on,4), [N 1])/baseMVA;
    D.Qmin = accumarray(gb(on), gen(on,5), [N 1])/baseMVA;

    % Sem limite nas barras PQ e na referencia
    semQ = tipo ~= 2;
    D.Qmax(semQ) =  Inf;
    D.Qmin(semQ) = -Inf;



    %% ____________________________________________________________________
    %%
    %%                   3. UNIDADES COM ESTATISMO
    %% ____________________________________________________________________

    %% a) Selecao
    if isfield(opc, 'estatismo')
        est = opc.estatismo(:) & on;
    else
        est = on & gen(:,9) > gen(:,10);
    end

    % Unidade em barra que deixou de ser PV ou referencia nao regula
    est  = est & tipo(gb) ~= 1;
    gidx = find(est);


    %% b) Dados por Unidade (p.u. na base do sistema)
    D.ger  = gb(gidx);
    D.gidx = gidx;
    D.Pug  = gen(gidx,2)/baseMVA;
    D.Pmax = gen(gidx,9)/baseMVA;
    D.Pmin = gen(gidx,10)/baseMVA;
    D.Reg  = opc.R*baseMVA./gen(gidx,9);



    %% ____________________________________________________________________
    %%
    %%                     4. MATRIZ DE ADMITANCIAS
    %% ____________________________________________________________________

    %% a) Ramos em Servico (convencao do MATPOWER)
    lig = branch(:,11) == 1;
    br  = branch(lig,:);
    f   = full(mapa(br(:,1)));
    t   = full(mapa(br(:,2)));
    ys  = 1./(br(:,3) + 1j*br(:,4));
    bc  = br(:,5);

    % Tap complexo (tap nulo = 1)
    a = br(:,9);
    a(a == 0) = 1;
    tc = a.*exp(1j*pi/180*br(:,10));

    Ytt = ys + 1j*bc/2;
    Yff = Ytt./(tc.*conj(tc));
    Yft = -ys./conj(tc);
    Ytf = -ys./tc;


    %% b) Montagem Esparsa com os Shunts de Barra
    Ysh0 = (bus(:,5) + 1j*bus(:,6))/baseMVA;
    Y = sparse([f; f; t; t], [f; t; f; t], [Yff; Yft; Ytf; Ytt], N, N);
    Y = Y + spdiags(Ysh0, 0, N, N);


    %% c) Separacao: Ysh = Soma das Linhas (solucao germinal V = 1)
    D.Ysh = full(sum(Y, 2));
    D.Ysr = Y - spdiags(D.Ysh, 0, N, N);


    %% d) Informacoes Adicionais
    D.num     = num;
    D.baseMVA = baseMVA;
    D.satur   = opc.satur;
end
