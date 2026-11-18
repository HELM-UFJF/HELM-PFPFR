function D = FPRP_Dados_MPC(mpc, usarQ)
%FPRP_DADOS_MPC  Estrutura de dados do PF-PFR a partir do formato estendido.
%
%   D = FPRP_Dados_MPC(mpc)
%   D = FPRP_Dados_MPC(mpc, usarQ)
%
%   Formato usado no tutorial de 3 barras e na reproducao do IEEE 39
%   barras: barras numeradas de 1 a N e uma unidade por barra, com o
%   estatismo nas colunas de mpc.bus:
%     1 barra, 2 tipo (1 PQ, 2 PV, 3 referencia), 3 Pd (MW), 4 Qd (Mvar),
%     6 Bsh (Mvar), 8 |V| especificado, 14 Pg (MW), 18 estatismo R na base
%     da maquina, 19 potencia nominal da maquina (MVA).
%   Ramos: 1 de, 2 para, 3 R, 4 X, 5 B total, 9 tap. Participam do
%   estatismo as barras PV e de referencia com R ~= 0. Os limites de P vem
%   de mpc.gen (colunas 9 e 10), se existir.
%
%   usarQ : true = limites de Q de mpc.bus (colunas 16 e 17) ou de mpc.gen
%           (colunas 4 e 5); padrao false (sem limites de Q)



    %% ____________________________________________________________________
    %%
    %%                            1. BARRAS
    %% ____________________________________________________________________

    %% a) Leitura
    if nargin < 2
        usarQ = false;
    end

    bus     = mpc.bus;
    branch  = mpc.branch;
    baseMVA = mpc.baseMVA;
    N       = size(bus, 1);


    %% b) Tipos, Tensoes e Injecoes
    D.N    = N;
    D.ref  = find(bus(:,2) == 3);
    D.tipo = ones(N, 1);
    D.tipo(bus(:,2) == 2) = 2;
    D.tipo(bus(:,2) == 3) = 0;
    D.Vesp = bus(:,8);
    D.Pd   = bus(:,3)/baseMVA;
    D.Qd   = bus(:,4)/baseMVA;
    D.Pg   = bus(:,14)/baseMVA;


    %% c) Unidades com Estatismo
    D.ger  = find((bus(:,2) == 3 | bus(:,2) == 2) & bus(:,18) ~= 0);
    D.Reg  = bus(D.ger,18).*baseMVA./bus(D.ger,19);
    D.Pmin = -Inf*ones(numel(D.ger), 1);
    D.Pmax =  Inf*ones(numel(D.ger), 1);

    if isfield(mpc, 'gen')
        for i = 1:numel(D.ger)
            g = find(mpc.gen(:,1) == D.ger(i));

            if ~isempty(g)
                D.Pmax(i) = mpc.gen(g(1),9)/baseMVA;
                D.Pmin(i) = mpc.gen(g(1),10)/baseMVA;
            end
        end
    end


    %% d) Limites de Q (a referencia nao e limitada)
    D.Qmin = -Inf*ones(N, 1);
    D.Qmax =  Inf*ones(N, 1);

    if usarQ
        for k = find(D.tipo == 2).'
            qa = -Inf;
            qb =  Inf;

            if size(bus, 2) >= 17 && bus(k,17) > bus(k,16)
                qa = bus(k,16)/baseMVA;
                qb = bus(k,17)/baseMVA;
            elseif isfield(mpc, 'gen')
                g = find(mpc.gen(:,1) == k);

                if ~isempty(g) && mpc.gen(g(1),4) > mpc.gen(g(1),5)
                    qa = mpc.gen(g(1),5)/baseMVA;
                    qb = mpc.gen(g(1),4)/baseMVA;
                end
            end

            D.Qmin(k) = qa;
            D.Qmax(k) = qb;
        end

        % Valores como -99999 e 99999 significam sem limite
        D.Qmin(D.Qmin < -500) = -Inf;
        D.Qmax(D.Qmax >  500) =  Inf;
    end



    %% ____________________________________________________________________
    %%
    %%                     2. MATRIZ DE ADMITANCIAS
    %% ____________________________________________________________________

    %% a) Ramos (tap no lado DE; susceptancia de linha nos dois lados)
    k   = branch(:,1);
    m   = branch(:,2);
    y   = 1./(branch(:,3) + 1j*branch(:,4));
    a   = branch(:,9);
    shl = 1j*branch(:,5)/2;

    Y = sparse([k; m; k; m], [k; m; m; k], [y./a.^2 + shl; y + shl; -y./a; -y./a], N, N);
    Y = Y + spdiags(1j*bus(:,6)/baseMVA, 0, N, N);


    %% b) Separacao: Ysh = Soma das Linhas (solucao germinal V = 1)
    D.Ysh = full(sum(Y, 2));
    D.Ysr = Y - spdiags(D.Ysh, 0, N, N);
end
