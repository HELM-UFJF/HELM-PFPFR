function [D, mpc] = Wang39_Dados(saida30)
%WANG39_DADOS  IEEE 39 barras no formato de FPRP_Dados_MPC (Wang et al.).
%
%   [D, mpc] = Wang39_Dados(saida30)
%
%   Converte o case39 do MATPOWER para o formato estendido de
%   FPRP_Dados_MPC, com a opcao de retirar o gerador da barra 30, como no
%   caso de Wang et al. (Energy Proceedings, 2020).
%     - tap nulo no MATPOWER significa tap = 1;
%     - |V| especificado das barras com gerador = Vg (mpc.gen, coluna 6);
%     - o estatismo de cada gerador e provisorio (R = 1 p.u.): a calibracao
%       e feita em Main_Wang39;
%     - os limites de P vem de mpc.gen (colunas 9 e 10).
%
%   saida30 : true = gerador da barra 30 fora de servico



    %% ____________________________________________________________________
    %%
    %%                     1. CONVERSAO DO CASO
    %% ____________________________________________________________________

    %% a) Dados do MATPOWER
    m  = case39();
    nb = size(m.bus, 1);


    %% b) Barras no Formato Estendido
    % Colunas: 1 barra, 2 tipo, 3 Pd, 4 Qd, 6 Bsh, 8 |V|, 14 Pg, 18 R,
    % 19 potencia nominal da maquina
    bus = zeros(nb, 19);
    bus(:,1) = m.bus(:,1);
    bus(:,2) = m.bus(:,2);
    bus(:,3) = m.bus(:,3);
    bus(:,4) = m.bus(:,4);
    bus(:,6) = m.bus(:,6);
    bus(:,8) = m.bus(:,8);

    for g = 1:size(m.gen, 1)
        k = m.gen(g,1);
        bus(k,8)  = m.gen(g,6);
        bus(k,14) = m.gen(g,2);
        bus(k,18) = 1;
        bus(k,19) = m.baseMVA;
    end


    %% c) Ramos (tap nulo = 1)
    br = m.branch(:,1:9);
    br(br(:,9) == 0, 9) = 1;

    mpc.baseMVA = m.baseMVA;
    mpc.bus     = bus;
    mpc.branch  = br;
    mpc.gen     = m.gen;



    %% ____________________________________________________________________
    %%
    %%                 2. SAIDA DO GERADOR DA BARRA 30
    %% ____________________________________________________________________

    if saida30
        mpc.bus(30,2)  = 1;
        mpc.bus(30,14) = 0;
        mpc.bus(30,18) = 0;
        mpc.gen(mpc.gen(:,1) == 30, :) = [];
    end



    %% ____________________________________________________________________
    %%
    %%                        3. ESTRUTURA D
    %% ____________________________________________________________________

    D = FPRP_Dados_MPC(mpc, false);
end
