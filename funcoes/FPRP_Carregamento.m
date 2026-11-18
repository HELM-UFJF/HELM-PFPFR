function S = FPRP_Carregamento(mpc, D, V)
%FPRP_CARREGAMENTO  Carregamento dos ramos em relacao a capacidade (rateA).
%
%   S = FPRP_Carregamento(mpc, D, V)
%
%   Calcula o fluxo aparente nos dois terminais de cada ramo em servico com
%   capacidade definida (rateA > 0), na convencao do MATPOWER, e devolve o
%   maior dos dois em % da capacidade.
%
%   mpc : caso no formato padrao do MATPOWER
%   D   : estrutura de FPRP_Dados_MP (numeracao interna)
%   V   : tensoes complexas
%   S   : carregamento de cada ramo (%), na ordem de mpc.branch filtrada



    %% a) Ramos em Servico com Capacidade
    br = mpc.branch(mpc.branch(:,11) == 1 & mpc.branch(:,6) > 0, :);
    [~, f] = ismember(br(:,1), D.num);
    [~, t] = ismember(br(:,2), D.num);


    %% b) Parametros (tap nulo = 1)
    ys = 1./(br(:,3) + 1j*br(:,4));
    a  = br(:,9);
    a(a == 0) = 1;
    tc  = a.*exp(1j*pi/180*br(:,10));
    Ytt = ys + 1j*br(:,5)/2;


    %% c) Correntes e Potencias nos Terminais
    If = Ytt./(tc.*conj(tc)).*V(f) - ys./conj(tc).*V(t);
    It = -ys./tc.*V(f) + Ytt.*V(t);
    Sf = abs(V(f).*conj(If));
    St = abs(V(t).*conj(It));


    %% d) Carregamento (%)
    S = 100*max(Sf, St)*mpc.baseMVA./br(:,6);
end
