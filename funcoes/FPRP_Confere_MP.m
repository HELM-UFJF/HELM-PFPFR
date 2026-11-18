function c = FPRP_Confere_MP(mpc, D, V, lam, Pg, kap, qst)
%FPRP_CONFERE_MP  Conferencia de uma solucao do PF-PFR pelo MATPOWER.
%
%   c = FPRP_Confere_MP(mpc, D, V, lam, Pg, kap, qst)
%
%   O estado final do PF-PFR e reproduzido no runpf, sem o estatismo:
%     - as unidades com estatismo recebem a geracao final;
%     - as cargas recebem o valor efetivo P0_L*(1 - kap*lambda), com Q
%       constante;
%     - as barras no limite de Q passam a PQ, com cada unidade no seu
%       proprio limite;
%     - o runpf e resolvido sem laco de limites.
%   A barra de referencia do MATPOWER fecha o balanco; se o PF-PFR estiver
%   certo, ela fornece a mesma potencia da unidade de referencia, as
%   tensoes coincidem e as barras que ficaram PV respeitam os limites de Q.
%
%   mpc, D : caso do cenario e a sua estrutura (FPRP_Dados_MP)
%   V, lam : solucao do PF-PFR
%   Pg     : geracao final das unidades com estatismo (ordem de D.gidx)
%   kap    : sensibilidade da carga a frequencia
%   qst    : estado de Q de cada barra (info.qstate: 0, +1 ou -1)
%
%   c.ok    : true se o runpf convergiu
%   c.dV    : max |V_MP - V| (p.u.), angulos referidos a referencia
%   c.dPref : P da referencia no MATPOWER menos o do PF-PFR (MW)
%   c.violQ : maior violacao de Q nas barras que ficaram PV (p.u.)



    %% ____________________________________________________________________
    %%
    %%                     1. CASO SEM ESTATISMO
    %% ____________________________________________________________________

    %% a) Geracao Final e Cargas Efetivas
    m = mpc;
    m.gen(D.gidx,2) = Pg(:)*m.baseMVA;
    m.bus(:,3) = m.bus(:,3)*(1 - kap*lam);


    %% b) Tensoes Especificadas Iguais as do PF-PFR
    on = find(m.gen(:,8) == 1);
    [~, ib] = ismember(m.gen(on,1), D.num);
    m.gen(on,6) = D.Vesp(ib);


    %% c) Tipos de Barra e Estado Final de Q
    tp = D.tipo;
    tp(D.tipo == 0) = 3;

    for k = find(qst(:) ~= 0).'
        gk = find(m.gen(:,1) == D.num(k) & m.gen(:,8) == 1);
        tp(k) = 1;

        if qst(k) > 0
            m.gen(gk,3) = m.gen(gk,4);
        else
            m.gen(gk,3) = m.gen(gk,5);
        end
    end
    m.bus(:,2) = tp;


    %% d) Fluxo de Potencia sem Laco de Limites
    op = mpoption('verbose', 0, 'out.all', 0, 'pf.enforce_q_lims', 0, 'pf.tol', 1e-11);
    r  = runpf(m, op);
    c.ok = r.success == 1;



    %% ____________________________________________________________________
    %%
    %%                          2. COMPARACAO
    %% ____________________________________________________________________

    %% a) Tensoes
    Vm   = r.bus(:,8).*exp(1j*pi/180*(r.bus(:,9) - r.bus(D.ref,9)));
    c.dV = max(abs(Vm - V(:)));


    %% b) Potencia da Referencia
    gr = find(m.gen(:,1) == D.num(D.ref) & m.gen(:,8) == 1);
    c.dPref = sum(r.gen(gr,2)) - sum(m.gen(gr,2));


    %% c) Limites de Q das Barras que Ficaram PV
    pv  = find(tp == 2);
    Qpv = zeros(numel(pv), 1);

    for i = 1:numel(pv)
        gk = m.gen(:,1) == D.num(pv(i)) & m.gen(:,8) == 1;
        Qpv(i) = sum(r.gen(gk,3))/m.baseMVA;
    end
    c.violQ = max([0; Qpv - D.Qmax(pv); D.Qmin(pv) - Qpv]);
end
