function [mpc, D, dP] = FPRP_Ajusta_Ref(mpc, opc)
%FPRP_AJUSTA_REF  Ajusta a referencia para a frequencia nominal no ponto base.
%
%   [mpc, D, dP] = FPRP_Ajusta_Ref(mpc, opc)
%
%   O PF-PFR do ponto base e resolvido pelo NR, sem limites de P. O desvio
%   lambda obtido corresponde a um desbalanco de beta_tot*lambda, que e
%   somado a geracao programada das unidades da referencia (repartido pelo
%   Pmax). Repete-se ate |lambda| < 1e-12. Num ponto vindo de um OPF, o
%   ajuste e da ordem do residuo do OPF.
%
%   mpc : caso base no formato padrao do MATPOWER
%   opc : as mesmas opcoes de FPRP_Dados_MP
%   dP  : ajuste total (MW)



    %% ____________________________________________________________________
    %%
    %%                    1. UNIDADES DA REFERENCIA
    %% ____________________________________________________________________

    %% a) Barra de Referencia (numeracao original)
    D    = FPRP_Dados_MP(mpc, opc);
    nref = D.num(D.ref);


    %% b) Unidades Ligadas na Referencia
    gr = find(mpc.gen(:,1) == nref & mpc.gen(:,8) == 1);
    w  = mpc.gen(gr,9)/sum(mpc.gen(gr,9));
    dP = 0;



    %% ____________________________________________________________________
    %%
    %%                            2. AJUSTE
    %% ____________________________________________________________________

    for it = 1:10

        %% a) PF-PFR sem Limites de P
        Dn = D;
        Dn.Pmin(:) = -Inf;
        Dn.Pmax(:) =  Inf;
        [~, lam, ~, inf0] = FPRP_Resolver(Dn, 'NR', 1e-12, 30);

        if ~inf0.conv
            error('FPRP_Ajusta_Ref: o NR nao convergiu no ponto base.');
        end

        if abs(lam) < 1e-12
            break
        end


        %% b) Desbalanco Equivalente Somado a Referencia (MW)
        bTot = sum(1./D.Reg) + sum(D.KL);
        dPi  = lam*bTot*mpc.baseMVA;
        mpc.gen(gr,2) = mpc.gen(gr,2) + w*dPi;
        dP = dP + dPi;
        D  = FPRP_Dados_MP(mpc, opc);
    end
end
