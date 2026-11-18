function [mpc, opc, nome] = FPRP_Caso(sis)
%FPRP_CASO  Ponto de operacao do estudo e opcoes de FPRP_Dados_MP.
%
%   [mpc, opc, nome] = FPRP_Caso(sis)
%
%   sis : '118'  - IEEE 118 barras (case118_fprp)
%         'pico' - ACTIVSg2000, hora 5368 (case_ACTIVSg2000_fprp_h5368)
%         'leve' - ACTIVSg2000, hora 8504 (case_ACTIVSg2000_fprp_h8504)
%
%   Regulam as unidades que podem aumentar a geracao: no IEEE 118, as
%   unidades reais (Pmax > 0); no ACTIVSg2000, as unidades a gas e a
%   carvao. Unidades despachadas em Pmax nao regulam em subfrequencia. No
%   ACTIVSg2000, a referencia passa a barra 2057 (maior unidade a
%   carvao), pois as unidades nucleares das barras 7098 e 7099 sao as
%   contingencias do estudo. A referencia e ajustada para 60 Hz no ponto
%   base (FPRP_Ajusta_Ref).



    %% a) Caso e Candidatas ao Estatismo
    opc = struct();

    switch sis
        case '118'
            mpc  = case118_fprp;
            cand = mpc.gen(:,9) > 0;
            nome = 'IEEE 118';

        case 'pico'
            mpc  = case_ACTIVSg2000_fprp_h5368;
            cand = ismember(mpc.genfuel(:), {'ng', 'coal'});
            opc.ref = 2057;
            nome = 'ACTIVSg2000, pico';

        case 'leve'
            mpc  = case_ACTIVSg2000_fprp_h8504;
            cand = ismember(mpc.genfuel(:), {'ng', 'coal'});
            opc.ref = 2057;
            nome = 'ACTIVSg2000, leve';

        otherwise
            error('FPRP_Caso: sistema desconhecido (%s).', sis);
    end


    %% b) Unidades com Estatismo (as que podem subir)
    opc.estatismo = cand & mpc.gen(:,2) < mpc.gen(:,9) - 1e-3;


    %% c) Referencia em 60 Hz no Ponto Base
    mpc = FPRP_Ajusta_Ref(mpc, opc);
end
