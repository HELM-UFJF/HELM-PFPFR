%                                                                       JMJ
%
% Grava_Caso.m
% Grava um caso do MATPOWER (bus, gen, branch e, opcionalmente, genfuel)
% num arquivo de funcao ../dados/nome.m, com as linhas de texto do cabecalho.
%
% Entradas: nome (sem extensao), base (mpc), comFuel (grava genfuel),
%           texto (celula com as linhas do cabecalho)

function Grava_Caso(nome, base, comFuel, texto)
    fid = fopen(fullfile('..', 'dados', [nome, '.m']), 'w');
    fprintf(fid, 'function mpc = %s\n', nome);
    fprintf(fid, '%%%s  %s\n', upper(nome), texto{1});
    for i = 2:numel(texto)
        fprintf(fid, '%%   %s\n', texto{i});
    end
    fprintf(fid, '\nmpc.version = ''2'';\n');
    fprintf(fid, 'mpc.baseMVA = %g;\n\n', base.baseMVA);
    Escreve_Matriz(fid, 'mpc.bus', base.bus);
    Escreve_Matriz(fid, 'mpc.gen', base.gen);
    Escreve_Matriz(fid, 'mpc.branch', base.branch);
    if comFuel
        fprintf(fid, 'mpc.genfuel = {\n');
        fprintf(fid, '\t''%s'';\n', base.genfuel{:});
        fprintf(fid, '};\n');
    end
    fprintf(fid, 'end\n');
    fclose(fid);
end
