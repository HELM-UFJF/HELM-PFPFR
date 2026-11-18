function Escreve_Matriz(fid, nome, M)
%ESCREVE_MATRIZ  Grava a matriz M no formato de caso do MATPOWER.


    %% a) Linhas da Matriz
    fprintf(fid, '%s = [\n', nome);

    for i = 1:size(M,1)
        fprintf(fid, '\t%s;\n', sprintf('%.10g\t', M(i,:)));
    end

    fprintf(fid, '];\n\n');
end
