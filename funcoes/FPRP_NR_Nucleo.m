function [V, lam, info] = FPRP_NR_Nucleo(Y, tipo, Vesp, ref, P, Q, beta, tol, nitmax, V0)
%FPRP_NR_NUCLEO  Newton-Raphson do PF-PFR para um estado fixo.
%
%   [V, lam, info] = FPRP_NR_Nucleo(Y, tipo, Vesp, ref, P, Q, beta, tol, nitmax, V0)
%
%   Incognitas: angulos theta (todas as barras, exceto a referencia),
%   modulos |V| das barras PQ e o desvio de frequencia lambda. Equacoes:
%   balanco de P em todas as barras, com a injecao P_k + beta_k*lambda, e
%   balanco de Q nas barras PQ. O Jacobiano e o classico, esparso, com uma
%   coluna a mais (-beta) para a frequencia; e refatorado a cada iteracao.
%
%   Y        : matriz nodal (rede e shunts), esparsa
%   tipo     : 0 referencia, 1 PQ, 2 PV
%   Vesp     : modulos de tensao especificados
%   ref      : barra de referencia
%   P, Q     : injecoes liquidas especificadas em lambda = 0 (p.u.)
%   beta     : rigidez de cada barra (1/R das unidades ativas + K_L)
%   tol      : residuo maximo (p.u.)
%   nitmax   : numero maximo de iteracoes
%   V0       : chute inicial opcional ([] = partida plana: theta = 0,
%              |V| = Vesp nas barras PV e na referencia, 1 nas PQ)
%
%   V, lam   : tensoes complexas e desvio de frequencia
%   info     : nit e n (iteracoes) e conv



    %% ____________________________________________________________________
    %%
    %%                           1. PREPARACAO
    %% ____________________________________________________________________

    %% a) Conjuntos de Barras
    N    = numel(tipo);
    P    = P(:);
    Q    = Q(:);
    beta = beta(:);
    PQ   = find(tipo(:) == 1);
    nRef = setdiff((1:N).', ref);
    nA   = numel(nRef);
    nV   = numel(PQ);


    %% b) Condicao Inicial
    Vm = Vesp(:);
    Vm(PQ) = 1;
    Va = zeros(N, 1);

    if nargin >= 10 && ~isempty(V0)
        Vm(PQ) = abs(V0(PQ));
        Va     = angle(V0(:));
    end

    lam  = 0;
    conv = false;
    nit  = 0;



    %% ____________________________________________________________________
    %%
    %%                           2. ITERACOES
    %% ____________________________________________________________________

    for it = 1:nitmax

        %% a) Residuos
        V = Vm.*exp(1j*Va);
        S = V.*conj(Y*V);

        % Residuo de P (todas as barras) e de Q (barras PQ)
        f1 = real(S) - (P + beta*lam);
        f2 = imag(S(PQ)) - Q(PQ);
        F  = [f1; f2];

        if max(abs(F)) < tol
            conv = true;
            break
        end


        %% b) Jacobiano Esparso
        % Derivadas de S em relacao a theta e a |V|
        Vd = spdiags(V, 0, N, N);
        Id = spdiags(Y*V, 0, N, N);
        En = spdiags(V./abs(V), 0, N, N);
        dS_dVa = 1j*Vd*conj(Id - Y*Vd);
        dS_dVm = Vd*conj(Y*En) + conj(Id)*En;

        J = [real(dS_dVa(:,nRef)),  real(dS_dVm(:,PQ)),  sparse(-beta);
             imag(dS_dVa(PQ,nRef)), imag(dS_dVm(PQ,PQ)), sparse(nV, 1)];


        %% c) Atualizacao
        dx = -J\F;
        Va(nRef) = Va(nRef) + dx(1:nA);
        Vm(PQ)   = Vm(PQ) + dx(nA+1:nA+nV);
        lam      = lam + dx(end);
        nit      = it;
    end



    %% ____________________________________________________________________
    %%
    %%                            3. SAIDAS
    %% ____________________________________________________________________

    info.nit  = nit;
    info.n    = nit;
    info.conv = conv;
end
