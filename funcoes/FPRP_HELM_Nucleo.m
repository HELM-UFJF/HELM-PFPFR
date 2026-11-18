function [V1, lam1, info] = FPRP_HELM_Nucleo(Ysr, Ysh, tipo, Vesp, P, Q, beta, ref, xi, nmax, nestag)
%FPRP_HELM_NUCLEO  HELM do PF-PFR para um estado fixo das unidades e barras.
%
%   [V1, lam1, info] = FPRP_HELM_Nucleo(Ysr, Ysh, tipo, Vesp, P, Q, beta, ref, xi, nmax, nestag)
%
%   O desvio de frequencia entra como mais uma serie, lambda(gamma), com
%   lambda(0) = 0. Cada barra com beta_k ~= 0 (unidade com estatismo ou
%   carga sensivel a frequencia) recebe o termo beta_k*lambda(gamma) na
%   equacao de potencia ativa. A equacao nodal da referencia e substituida
%   pela linha de controle de frequencia, que e o balanco de potencia ativa
%   do sistema:
%       Re{sum_km conj(Ysr_km) V_k V_m* + gamma sum_k conj(Ysh_k) |V_k|^2}
%           = gamma sum_k P_k + beta_tot lambda(gamma)
%   A matriz de recursao, (2N+1)x(2N+1) e esparsa, nao depende da ordem n:
%   e fatorada uma unica vez. As series sao avaliadas em gamma = 1 pelo
%   algoritmo epsilon de Wynn incremental (Wynn_Passo), com todas as series
%   empilhadas num unico vetor.
%
%   Ysr, Ysh : matriz nodal sem shunts (linhas de soma nula) e shunts
%   tipo     : 0 referencia (Vtheta), 1 PQ, 2 PV
%   Vesp     : modulos de tensao especificados
%   P, Q     : injecoes liquidas especificadas em lambda = 0 (p.u.)
%   beta     : rigidez de cada barra (1/R das unidades ativas + K_L)
%   ref      : barra de referencia
%   xi       : tolerancia do criterio de parada
%   nmax     : ordem maxima das series
%   nestag   : ordens sem melhora ate a parada por estagnacao (padrao 15)
%
%   V1, lam1 : tensoes e desvio de frequencia em gamma = 1
%   info     : V, lam, Qh (coeficientes), n (ordem da estimativa), nrodadas
%              (ordens calculadas), conv, estag e hist (erro por ordem)



    %% ____________________________________________________________________
    %%
    %%                           1. PREPARACAO
    %% ____________________________________________________________________

    %% a) Dimensoes e Conjuntos de Barras
    if nargin < 11 || isempty(nestag)
        nestag = 15;
    end

    N   = numel(Ysh);
    nu  = 2*N + 1;
    Ysr = sparse(Ysr);
    Ysh = Ysh(:);
    P   = P(:);
    Q   = Q(:);
    beta = beta(:);
    G   = real(Ysr);
    B   = imag(Ysr);

    % Barras PV (a referencia tem linhas proprias)
    PV  = find(tipo(:) == 2);
    PV  = PV(PV ~= ref);

    % Barras sensiveis a frequencia (fora a referencia)
    H   = find(beta ~= 0);
    H   = H(H ~= ref);

    % Menor ordem em que a parada e avaliada
    nmin = 4;


    %% b) Coeficientes da Linha de Controle de Frequencia
    % d_k = conj(sum_m Ysr_km) + sum_m Ysr_mk; nulos para Y simetrica
    d   = conj(full(sum(Ysr, 2))) + full(sum(Ysr, 1)).';
    cre = real(d);
    cim = -imag(d);



    %% ____________________________________________________________________
    %%
    %%                     2. MATRIZ DE RECURSAO
    %% ____________________________________________________________________

    %% a) Blocos da Rede
    % Linhas 2k-1 e 2k: partes real e imaginaria da equacao da barra k;
    % colunas 2m-1 e 2m: Re{V_m} (ou Q_m, se m e PV) e Im{V_m}
    [ik, im, g] = find(G);
    [jk, jm, b] = find(B);
    Ia = [2*ik-1; 2*jk-1; 2*jk; 2*ik];
    Ja = [2*im-1; 2*jm; 2*jm-1; 2*im];
    Va = [g; -b; b; g];


    %% b) Colunas das Barras PV e Linhas da Referencia
    % A coluna de Re{V} das barras PV vira a coluna de Q (removida aqui)
    ePV = false(nu, 1);
    ePV(2*PV-1) = true;
    manter = ~ePV(Ja);

    % As linhas da referencia sao substituidas por linhas de identidade
    manter = manter & Ia ~= 2*ref-1 & Ia ~= 2*ref;
    Ia = Ia(manter);
    Ja = Ja(manter);
    Va = Va(manter);


    %% c) Coluna de Q das Barras PV
    % Unico elemento nao nulo: 1 na equacao imaginaria da propria barra
    Ia = [Ia; 2*PV];
    Ja = [Ja; 2*PV-1];
    Va = [Va; ones(numel(PV), 1)];


    %% d) Coluna de lambda
    % -beta_k na equacao real de cada barra sensivel a frequencia
    Ia = [Ia; 2*H-1];
    Ja = [Ja; nu*ones(numel(H), 1)];
    Va = [Va; -beta(H)];


    %% e) Linhas de Identidade da Referencia
    Ia = [Ia; 2*ref-1; 2*ref];
    Ja = [Ja; 2*ref-1; 2*ref];
    Va = [Va; 1; 1];


    %% f) Linha de Controle de Frequencia
    % Re{V} das barras PV e dado, nao incognita
    cpar = cre;
    cpar(PV) = 0;
    Ia = [Ia; nu*ones(2*N+1, 1)];
    Ja = [Ja; (1:2:2*N).'; (2:2:2*N).'; nu];
    Va = [Va; cpar; cim; -sum(beta)];


    %% g) Montagem e Fatoracao (uma unica vez)
    A = sparse(Ia, Ja, Va, nu, nu);
    [L, U, Pm, Qc] = lu(A);



    %% ____________________________________________________________________
    %%
    %%                       3. SOLUCAO GERMINAL
    %% ____________________________________________________________________

    %% a) Armazenamento (coluna n+1 guarda a ordem n)
    V   = zeros(N, nmax+1);
    W   = zeros(N, nmax+1);
    Qh  = zeros(N, nmax+1);
    Ic  = zeros(N, nmax+1);
    lam = zeros(1, nmax+1);
    V(:,1) = 1;
    W(:,1) = 1;
    Sesp   = P - 1j*Q;


    %% b) Controle da Parada
    hist    = zeros(1, 0);
    conv    = false;
    estag   = false;
    Vp      = V(:,1);
    lamp    = 0;
    Vp_ant  = [];
    lamp_ant = 0;
    errbest = Inf;
    nbest   = 0;
    Vbest   = Vp;
    lambest = 0;


    %% c) Tabela Epsilon de Wynn (somas parciais de ordem 0)
    sS = [V(:,1); lam(1)];
    [wS, ~] = Wynn_Passo([], sS);



    %% ____________________________________________________________________
    %%
    %%                     4. RECURSAO DAS SERIES
    %% ____________________________________________________________________

    for n = 1:nmax

        %% a) Parte Real Conhecida das Barras PV (|V|^2 = Vesp^2)
        Vre = zeros(N, 1);

        if n == 1
            Vre(PV) = (Vesp(PV).^2 - 1)/2;
        else
            Vre(PV) = -0.5*real(sum(V(PV,2:n).*conj(V(PV,n:-1:2)), 2));
        end


        %% b) Lado Direito das Barras
        R = Sesp.*conj(W(:,n)) - Ysh.*V(:,n);

        % Barras PV: Q e incognita (serie)
        if n > 1
            sq = sum(Qh(PV,2:n).*conj(W(PV,n:-1:2)), 2);
        else
            sq = zeros(numel(PV), 1);
        end
        R(PV) = P(PV).*conj(W(PV,n)) - 1j*sq - Ysh(PV).*V(PV,n);

        % Barras sensiveis a frequencia: convolucao de lambda com W*
        if n > 1 && ~isempty(H)
            R(H) = R(H) + beta(H).*(conj(W(H,n:-1:2))*lam(2:n).');
        end

        rhs = zeros(nu, 1);
        rhs(1:2:2*N) = real(R) - G(:,PV)*Vre(PV);
        rhs(2:2:2*N) = imag(R) - B(:,PV)*Vre(PV);


        %% c) Lado Direito da Referencia
        % V_ref(gamma) = 1 + gamma*(Vesp - 1)
        rhs(2*ref-1) = (n == 1)*(Vesp(ref) - 1);
        rhs(2*ref)   = 0;


        %% d) Lado Direito da Linha de Controle de Frequencia
        % Perdas nos shunts e na rede; as da rede usam a serie de corrente
        % I[m] = Ysr*V[m]: sum_p V[p].'*conj(I[n-p])
        rc = (n == 1)*sum(P) - real(sum(conj(Ysh).*sum(V(:,1:n).*conj(V(:,n:-1:1)), 2)));

        if n > 1
            rc = rc - real(sum(sum(V(:,2:n).*conj(Ic(:,n:-1:2)))));
        end
        rhs(nu) = rc - sum(cre(PV).*Vre(PV));


        %% e) Solucao da Ordem n
        x  = Qc*(U\(L\(Pm*rhs)));
        Vn = x(1:2:2*N) + 1j*x(2:2:2*N);
        Vn(PV) = Vre(PV) + 1j*x(2*PV);

        V(:,n+1)   = Vn;
        Ic(:,n+1)  = Ysr*Vn;
        Qh(PV,n+1) = x(2*PV-1);
        lam(n+1)   = x(nu);


        %% f) Inverso da Tensao
        % W[n] = -sum_{m=0}^{n-1} W[m] V[n-m]
        W(:,n+1) = -sum(W(:,1:n).*V(:,n+1:-1:2), 2);


        %% g) Nova Antidiagonal da Tabela Epsilon
        sS = sS + [V(:,n+1); lam(n+1)];
        [wS, eS] = Wynn_Passo(wS, sS);



        %% ________________________________________________________________
        %%
        %%                    5. CRITERIO DE PARADA
        %% ________________________________________________________________

        %% a) Ordens Iniciais
        if n < nmin
            hist(end+1) = Inf;
            continue
        end


        %% b) Estimativas em gamma = 1
        Vp   = eS(1:N);
        lamp = real(eS(N+1));

        if isempty(Vp_ant)
            Vp_ant   = Vp;
            lamp_ant = lamp;
            hist(end+1) = Inf;
            continue
        end


        %% c) Diferenca entre Estimativas Consecutivas
        err = max(max(abs(Vp - Vp_ant)), abs(lamp - lamp_ant));
        hist(end+1) = err;
        Vp_ant   = Vp;
        lamp_ant = lamp;

        if err < xi
            conv = true;
            break
        end


        %% d) Melhor Estimativa e Estagnacao
        if err < errbest
            errbest = err;
            nbest   = n;
            Vbest   = Vp;
            lambest = lamp;
        elseif n - nbest >= nestag
            estag = true;
            break
        end
    end



    %% ____________________________________________________________________
    %%
    %%                            6. SAIDAS
    %% ____________________________________________________________________

    %% a) Solucao
    % Convergiu: ultima estimativa; senao, a de menor erro
    if conv || nbest == 0
        V1   = Vp;
        lam1 = lamp;
        info.n = n;
    else
        V1   = Vbest;
        lam1 = lambest;
        info.n = nbest;
    end


    %% b) Informacoes
    info.V        = V(:,1:n+1);
    info.lam      = lam(1:n+1);
    info.Qh       = Qh(:,1:n+1);
    info.nrodadas = n;
    info.conv     = conv;
    info.estag    = estag;
    info.hist     = hist;
end
