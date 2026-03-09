% In generale, con N flavour e filling n il valore di aspettazione 
% <n_1 * n_2 * … * n_q> = n!/N! * (N-q)!/(n-q)!
% Dove i pedici sono i flavour
N = 4;
eucli = zeros(N,1);
bures = zeros(N,1);
mutua = zeros(N,1);
corre = zeros(N,1);
for n=1:N-1
   %dens = factorial(n)/factorial(N) * factorial(N-1)/factorial(n-1);
   %assert(abs(dens-n/N)<1e-12,"Wrong density!");
   dens = n/N;
   if n<2
      docc = 0;
   else
      docc = factorial(n)/factorial(N) * factorial(N-2)/factorial(n-2);
   end
   RDM1 = diag([dens,1-dens]);
   RDM2 = diag([1-2*dens+docc,dens-docc,dens-docc,docc]);
   PROD = kron(RDM1,RDM1); % Closest uncorrelated state
   eucli(n) = norm(RDM2-PROD,'fro'); % Frobenius (Euclidean) distance
   bures(n) = wootters(sqrt(RDM2.*PROD)); % Fidelity of classical states
   mutua(n) = 2*shannon(RDM1) - shannon(RDM2); % Mutual Information
   corre(n) = (docc-dens*dens)^2 / dens^2;
end

plot(0:N,[0;eucli],'s:','Linewidth',1.5); 
hold on
plot(0:N,[0;bures],'d:','Linewidth',1.5);
plot(0:N,[0;mutua*1e5],'*:','Linewidth',1.5);
plot(0:N,[0;corre*1e6],'^:','Linewidth',1.5);

legend(["Euclidean","Bures","Mutual Info $\times 10^5$","$|\langle n_\alpha n_\beta\rangle - \langle n_\alpha\rangle\langle n_\beta\rangle|^2 \times 10^6$"],'Interpreter','latex','Location','north')
xlabel("$n$",'Interpreter','latex')
ylabel("Correlation degree (different units)",'Interpreter','latex')
title(sprintf("SU(%d)",N))

%% definitions for diagonal density matrices

function D = wootters(A)
    try
      D = acos(sum(sum(A)));
    catch
      disp Not a matrix eheh
    end
end

function S = shannon(A)
    a = diag(A); 
    if not(isvector(a))
        error("Pass the diagonal matrix, not only the diagonal");
    end
    S = 0;
    for i = 1:length(a)
        if not(a(i)==0)
           S = S - a(i)*log2(a(i));
        end
    end
end