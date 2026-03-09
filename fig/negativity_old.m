%% Compute negativity for a two-orbital RDM in usual EDIPACK ordering
function [E,N] = negativity(RDMij)
    % Preprocess the RDM
    TiRDM = partial_transpose(RDMij);
    % Get the eigenvalues
    p = eig(TiRDM);
    % Test for negativity
    N = -sum(p(p<0));
    % Measure the entanglement
    E = log2(2*N+1);
 end
 %% Partial transpose on the "left" single site
 function Ti_RDM_ij = partial_transpose(RDM_ij)
    % HARDCODED MAGIC NUMBERS (we have a simple dimer)
    Nlat = 2;
    Norb = 1;
    % Let's assert the RDM has the right dimensions
    Nrdm = size(RDM_ij,1);
    Nlso = Nlat*Norb*2;
    Nrdm = 2^Nlso;
    assert(all(size(RDM_ij)==[Nrdm,Nrdm]))
    % Init the partial-transposed RDM to size(RDM)
    Ti_RDM_ij = zeros(Nrdm,Nrdm);
    % Rotate the RDM_ij to the proper block form
    for i = 1:Nrdm
       ket = fliplr(dec2bin(i-1,Nlso));
       kup = ket(1:Nlat*Norb);
       kdw = ket(Nlat*Norb+1:end);
       for ik = 1:Nlso
          if(mod(ik,2)==1)
                ket(ik) = kup(round((ik+1)/2));
          else
                ket(ik) = kdw(round(ik/2));
          end
       end
       newI = bin2dec(fliplr(ket))+1;
       for j = 1:Nrdm
             bra = fliplr(dec2bin(j-1,Nlso));
             bup = bra(1:Nlat*Norb);
             bdw = bra(Nlat*Norb+1:end);
             for ik = 1:Nlso
                if(mod(ik,2)==1)
                   bra(ik) = bup(round((ik+1)/2));
                else
                   bra(ik) = bdw(round(ik/2));
                end
             end
             newJ = bin2dec(fliplr(bra))+1;
             % (i,j) ---> (newI,newJ)
             Ti_RDM_ij(newI,newJ) = RDM_ij(i,j);
       end
    end
    % Finally, the partial transpose on each "left site" block
    N_1site = 4;
    for i = 1:Nrdm/N_1site
       for j = 1:Nrdm/N_1site
          block = Ti_RDM_ij(1+(i-1)*N_1site:i*N_1site,1+(j-1)*N_1site:j*N_1site);
          Ti_RDM_ij(1+(i-1)*N_1site:i*N_1site,1+(j-1)*N_1site:j*N_1site) = block';
       end
    end
 end