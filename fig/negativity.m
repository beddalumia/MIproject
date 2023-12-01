set(0,'defaulttextinterpreter','latex')

% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/4sites2replicas')

Ulist = QcmP.post.get_list('U');

Uloc=input('Which value of U? ','s');

while not(any(str2double(Uloc)==Ulist))
     disp("U must be in the following set:")
     disp(Ulist);
     Uloc=input('Which value of U? ','s');
end

cd(sprintf('U=%f',str2double(Uloc)));

RDMij = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
TiRDM = partial_transpose(RDMij);

p = eig(TiRDM);
N = -sum(p(p<0)) % Negativity
E = log2(2*N+1)  % Logarithmic Negativity

cd(HERE)

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
