set(0,'defaulttextinterpreter','latex')

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/4sites2replicas/')

Uloc = load('U_list.txt');

s1 = QcmP.post.eentropy_line('U','1sites');
s2 = QcmP.post.eentropy_line('U','2sites');

MI =  2.*s1 - s2;

print_basis

[pSSR,nSSR] = build_SSRs();

[logN,N] = get_negativities();

cd(HERE)

%% Actual graphics

QcmP.plot.import_colorlab

% Local entropy as a upper bound for all 1-site correlation
plot(Uloc,s1,':','LineWidth',1.5,'Color',str2rgb('matlab4'))
hold on
% Upper bound on bond entanglement (bond mutual information)
plot(Uloc,MI,'-.','LineWidth',1.5,'Color',str2rgb('pyplot3'))
% Upper bound on bond entanglement (bond logarithmic negativity)
plot(Uloc,logN,'--','LineWidth',2.5,'Color',str2rgb('goldenrod'))
% Filling area between the bounds (entanglement interval...)
fill([Uloc;flipud(Uloc)],[pSSR;flipud(logN)],str2rgb("light khaki"),...
    'EdgeColor','none')
% Lower bound on bond entanglement (meaningful entanglement)
plot(Uloc,nSSR,'-','LineWidth',1.5,'Color',str2rgb('Neon Blue'))
plot(Uloc,pSSR,'--','LineWidth',1.5,'Color',str2rgb('Hot Pink'))

xlim([0,8])
xlabel("$U/D$")
ylim([-0.05,2.05])
ylabel("Units of $\log(2)$")
%set(gca,'FontSize',15)

legend(["$s_i$",...
        "$I_{\langle ij \rangle}$",...
        "$N_{\langle ij \rangle}$",...
        "$E_{\langle ij \rangle}$",...
        "$E_{\langle ij \rangle}^\mathrm{N-SSR}$",...
        "$E_{\langle ij \rangle}^\mathrm{P-SSR}$",],...
    "Interpreter",'latex','Location','northeast')
legend('boxoff')

%% Export to TikZ
addpath([HERE,'/../lib/m2tex/src']);
matlab2tikz('bond_entanglement.tex','strict',true,...
    'width','0.3\textwidth','height','0.7\textwidth')
rmpath([HERE,'/../lib/m2tex/src']);

%% Utilities

function [E,N] = get_negativities()

    [mold,UDIR] = QcmP.post.get_list('U');
 
    RDMs = cell(size(mold)); 
    E = zeros(size(mold));
    N = zeros(size(mold));
 
    for i = 1:length(mold)
       cd(UDIR(i))
       RDMs{i} = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
       [E(i),N(i)] = negativity(RDMs{i});
       cd('..')
    end
 
 end

function [pE,nE] = build_SSRs()

   [mold,UDIR] = QcmP.post.get_list('U');

   RDMs = cell(size(mold)); 
   nRDMs = RDMs; pRDMs = RDMs;
   pE = zeros(size(mold));
   nE = zeros(size(mold));

   for i = 1:length(mold)
      cd(UDIR(i))
      RDMs{i} = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
      [pRDMs{i},nRDMs{i}] = filter_RDM(RDMs{i});
      nE(i) = build_nSSR(nRDMs{i});
      pE(i) = build_pSSR(pRDMs{i});
      cd('..')
   end

end

function E = build_nSSR(RDM)
      % N-SSR computation, assuming global singlet
      t = max(RDM(7,7),RDM(10,10));
      r = min(RDM(7,7),RDM(10,10)) + RDM(4,4) + RDM(13,13);
      E = 0;
      if r < t
         E = E + r * log2(2*r/(r+t)) + t * log2(2*t/(r+t));
      end
end

function E = build_pSSR(RDM)
   % P-SSR computation, assuming particle-hole symmetry
   t = max(RDM(7,7),RDM(10,10));
   r = min(RDM(7,7),RDM(10,10)) + RDM(4,4) + RDM(13,13);
   T = max(RDM(6,6),RDM(11,11));
   R = min(RDM(6,6),RDM(11,11)) + RDM(1,1) + RDM(16,16);
   E = 0;
   if r < t
      E = E + r * log2(2*r/(r+t)) + t * log2(2*t/(r+t));
      E = E + R * log2(2*R/(R+T)) + T * log2(2*T/(R+T));
   end
end

function [pRDM,nRDM] = filter_RDM(RDM)
   % Global singlet assertion
   assert(all(abs(RDM(4,4)-RDM(13,13))<1e-12))
   % Particle-hole symmetry assertion
   assert(all(abs(RDM(1,1)-RDM(16,16))<1e-4))
   % Trim all off-diagonals
   nRDM = diag(diag(RDM));
   % Restore the spin-flip terms
   nRDM(7,10) = RDM(7,10);
   nRDM(10,7) = RDM(10,7);
   % Rotate the spin-flip terms
   aa = nRDM(7,7);
   bb = nRDM(10,10);
   ab = nRDM(7,10);
   ba = nRDM(10,7);
   AB = [aa,ab;ba,bb];
   rot = [1,-1; 1,1]./sqrt(2);
   new = rot*AB*rot';
   nRDM(7,7) = new(1,1);
   nRDM(10,10) = new(2,2);
   pRDM = nRDM;
   % Restore the pair-hopping terms
   pRDM(6,11) = RDM(6,11);
   pRDM(11,6) = RDM(11,6);
   % Rotate the pair-hopping terms
   aa = pRDM(6,6);
   bb = pRDM(11,11);
   ab = pRDM(6,11);
   ba = pRDM(11,6);
   AB = [aa,ab;ba,bb];
   rot = [1,-1; 1,1]./sqrt(2);
   new = rot*AB*rot';
   pRDM(6,6) = new(1,1);
   pRDM(11,11) = new(2,2);
end

%% FROM ED_SETUP:
    % |imp_up>|bath_up> * |imp_dw>|bath_dw>        <- 2*Nlat*Norb bits
    % |imp_sigma> = | (1…Norb)_1...(1…Norb)_Nlat > <--- Nlat*Norb bits
    % lso indices are: io = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
    function ket = build_ket(state)
      %% BUILD_KET : Puts together a pretty label for a pure state component
      %  
      %  >> ket = build_ket(state)
      %
      %  state :: integer representation of a basis state (bits are spins)
      %
      %  This depends entirely on the Fock basis conventions choosen in all
      %  ED-based solvers from QcmPlab.
      %
      %  Copyright 2022 Gabriele Bellomia
      %
      Nlat = 2;
      Norb = 1;
      for ilat = 1:Nlat
          for ispin = 1:2
              shift = (ilat-1)*Norb + (ispin-1)*Norb*Nlat;
              index = shift+(1:Norb);
              vec(index) = bitget(state,index);
          end
      end
      kup = num2str(vec(1:Norb*Nlat));
      kdw = num2str(vec(Norb*Nlat+1:end));
      ket = ['| ',strrep(kup,'1','↑'),' 〉⊗ ',...
          '| ',strrep(kdw,'1','↓'), ' 〉'];
      ket = strrep(ket,'0','•');
  end
  %
  function print_basis()
   for state = 0:1:15
       label = build_ket(state);
       fprintf('%d\t',state+1)
       disp(label)
   end
  end

%% Compute negativity
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
 