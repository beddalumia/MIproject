set(0,'defaulttextinterpreter','latex')

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/4sites2replicas/')

Uloc = load('U_list.txt');
S1 = QcmP.post.eentropy_line('U','1sites');
S2 = QcmP.post.eentropy_line('U','2sites');

Ifull =  2.*S1 - S2;

print_basis

[pSSR_I,nSSR_I] = build_SSRs_I(S1,S2);

[pSSR_S1,nSSR_S1] = build_SSRs_S(S1);

for i = 1:length(Uloc)
    cd(sprintf("U=%f",Uloc(i)))
    rdm = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
    neg = negativity(rdm,false);
    Nfull(i) = neg;
    cd('..')
end

[pSSR_N,nSSR_N] = build_SSRs_N(Nfull);

cd(HERE)

%% Actual graphics

QcmP.plot.import_colorlab

figure("Name",'Plaquette Reducing Factors')

tiledlayout(2,3,'TileSpacing','compact')

nexttile
plot(Uloc,S1,'LineWidth',1.5,'Color',str2rgb('matlab4'))
hold on
plot(Uloc,nSSR_S1,'-','LineWidth',1.5,'Color',str2rgb('Neon Blue'))
plot(Uloc,pSSR_S1,'--','LineWidth',1.5,'Color',str2rgb('Hot Pink'))
xlim([0,8])
%ylim([1,2])
legend(["$s_i$","$E_i^\text{N-SSR}$","$E_i^\text{P-SSR}$"],'Location','southwest','Interpreter','latex')
legend('boxoff')
ylabel("Units of $\log(2)$")
xlabel(""); %xticklabels([])

nexttile
plot(Uloc,Ifull,'LineWidth',1.5,'Color',str2rgb('pyplot3'))
hold on
plot(Uloc,nSSR_I,'-','LineWidth',1.5,'Color',str2rgb('Neon Blue'))
plot(Uloc,pSSR_I,'--','LineWidth',1.5,'Color',str2rgb('Hot Pink'))
xlim([0,8])
ylim([0,1.1])
legend(["$I_{\langle ij \rangle}$","$I_{\langle ij \rangle}^\text{N-SSR}$","$I_{\langle ij \rangle}^\text{P-SSR}$"],'Location','southeast','Interpreter','latex')
legend('boxoff')
xlabel(""); %xticklabels([])

nexttile
plot(Uloc,Nfull,'LineWidth',1.5,'Color',[0.5882352941176471,0.1411764705882353,0.396078431372549])
hold on
plot(Uloc,nSSR_N,'-','LineWidth',1.5,'Color',str2rgb('Neon Blue'))
plot(Uloc,pSSR_N,'--','LineWidth',1.5,'Color',str2rgb('Hot Pink'))
xlim([0,8])
ylim([0.0,1.1])
legend(["$\mathcal{N}_{\langle ij \rangle}$","$\mathcal{N}_{\langle ij \rangle}^\text{N-SSR}$","$\mathcal{N}_{\langle ij \rangle}^\text{P-SSR}$"],'Location','southeast','Interpreter','latex')
legend('boxoff')
xlabel(""); %xticklabels([])

nexttile
plot(Uloc,S1./nSSR_S1,'-','LineWidth',1.5,'Color',str2rgb('Neon Blue'))
hold on
plot(Uloc,S1./pSSR_S1,'--','LineWidth',1.5,'Color',str2rgb('Hot Pink'))
xlim([0,8])
ylim([1,50])
xlabel("$U/D$")
ylabel("Superselection Factor")
legend(["N-SSR","P-SSR"],'Location','northwest')
legend('boxoff')

nexttile
plot(Uloc,Ifull./nSSR_I,'-','LineWidth',1.5,'Color',str2rgb('Neon Blue'))
hold on
plot(Uloc,Ifull./pSSR_I,'--','LineWidth',1.5,'Color',str2rgb('Hot Pink'))
xlim([0,8])
ylim([1,7])
xlabel("$U/D$")
%ylabel("Superselection Divisor")
legend(["N-SSR","P-SSR"],'Location','northeast')
legend('boxoff')

nexttile
plot(Uloc,Nfull'./nSSR_N,'-','LineWidth',1.5,'Color',str2rgb('Neon Blue'))
hold on
plot(Uloc,Nfull'./pSSR_N,'--','LineWidth',1.5,'Color',str2rgb('Hot Pink'))
xlim([0,8])
ylim([1,16])
xlabel("$U/D$")
%ylabel("Superselection Divisor")
legend(["N-SSR","P-SSR"],'Location','northeast')
legend('boxoff')

%% Export to TikZ
addpath([HERE,'/../lib/m2tex/src']);
matlab2tikz('superselection.tex','strict',true,...
    'width','0.9\textwidth','height','0.6\textwidth')
rmpath([HERE,'/../lib/m2tex/src']);

%% Utilities

function [pEE,nEE] = build_SSRs_S(fullEE)

   [pmold,UDIR] = QcmP.post.get_list('U');

   p1 = zeros(size(pmold));
   p2 = zeros(size(pmold));
   p3 = zeros(size(pmold));
   p4 = zeros(size(pmold));

   for i = 1:length(pmold)
      cd(UDIR(i))
      ptmp = load('probabilities_1sites.dat');
      p1(i) = ptmp(1);
      p2(i) = ptmp(2);
      p3(i) = ptmp(3);
      p4(i) = ptmp(4);
      cd('..')
   end

   pEE = (p1+p4).*log2(p1+p4) + (p2+p3).*log2(p2+p3) + fullEE;
   nEE = (p2+p3).*log2(p2+p3) - p2.*log2(p2) - p3.*log2(p3);

end


function [pMI,nMI] = build_SSRs_I(S1,S2)

   [mold,UDIR] = QcmP.post.get_list('U');

   RDMs = cell(size(mold)); 
   nRDMs = RDMs; pRDMs = RDMs;
   pS2 = zeros(size(mold));
   nS2 = zeros(size(mold));

   for i = 1:length(mold)
      cd(UDIR(i))
      RDMs{i} = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
      [pRDMs{i},nRDMs{i}] = filter_RDM(RDMs{i});
      pS2(i) = vonNeumann(pRDMs{i});
      nS2(i) = vonNeumann(nRDMs{i});
      assert(abs(S2(i)-vonNeumann(RDMs{i}))<1e-12)
      cd('..')
   end

   pMI = 2.*S1 - pS2;
   nMI = 2.*S1 - nS2;

end

function [pN,nN] = build_SSRs_N(N)

   [mold,UDIR] = QcmP.post.get_list('U');

   RDMs = cell(size(mold)); 
   nRDMs = RDMs; pRDMs = RDMs;
   pN = zeros(size(mold));
   nN = zeros(size(mold));

   for i = 1:length(mold)
      cd(UDIR(i))
      RDMs{i} = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
      [pRDMs{i},nRDMs{i}] = filter_RDM(RDMs{i});
      pN(i) = negativity(pRDMs{i});
      nN(i) = negativity(nRDMs{i});
      assert(abs(N(i)-negativity(RDMs{i}))<1e-12)
      cd('..')
   end

end

function [pRDM,nRDM] = filter_RDM(RDM)
   % Trim all off-diagonals
   nRDM = diag(diag(RDM));
   % Restore the spin-flip terms
   nRDM(7,10) = RDM(7,10);
   nRDM(10,7) = RDM(10,7);
   pRDM = nRDM;
   % Restore the pair-hopping terms
   pRDM(6,11) = RDM(6,11);
   pRDM(11,6) = RDM(11,6);
end

function E = vonNeumann(RDM)
   p = eig(RDM);
   E = -sum(p.*log2(p));
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
function [E,N] = negativity(RDMij,is_asci)
    if(nargin<2)
        is_asci=false
    end
    % Preprocess the RDM
    TiRDM = partial_transpose(RDMij,is_asci);
    % Get the eigenvalues
    p = eig(TiRDM);
    % Test for negativity
    N = -sum(p(p<0));
    % Measure the entanglement
    E = log2(2*N+1);
 end
 %% Partial transpose on the "left" single site
 function Ti_RDM_ij = partial_transpose(RDM_ij,is_asci)
    % HARDCODED MAGIC NUMBERS (we have a simple dimer)
    Nlat = 2;
    Norb = 1;
    if(nargin<2)
        is_asci=false
    end
    if(is_asci)
        % HARDCODED ASCI STATES
        dets = ["0110","1001","0101","1010","0010","0001","1110","1101","1100","0011","0100","1000","0111","1011","0000","1111"];
        indx = bin2dec(dets);
        % FIX STUPID ASCI RDM
        RDM_ij = RDM_ij(indx+1,indx+1);
    end
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