set(0,'defaulttextinterpreter','latex')

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

%cd('../../Data/CDMFT/Carlos/EEs/')

cd(HERE)

cd('../../Data/CDMFT/4sites2replicas/')

Uloc = load('U_list.txt');
S1 = postDMFT.eentropy_line('1sites');
S2 = postDMFT.eentropy_line('2sites');

Ifull =  2.*S1 - S2;

print_basis

[pSSR_I,nSSR_I] = build_SSRs_I(S1,S2);

[pSSR_S1,nSSR_S1] = build_SSRs_S(S1);

cd(HERE)

%% Actual graphics

plotDMFT.import_colorlab

figure("Name",'Plaquette Reducing Factors')

subplot(2,2,1,'align')
plot(Uloc,S1,'LineWidth',1.5,'Color',str2rgb('matlab4'))
xlim([0,3])
ylim([1,2])
ylabel("Units of $\log(2)$")
xlabel(""); xticklabels([])

subplot(2,2,2,'align')
plot(Uloc,Ifull,'LineWidth',1.5,'Color',str2rgb('matlab1'))
%xlim([0,3])
ylim([0.48,1.05])
ylabel("Units of $\log(2)$")
xlabel(""); xticklabels([])

subplot(2,2,3,'align')
plot(Uloc,S1./pSSR_S1)
hold on
plot(Uloc,S1./nSSR_S1)
xlim([0,3])
ylim([0,11])
xlabel("$U/D$")

subplot(2,2,4,'align')
plot(Uloc,Ifull./pSSR_I)
hold on
plot(Uloc,Ifull./nSSR_I)
%xlim([0,3])
%ylim([0,11])
xlabel("$U/D$")



%% Utilities

function [pEE,nEE] = build_SSRs_S(fullEE)

   [pmold,UDIR] = postDMFT.get_list('U');

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

   [mold,UDIR] = postDMFT.get_list('U');

   RDMs = cell(size(mold)); 
   nRDMs = RDMs; pRDMs = RDMs;
   pS2 = zeros(size(mold));
   nS2 = zeros(size(mold));

   for i = 1:length(mold)
      cd(UDIR(i))
      RDMs{i} = postDMFT.get_Hloc('reduced_density_matrix_2sites.dat');
      [pRDMs{i},nRDMs{i}] = filter_RDM(RDMs{i});
      pS2(i) = vonNeumann(pRDMs{i});
      nS2(i) = vonNeumann(nRDMs{i});
      assert(abs(S2(i)-vonNeumann(RDMs{i}))<1e-12)
      cd('..')
   end

   pMI = 2.*S1 - pS2;
   nMI = 2.*S1 - nS2;

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
   E = -sum(real(p.*log2(p)));
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