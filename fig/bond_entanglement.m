set(0,'defaulttextinterpreter','latex')

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/4sites2replicas/')

Uloc = load('U_list.txt');

s1 = postDMFT.eentropy_line('1sites');
s2 = postDMFT.eentropy_line('2sites');

MI =  2.*s1 - s2;

print_basis

[pSSR,nSSR] = build_SSRs();

cd(HERE)

%% Actual graphics

plotDMFT.import_colorlab

% Local entropy as a upper bound for all 1-site correlation
plot(Uloc,s1,':','LineWidth',1.5,'Color',str2rgb('matlab4'))
hold on
% Upper bound on bond entanglement (bond mutual information)
plot(Uloc,MI,'-.','LineWidth',1.5,'Color',str2rgb('pyplot3'))
% Filling area between the bounds (entanglement interval...)
fill([Uloc;flipud(Uloc)],[nSSR;flipud(MI)],str2rgb("light khaki"),...
    'EdgeColor','none')
% Lower bound on bond entanglement (meaningful entanglement)
plot(Uloc,nSSR,'-','LineWidth',1.5,'Color',str2rgb('pea'))
plot(Uloc,pSSR,'--','LineWidth',1.5,'Color',str2rgb('peach'))

xlim([-0.01,8.01])
xlabel("$U/D$")
ylabel("Units of $\log(2)$")
%set(gca,'FontSize',15)

legend(["local entanglement entropy $s_i$",...
        "total NN correlation $s_i+s_j-s_{\langle ij \rangle}$",...
        "unrestricted NN entanglement",...
        "NN entanglement under N-SSR",...
        "NN entanglement under P-SSR",],...
    "Interpreter",'latex','Location','northeast')
legend('boxoff')

matlab2tikz('bond_entanglement.tex','strict',true,...
    'width','0.45\textwidth','height','0.6\textwidth')

%% Utilities

function [pE,nE] = build_SSRs()

   [mold,UDIR] = postDMFT.get_list('U');

   RDMs = cell(size(mold)); 
   nRDMs = RDMs; pRDMs = RDMs;
   pE = zeros(size(mold));
   nE = zeros(size(mold));

   for i = 1:length(mold)
      cd(UDIR(i))
      RDMs{i} = postDMFT.get_Hloc('reduced_density_matrix_2sites.dat');
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
   rot = [1,-1; 1,1];
   new = rot*AB*rot';
   nRDM(7,7) = new(1,1);
   nRDM(10,10) = new(2,2);
   %
   pRDM = nRDM;
   % Restore the pair-hopping terms
   pRDM(6,11) = RDM(6,11);
   pRDM(11,6) = RDM(11,6);
   % Global singlet assertion
   assert(all(abs(RDM(4,4)-RDM(13,13))<1e-12))
   % Particle-hole symmetry assertion
   assert(all(abs(RDM(1,1)-RDM(16,16))<1e-4))
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

