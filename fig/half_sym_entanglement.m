set(0,'defaulttextinterpreter','latex')
%   addpath ../lib/qcmp-lab/

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/4sites2replicas')

[u,udir] = QcmP.post.get_list('U');

%s1 = QcmP.post.eentropy_line('U','1sites');
%s2 = QcmP.post.eentropy_line('U','1sites');
%s12 = QcmP.post.eentropy_line('U','2sites');
QcmP.post.observables_line('U','site001');
dens = load('dens__site001.txt');
%docc = load('docc_1_site001.txt');

RDM1 = cell(size(u));
RDM2 = cell(size(u));
RDM12 = cell(size(u));
s1 = zeros(size(u));
s2 = zeros(size(u));
s12 = zeros(size(u));
z = zeros(size(u));
for i = 1:length(u)
   cd(udir(i))
    try
        RDM1{i} = load('reduced_density_matrix_i1l1.dat');
        RDM2{i} = load('reduced_density_matrix_i2l1.dat');
        RDM4{i} = load('reduced_density_matrix_i4l1.dat');
        RDM12{i} = QcmP.post.get_Hloc('reduced_density_matrix_i1l1_i2l1.dat');
        RDM14{i} = QcmP.post.get_Hloc('reduced_density_matrix_i1l1_i4l1.dat');
    catch
        RDM1{i} = load('reduced_density_matrix_1sites.dat');
        RDM2{i} = load('reduced_density_matrix_1sites.dat');
        RDM12{i} = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
    end
   s1(i) = vonNeumann(RDM1{i});
   s2(i) = vonNeumann(RDM2{i});
   s12(i) = vonNeumann(RDM12{i});
   s14(i) = vonNeumann(RDM14{i});
   z(i) = load('zeta_last_site001.ed');
   cd('..')
end

MI =  s1 + s2 - s12;

print_basis

[pSSR,nSSR] = build_SSRs(RDM12);
[pSSR_,nSSR_] = build_SSRs(RDM14);

[logN,N] = get_negativities(RDM12);

[En,Ed] = symmetry_resolved(RDM12);
[En_,Ed_] = symmetry_resolved(RDM14);

[N0,N1,N2] = imbalance_negativities(RDM12);
[N0_,N1_,N2_] = imbalance_negativities(RDM14);

% N-SSR == En?
assert(all(abs(nSSR - En)<10^-3))
% P-SSR == N-SSR + Ed?
assert(all(abs(pSSR - nSSR - Ed)<10^-3))

cd(HERE)

%% Actual graphics

QcmP.plot.import_colorlab

tiledlayout(1,3)

%nexttile

% Upper bound on REE 
%plot(u,MI,':s','LineWidth',0.5,'Color',[0.49400,0.18400,0.55600])
%hold on
%fm = fill([u;flipud(u)],[MI;flipud(pSSR)],[1,0.698039215686274,0.815686274509804],...
%    'EdgeColor','none'); fm.FaceAlpha=0.5;
% Lower bound on REE (SSR entanglement)
%plot(u,nSSR,'.:','LineWidth',0.5,'Color',str2rgb('Neon Blue'))
%plot(u,pSSR,':o','LineWidth',0.5,'Color',str2rgb('Hot Pink'))
% Scatter
%plot(u,MI,'s','LineWidth',0.5,'Color',[0.49400,0.18400,0.55600])
%plot(u,nSSR,'.','LineWidth',0.5,'Color',str2rgb('Neon Blue'))
%plot(u,pSSR,'o','LineWidth',0.5,'Color',str2rgb('Hot Pink'))
% Axes
%xlabel("$\delta = 1-n$")
%ylabel("[bit]")
%ylim([0,1.05]);
%legend(["$I_{\langle ij \rangle}$",...
%        "$E_{\langle ij \rangle}$",...
%        "$E_{\langle ij \rangle}^\mathrm{N-SSR}$",...
%        "$E_{\langle ij \rangle}^\mathrm{P-SSR}$",],...
%    "Interpreter",'latex','Location','northwest')
%legend('boxoff')

%plot(u,dens,':s','LineWidth',0.5,'Color',hex2rgb('#211957'));
%hold on
%plot(u,docc,':d','LineWidth',0.5,'Color',hex2rgb('#BE527C'));
%plot(u,z,':s','LineWidth',0.5,'Color',hex2rgb('#211957'))

nexttile

% Symmetry-resolved entanglement
plot(u(u>0),En(u>0),':+','LineWidth',0.5,'Color',str2rgb('sky')); hold on
plot(u(u>0),Ed(u>0),':x','LineWidth',0.5,'Color',str2rgb('strawberry'))

% Axes
xlabel("$U/D$")
ylabel("[bit]")
%ylim([0,1.05]);
legend(["$E^{\uparrow\downarrow}_{\langle ij \rangle}$",...
        "$E^\mathrm{hd}_{\langle ij \rangle}$"],...
    "Interpreter",'latex','Location','east','Orientation','horizontal')
legend('boxoff')
%rmpath ../lib/qcmp-lab/

nexttile

% Symmetry-resolved negativity
plot(u(u>0),log2(2*N0(u>0)+1),':+','LineWidth',0.5,'Color',str2rgb('cornflower')); hold on
plot(u(u>0),log2(2*N1(u>0)+1),':x','LineWidth',0.5,'Color',str2rgb('amber'))
plot(u(u>0),log2(2*N2(u>0)+1),':x','LineWidth',0.5,'Color',str2rgb('rose pink'))

% Axes
xlabel("$U/D$")
ylabel("[bit]")
%ylim([0,1.05]);
legend(["$\mathcal{N}^\mathrm{F0}_{\langle ij \rangle}$",...
        "$\mathcal{N}^\mathrm{F1}_{\langle ij \rangle}$",...
        "$\mathcal{N}^\mathrm{F2}_{\langle ij \rangle}$"],...
    "Interpreter",'latex','Location','east','Orientation','horizontal')
legend('boxoff')
%rmpath ../lib/qcmp-lab/

%% Export to TikZ
addpath([HERE,'/../lib/m2tex/src']);
%matlab2tikz('doped_entanglement.tex','strict',true,...
%    'width','0.9\textwidth','height','1.5\textwidth')
rmpath([HERE,'/../lib/m2tex/src']);

nexttile

% Full entanglement
plot(u(u>0),MI(u>0),':s','LineWidth',0.5,'Color',str2rgb('matlab4')); hold on
plot(u(u>0),logN(u>0),':d','LineWidth',0.5,'Color',str2rgb('lilac'));

% Axes
xlabel("$U/D$")
ylabel("[bit]")
%ylim([0,1.05]);
legend(["$I_{\langle ij \rangle}$",...
        "$\mathcal{N}^\mathrm{F}_{\langle ij \rangle}$"],...
    "Interpreter",'latex','Location','south','Orientation','horizontal')
legend('boxoff')
ylim([0,max(MI)])
%rmpath ../lib/qcmp-lab/

%% Export to TikZ
addpath([HERE,'/../lib/m2tex/src']);
%matlab2tikz('doped_entanglement.tex','strict',true,...
%    'width','0.9\textwidth','height','1.5\textwidth')
rmpath([HERE,'/../lib/m2tex/src']);

%% Utilities

function E = vonNeumann(RDM)
   p = eig(RDM);
   E = -sum(p.*log2(p));
end

function [E,N] = get_negativities(RDMs)

    [mold,UDIR] = QcmP.post.get_list('U');
 
    E = zeros(size(mold));
    N = zeros(size(mold));
 
    for i = 1:length(mold)
       cd(UDIR(i))
       [~,E(i),~,N(i)] = negativity(RDMs{i});
       cd('..')
    end
 
 end

function [pE,nE] = build_SSRs(RDMs)

   [mold,UDIR] = QcmP.post.get_list('U');

   pE = zeros(size(mold));
   nE = zeros(size(mold));

   for i = 1:length(mold)
      cd(UDIR(i))
      [pE(i),nE(i)] = build_SSR(RDMs{i});
      cd('..')
   end

end

function [x,y] = symmetry_resolved(RDMs)

   [mold,UDIR] = QcmP.post.get_list('U');
   x = zeros(size(mold));
   y = zeros(size(mold));

   for i = 1:length(mold)
      cd(UDIR(i))
      [x(i),y(i)] = build_symE(RDMs{i});
      cd('..')
   end

end

function [x,y,z] = imbalance_negativities(RDMs)

   [mold,UDIR] = QcmP.post.get_list('U');
   x = zeros(size(mold));
   y = zeros(size(mold));
   z = zeros(size(mold));

   for i = 1:length(mold)
      cd(UDIR(i))
      [x(i),y(i),z(i)] = sym_negativity(RDMs{i});
      cd('..')
   end

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
 
