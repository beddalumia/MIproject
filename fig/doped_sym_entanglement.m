set(0,'defaulttextinterpreter','latex')
%   addpath ../lib/qcmp-lab/

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/DOPED_MOTT/Uloc2.3')

[mu,mudir] = QcmP.post.get_list('xmu');

%s1 = QcmP.post.eentropy_line('xmu','1sites');
%s2 = QcmP.post.eentropy_line('xmu','1sites');
%s12 = QcmP.post.eentropy_line('xmu','2sites');
QcmP.post.observables_line('xmu','site001');
dens = load('dens__site001.txt');
%docc = load('docc_1_site001.txt');

RDM1 = cell(size(mu));
RDM2 = cell(size(mu));
RDM12 = cell(size(mu));
s1 = zeros(size(mu));
s2 = zeros(size(mu));
s12 = zeros(size(mu));
z = zeros(size(mu));
for i = 1:length(mu)
   cd(mudir(i))
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
assert(all(abs(nSSR - En)<10^-5))
% P-SSR == N-SSR + Ed?
assert(all(abs(pSSR - nSSR - Ed)<10^-5))

writematrix(En)
writematrix(Ed)
writematrix(mu)

cd(HERE)


%stop

%% Actual graphics

QcmP.plot.import_colorlab

tiledlayout(1,3)

%nexttile

% Upper bound on REE 
%plot(1-dens,MI,':s','LineWidth',0.5,'Color',[0.49400,0.18400,0.55600])
%hold on
%fm = fill([mu;flipud(mu)],[MI;flipud(pSSR)],[1,0.698039215686274,0.815686274509804],...
%    'EdgeColor','none'); fm.FaceAlpha=0.5;
% Lower bound on REE (SSR entanglement)
%plot(1-dens,nSSR,'.:','LineWidth',0.5,'Color',str2rgb('Neon Blue'))
%plot(1-dens,pSSR,':o','LineWidth',0.5,'Color',str2rgb('Hot Pink'))
% Scatter
%plot(1-dens,MI,'s','LineWidth',0.5,'Color',[0.49400,0.18400,0.55600])
%plot(1-dens,nSSR,'.','LineWidth',0.5,'Color',str2rgb('Neon Blue'))
%plot(1-dens,pSSR,'o','LineWidth',0.5,'Color',str2rgb('Hot Pink'))
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

%plot(1-dens,dens,':s','LineWidth',0.5,'Color',hex2rgb('#211957'));
%hold on
%plot(1-dens,docc,':d','LineWidth',0.5,'Color',hex2rgb('#BE527C'));
%plot(1-dens,z,':s','LineWidth',0.5,'Color',hex2rgb('#211957'))

nexttile

% Symmetry-resolved entanglement
plot(1-dens,En,':+','LineWidth',0.5,'Color',str2rgb('sky')); hold on
plot(1-dens,Ed,':x','LineWidth',0.5,'Color',str2rgb('strawberry'))

% Axes
xlabel("$\delta = 1-n$")
ylabel("[bit]")
%ylim([0,1.05]);
legend(["$E^{\uparrow\downarrow}_{\langle ij \rangle}$",...
        "$E^\mathrm{hd}_{\langle ij \rangle}$"],...
    "Interpreter",'latex','Location','northeast','Orientation','horizontal')
legend('boxoff')
xlim([-0.02,0.4])
%rmpath ../lib/qcmp-lab/

nexttile

% Symmetry-resolved negativity
plot(1-dens,log2(2*N0+1),':+','LineWidth',0.5,'Color',str2rgb('cornflower')); hold on
plot(1-dens,log2(2*N1+1),':o','LineWidth',0.5,'Color',str2rgb('amber'))
plot(1-dens,log2(2*N2+1),':x','LineWidth',0.5,'Color',str2rgb('rose pink'))

% Axes
xlabel("$\delta = 1-n$")
ylabel("[bit]")
%ylim([0,1.05]);
legend(["$\mathcal{N}^\mathrm{F0}_{\langle ij \rangle}$",...
        "$\mathcal{N}^\mathrm{F1}_{\langle ij \rangle}$",...
        "$\mathcal{N}^\mathrm{F2}_{\langle ij \rangle}$"],...
    "Interpreter",'latex','Location','east','Orientation','horizontal')
legend('boxoff')
xlim([-0.02,0.4])

%% Export to TikZ
addpath([HERE,'/../lib/m2tex/src']);
%matlab2tikz('doped_entanglement.tex','strict',true,...
%    'width','0.9\textwidth','height','1.5\textwidth')
rmpath([HERE,'/../lib/m2tex/src']);

nexttile

% Full entanglement
plot(1-dens,MI,':s','LineWidth',0.5,'Color',str2rgb('matlab4')); hold on
plot(1-dens,logN,':d','LineWidth',0.5,'Color',str2rgb('lilac'));

% Axes
xlabel("$\delta = 1-n$")
ylabel("[bit]")
%ylim([0,1.05]);
legend(["$I_{\langle ij \rangle}$",...
        "$\mathcal{N}^\mathrm{F}_{\langle ij \rangle}$"],...
    "Interpreter",'latex','Location','northeast','Orientation','horizontal')
legend('boxoff')
xlim([-0.02,0.4])
ylim([0,1])
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

    [mold,UDIR] = QcmP.post.get_list('xmu');
 
    E = zeros(size(mold));
    N = zeros(size(mold));
 
    for i = 1:length(mold)
       cd(UDIR(i))
       [~,E(i),~,N(i)] = negativity(RDMs{i});
       cd('..')
    end
 
 end

function [pE,nE] = build_SSRs(RDMs)

   [mold,UDIR] = QcmP.post.get_list('xmu');

   pE = zeros(size(mold));
   nE = zeros(size(mold));

   for i = 1:length(mold)
      cd(UDIR(i))
      [pE(i),nE(i)] = build_SSR(RDMs{i});
      cd('..')
   end

end

function [x,y] = symmetry_resolved(RDMs)

   [mold,UDIR] = QcmP.post.get_list('xmu');
   x = zeros(size(mold));
   y = zeros(size(mold));

   for i = 1:length(mold)
      cd(UDIR(i))
      [x(i),y(i)] = build_symE(RDMs{i});
      cd('..')
   end

end

function [x,y,z] = imbalance_negativities(RDMs)

   [mold,UDIR] = QcmP.post.get_list('xmu');
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
 
