set(0,'defaulttextinterpreter','latex')

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/doped/Uloc2.3')

[mu,mudir] = QcmP.post.get_list('xmu')

si = QcmP.post.eentropy_line('xmu','1sites');
sj = QcmP.post.eentropy_line('xmu','1sites');
sij = QcmP.post.eentropy_line('xmu','2sites');

RDMi = cell(size(mu)); 
RDMj = cell(size(mu)); 
RDMij = cell(size(mu));
si = zeros(size(mu));
sj = zeros(size(mu));
sij = zeros(size(mu));
for i = 1:length(mu)
   cd(mudir(i))
   RDMi{i} = load('reduced_density_matrix_1sites.dat');
   RDMj{i} = load('reduced_density_matrix_1sites.dat');
   RDMij{i} = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
   si(i) = vonNeumann(RDMi{i});
   sj(i) = vonNeumann(RDMj{i});
   sij(i) = vonNeumann(RDMij{i});
   cd('..')
end

MI =  si + sj - sij;

print_basis

[pSSR,nSSR] = build_SSRs(RDMij);

[logN,N] = get_negativities(RDMij);

cd(HERE)

%% Actual graphics

QcmP.plot.import_colorlab

tiledlayout(1,2)

nexttile

% Upper bound on REE 
plot(mu,MI,':s','LineWidth',0.5,'Color',[0.49400,0.18400,0.55600])
hold on
fm = fill([mu;flipud(mu)],[MI;flipud(pSSR)],[1,0.698039215686274,0.815686274509804],...
    'EdgeColor','none'); fm.FaceAlpha=0.5;
% Lower bound on REE (SSR entanglement)
plot(mu,nSSR,'.:','LineWidth',0.5,'Color',str2rgb('Neon Blue'))
plot(mu,pSSR,':o','LineWidth',0.5,'Color',str2rgb('Hot Pink'))
% Scatter
plot(mu,MI,'s','LineWidth',0.5,'Color',[0.49400,0.18400,0.55600])
plot(mu,nSSR,'.','LineWidth',0.5,'Color',str2rgb('Neon Blue'))
plot(mu,pSSR,'o','LineWidth',0.5,'Color',str2rgb('Hot Pink'))
% Axes
xlabel("$\mu/D$")
ylabel("[bit]")
ylim([0,1.05]);
legend(["$I_{\langle ij \rangle}$",...
        "$E_{\langle ij \rangle}$",...
        "$E_{\langle ij \rangle}^\mathrm{N-SSR}$",...
        "$E_{\langle ij \rangle}^\mathrm{P-SSR}$",],...
    "Interpreter",'latex','Location','northwest')
legend('boxoff')

nexttile

% Upper bound on REE 
% Upper bound on distillable entanglement
plot(mu,logN,':d','LineWidth',0.5,'Color',[[0.5882352941176471,0.1411764705882353,0.396078431372549]])
hold on
% Filling area above zero (distillable entanglement...)
fn = fill([mu;flipud(mu)],[logN;zeros(size(logN))],[1.00000,0.72549,0.00000],...
    'EdgeColor','none'); fn.FaceAlpha=0.3
% Scatter
plot(mu,logN,'d','LineWidth',0.5,'Color',[[0.5882352941176471,0.1411764705882353,0.396078431372549]])
% Axes
xlabel("$\mu/D$")
ylabel("[bit]")
ylim([0,1.05]);
legend(["$N_{\langle ij \rangle}$",...
        "$E^{\mathrm{D}}_{\langle ij \rangle}$"],...
    "Interpreter",'latex','Location','northwest')
legend('boxoff')

%% Export to TikZ
addpath([HERE,'/../lib/m2tex/src']);
matlab2tikz('doped_entanglement.tex','strict',true,...
    'width','0.9\textwidth','height','1.5\textwidth')
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
       [E(i),N(i)] = negativity(RDMs{i});
       cd('..')
    end
 
 end

function [pE,nE] = build_SSRs(RDMs)

   [mold,UDIR] = QcmP.post.get_list('xmu');

   nRDMs = RDMs; pRDMs = RDMs;
   pE = zeros(size(mold));
   nE = zeros(size(mold));

   for i = 1:length(mold)
      cd(UDIR(i))
      size(RDMs{i})
      [pE(i),nE(i)] = build_SSR(RDMs{i});
      %[pRDMs{i},nRDMs{i}] = filter_RDM(RDMs{i});
      %nE(i) = build_nSSR(nRDMs{i},1);
      %pE(i) = build_pSSR(pRDMs{i},1);
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
 