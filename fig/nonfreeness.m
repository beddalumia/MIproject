set(0,'defaulttextinterpreter','latex')

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

% cd('../../Data/CDMFT/Carlos/LATEST/')
% 
% temp = readtable('SingleSite2.dat');
% Uloc = temp.Var1;
% S1x2 = temp.Var2;
% S2x2 = temp.Var3;
% S3x2 = temp.Var4;
% S4x2 = temp.Var5;
% %temp = readtable('Dimer.dat');
% %temp = readtable('Plaquette.dat');
% 
% cd(HERE)

cd('../../Data/CDMFT/2sites5replicas/')

U2 = load('U_list.txt');
S2_1s = QcmP.post.eentropy_line('U','1sites');
S2_2s = QcmP.post.eentropy_line('U','2sites');

N1_dimer = zeros(size(U2));
N2_dimer = zeros(size(U2));
for i=1:length(U2)
    cd(sprintf("U=%f",U2(i)))
    RDM1 = load('reduced_density_matrix_rank4.dat');
    RDM2 = QcmP.post.get_Hloc('reduced_density_matrix_rank16.dat');
    cd ..
    N1_dimer(i) = get_nonfreeness(RDM1);
    N2_dimer(i) = get_nonfreeness(RDM2);
end
assert(all(abs(2-S2_1s-N1_dimer)<1e-3))

cd(HERE)

cd('../../Data/CDMFT/4sites2replicas/')

U4 = load('U_list.txt');
S4_1s = QcmP.post.eentropy_line('U','1sites');
S4_2s = QcmP.post.eentropy_line('U','2sites');
S4_3s = QcmP.post.eentropy_line('U','3sites');
S4_4s = QcmP.post.eentropy_line('U','4sites');

N1_plaquette = zeros(size(U4));
N2_plaquette = zeros(size(U4));
N3_plaquette = zeros(size(U4));
N4_plaquette = zeros(size(U4));
for i=1:length(U4)
    cd(sprintf("U=%f",U4(i)))
    RDM1 = load('reduced_density_matrix_1sites.dat');
    RDM2 = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
    RDM3 = QcmP.post.get_Hloc('reduced_density_matrix_3sites.dat');
    RDM4 = QcmP.post.get_Hloc('reduced_density_matrix_4sites.dat');
    cd ..
    N1_plaquette(i) = get_nonfreeness(RDM1);
    N2_plaquette(i) = get_nonfreeness(RDM2);
    %N3_plaquette(i) = get_nonfreeness(RDM3);
    %N4_plaquette(i) = get_nonfreeness(RDM4);
end
assert(all(abs(2-S4_1s-N1_plaquette)<1e-3))

cd(HERE)

%% Actual graphics

QcmP.plot.import_colorlab

figure("Name",'Nonfreness scaling')
plot(U2*4,2-S2_1s,'-+','LineWidth',1.5,'Color',str2rgb('lilac'))
hold on
plot(U4*4,2-S4_1s,'-*','LineWidth',1.5,'Color',str2rgb('matlab4'))
plot(U2*4,N2_dimer,'-+','LineWidth',1.5,'Color',str2rgb('red'))
plot(U4*4,N2_plaquette,'-*','LineWidth',1.5,'Color',str2rgb('light green'))
plot(U4*4,N3_plaquette,'-*','LineWidth',1.5,'Color',str2rgb('matlab3'))
plot(U4*4,N4_plaquette,'-*','LineWidth',1.5,'Color',str2rgb('red'))
% ylim([1,2])
% plot(Uloc,S1x2,'o','LineWidth',1.5,'MarkerSize',10,'Color',str2rgb('pyplot2'))
% plot(Uloc,S2x2,'x','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('00A4EF'))
% plot(Uloc,S3x2,'s','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('FFB900'))
% plot(Uloc,S4x2,'^','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('7FBA00'))
xline(0.5,':')
xlabel("$U/t$")
ylabel("Units of $\log(2)$")
title('Nonfreness scaling')
set(gca,'FontSize',15)
legend(["1-site, from $1\times2$, $N_\mathrm{bath}=10$ (CDMFT/ED)",...
        "1-site, from $2\times2$, $N_\mathrm{bath}=08$ (CDMFT/ED)",...
        "2-site, from $1\times2$, $N_\mathrm{bath}=10$ (CDMFT/ED)",...
        "2-site, from $2\times2$, $N_\mathrm{bath}=08$ (CDMFT/ED)",...
        "3-site, from $2\times2$, $N_\mathrm{bath}=08$ (CDMFT/ED)",...
        "4-site, from $2\times2$, $N_\mathrm{bath}=08$ (CDMFT/ED)",...
        "$U=2t$"],...
        'Interpreter','latex'); legend('boxoff')

%% Export to TikZ
% addpath([HERE,'/../lib/m2tex/src']);
% matlab2tikz('local_scaling.tex','strict',true,'noSize',true)
% rmpath([HERE,'/../lib/m2tex/src']);


figure("Name",'Nonfreness per site')
plot(U2*4,2-S2_1s,'-+','LineWidth',1.5,'Color',str2rgb('lilac'))
hold on
plot(U4*4,2-S4_1s,'-*','LineWidth',1.5,'Color',str2rgb('matlab4'))
plot(U2*4,N2_dimer/2,'-+','LineWidth',1.5,'Color',str2rgb('red'))
plot(U4*4,N2_plaquette/2,'-*','LineWidth',1.5,'Color',str2rgb('light green'))
plot(U4*4,N3_plaquette/3,'-*','LineWidth',1.5,'Color',str2rgb('matlab3'))
plot(U4*4,N4_plaquette/4,'-*','LineWidth',1.5,'Color',str2rgb('red'))
% ylim([1,2])
% plot(Uloc,S1x2,'o','LineWidth',1.5,'MarkerSize',10,'Color',str2rgb('pyplot2'))
% plot(Uloc,S2x2,'x','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('00A4EF'))
% plot(Uloc,S3x2,'s','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('FFB900'))
% plot(Uloc,S4x2,'^','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('7FBA00'))
xline(0.5,':')
xlabel("$U/t$")
ylabel("Units of $\log(2)$")
title('Nonfreness per site')
set(gca,'FontSize',15)
legend(["1-site, from $1\times2$, $N_\mathrm{bath}=10$ (CDMFT/ED)",...
        "1-site, from $2\times2$, $N_\mathrm{bath}=08$ (CDMFT/ED)",...
        "2-site, from $1\times2$, $N_\mathrm{bath}=10$ (CDMFT/ED)",...
        "2-site, from $2\times2$, $N_\mathrm{bath}=08$ (CDMFT/ED)",...
        "3-site, from $2\times2$, $N_\mathrm{bath}=08$ (CDMFT/ED)",...
        "4-site, from $2\times2$, $N_\mathrm{bath}=08$ (CDMFT/ED)",...
        "$U=2t$"],...
        'Interpreter','latex'); legend('boxoff')

%% contains

function N = get_nonfreeness(RDM)
% Compute Nonfreeness as the relative entropy between the given state RDM
% and the closest free state (Gaussian density matrix). 
%
%    >> N = S(❬c†_i c_j❭) + S(1-❬c†_i c_j❭) - S(RDM)
%
% Where ❬c†_i c_j❭ is evaluated on the given RDM.

    % Get number of single-body modes
    N = size(RDM,1); assert(N==size(RDM,2),"Not a square matrix!");
    nmodes = round(log(N)/log(4)); % N = 4^nmodes

    % Pre-compute all the one-body Slater-Condon excitations
    sc_matrix = SlaterCondon(nmodes);

    % Compute the /spinful/ one-body density matrix ❬c†_i c_j❭
    obdm = zeros(2*nmodes,2*nmodes);
    for n = 1:2*nmodes
        for m = 1:2*nmodes
            obdm(n,m) = trace(RDM * squeeze(sc_matrix(n,m,:,:)));
        end
    end
    
    % Get the entropies, hence the nonfreeness
    f = eig(obdm); p = eig(RDM);
    N = Shannon(f) + Shannon(1-f) - Shannon(p);

end

function S = Shannon(p)
% Input: a probability vector
% Output: its Shannon entropy
q = p(p>0);
S = -sum(q.*log2(q));
end

function smatrix = SlaterCondon(nmodes)
    % SLATERCONDON : Implementation of Slater-Condon rules for fermions
    %                It pre-computes all the ❬istate| cdg_is c_js |jstate❭
    %                matrix elements, storing them in a 4D array.
    %                This can be used to build many-body representations of
    %                one-body operators. NB: it assumes Nup,Ndw conservation.
    %
    %  >> smatrix = SlaterCondon(nmodes :: number of single-fermion modes)
    %     smatrix :: 4D array [2*nmodes,2*nmodes,4^nmodes,4^nmodes]
    %     
    % ! This can be made way faster by implementing Slater-Condon
    %   rules in an efficient machine-tuned way, as discussed in
    %   https://arxiv.org/abs/1311.6244 (hal-01539072)
    %
    N = 4^nmodes;
    smatrix = zeros(2*nmodes,2*nmodes,N,N);
    for istate = 0:1:N-1 % nmode-orbital states
        for jstate = 0:1:N-1 % nmode-orbital states
            % ∑_s ∑_ij ❬istate| cdg_is c_js |jstate❭
            for ispin=1:2
                for imode=1:nmodes
                    for jmode=1:nmodes
                        % Apply cdg_is to ❬istate|
                        ibra = build_ket(istate,nmodes);
                        if ibra(imode+(ispin-1)*nmodes)==0
                            continue
                        end
                        ibra(imode+(ispin-1)*nmodes) = ibra(imode+(ispin-1)*nmodes)-1;
                        isign = (-1)^(sum(ibra(1:imode+(ispin-1)*nmodes)));
                        % Apply c_js to |jstate❭
                        jket = build_ket(jstate,nmodes);
                        if jket(jmode+(ispin-1)*nmodes)==0
                            continue
                        end
                        jket(jmode+(ispin-1)*nmodes) = jket(jmode+(ispin-1)*nmodes)-1;
                        jsign = (-1)^(sum(jket(1:jmode+(ispin-1)*nmodes)));
                        % Overlap ❬istate| cdg_is c_js |jstate❭
                        if isequal(ibra,jket)
                            smatrix(imode+(ispin-1)*nmodes,jmode+(ispin-1)*nmodes,istate+1,jstate+1) = isign * jsign;
                        end
                    end
                end
            end
        end
    end
end
% function vec = build_ket(state,nmodes)
%   % BUILD_KET : Puts together a bit representation of a Slater determinant
%   %  
%   %  >> [vec,ket] = build_ket(state)
%   %
%   %  state :: integer representation of a basis state (bits are occupation numbers)
%   %
%   % The ordering of the bits follows the convention of EDIpack: |up❭ ⊗ |dw❭
%   %
%   vec = zeros(4^nmodes,1);
%   for imode = 1:nmodes
%       for ispin = 1:2
%           index = imode + (ispin-1)*nmodes;
%           vec(index) = bitget(state,index);
%       end
%   end
% end
 % FROM ED_SETUP:
    % |imp_up>|bath_up> * |imp_dw>|bath_dw>        <- 2*Nlat*Norb bits
    % |imp_sigma> = | (1…Norb)_1...(1…Norb)_Nlat > <--- Nlat*Norb bits
    % lso indices are: io = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
    function [vec,ket] = build_ket(state,Nlat)
      % BUILD_KET : Puts together a bit representation of a Slater determinant
      %             (and, if asked, a pretty string to print it to screen...)
      %  
      %  >> [vec,ket] = build_ket(state)
      %
      %  state :: integer representation of a basis state (bits are occupation numbers)
      %
      %  This depends entirely on the Fock basis conventions choosen in all
      %  ED-based solvers from QcmPlab (and many other codes, actually)
      %
      %  Copyright 2022 Gabriele Bellomia
      %
      %Nlat = 2;
      Norb = 1;
      for ilat = 1:Nlat
          for ispin = 1:2
              shift = (ilat-1)*Norb + (ispin-1)*Norb*Nlat;
              index = shift+(1:Norb);
              vec(index) = bitget(state,index);
          end
      end
      if nargout>1
          kup = num2str(vec(1:Norb*Nlat));
          kdw = num2str(vec(Norb*Nlat+1:end));
          ket = ['| ',strrep(kup,'1','↑'),' 〉⊗ ',...
              '| ',strrep(kdw,'1','↓'), ' 〉'];
          ket = strrep(ket,'0','•');
      end
  end
  function [nRDM,pRDM] = SSR_filter(RDM)
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

