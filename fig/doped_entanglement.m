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
      [pRDMs{i},nRDMs{i}] = filter_RDM(RDMs{i});
      nE(i) = build_nSSR(nRDMs{i},1);
      pE(i) = build_pSSR(pRDMs{i},1);
      cd('..')
   end

end

function E = build_nSSR(RDM, general)
   % N-SSR computation, assuming global singlet
   t = max(RDM(7,7),RDM(10,10));
   r = min(RDM(7,7),RDM(10,10)) + RDM(4,4) + RDM(13,13);
   E = 0;
   if r < t
      E = E + r * log2(2*r/(r+t)) + t * log2(2*t/(r+t));
   end
   if(general)
      disp('GENERAL COMPUTATION FOR N-SSR')
      p = diag(RDM);
      q = diag(RDM);
      % N-SSR computation, general formula
      p(7)  = max(RDM(7,7),RDM(10,10));
      p(10) = min(RDM(7,7),RDM(10,10));
      a1 = p(7)+p(10)+p(4)+p(13);
      b1 = a1^2-(p(4)-p(13))^2;
      c1 = (p(7)-p(10))*a1;
      d1 = (p(4)+p(13))^2*(p(7)-p(10))^2+8*p(4)*p(13)*(2*p(4)*p(13)+(p(4)+p(13))*(p(7)-p(10))+2*p(7)*p(10));
      q(7) = (b1+c1+sqrt(d1))/(4*(a1-p(10)));
      q(10) = (b1-c1-sqrt(d1))/(4*(a1-p(7)));
      q(4) = p(4) + (p(7)+p(10)-q(7)-q(10))/2;
      q(13) = p(13) + (p(7)+p(10)-q(7)-q(10))/2;
      % Final relative entropy
      E = 0; indices = [4,7,10,13];
      for i=indices
         if(p==0)
            continue
         end
         E = E + p(i)*log2(p(i)/q(i));
      end
   end
end

function E = build_pSSR(RDM,general)
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
   if(general)
   disp('GENERAL COMPUTATION FOR P-SSR')
      p = diag(RDM);
      q = diag(RDM);
      % P-SSR computation, general formula (N-SSR part)
      p(7)  = max(RDM(7,7),RDM(10,10));
      p(10) = min(RDM(7,7),RDM(10,10));
      a1 = p(7)+p(10)+p(4)+p(13);
      b1 = a1^2-(p(4)-p(13))^2;
      c1 = (p(7)-p(10))*a1;
      d1 = (p(4)+p(13))^2*(p(7)-p(10))^2+8*p(4)*p(13)*(2*p(4)*p(13)+(p(4)+p(13))*(p(7)-p(10))+2*p(7)*p(10));
      q(7) = (b1+c1+sqrt(d1))/(4*(a1-p(10)));
      q(10) = (b1-c1-sqrt(d1))/(4*(a1-p(7)));
      q(4) = p(4) + (p(7)+p(10)-q(7)-q(10))/2;
      q(13) = p(13) + (p(7)+p(10)-q(7)-q(10))/2;
      % P-SSR computation, general formula (pure P-SSR)
      p(6)  = max(RDM(6,6),RDM(11,11));
      p(11) = min(RDM(6,6),RDM(11,11));
      a2 = p(6)+p(11)+p(1)+p(16);
      b2 = a2^2-(p(1)-p(16))^2;
      c2 = (p(6)-p(11))*a2;
      d2 = (p(1)+p(16))^2*(p(6)-p(11))^2+8*p(1)*p(16)*(2*p(1)*p(16)+(p(1)+p(16))*(p(6)-p(11))+2*p(6)*p(11));
      q(6) = (b2+c2+sqrt(d2))/(4*(a2-p(11)));
      q(11) = (b2-c2-sqrt(d2))/(4*(a2-p(6)));
      q(1) = p(1) + (p(6)+p(11)-q(6)-q(11))/2;
      q(16) = p(16) + (p(6)+p(11)-q(6)-q(11))/2;
      % Final relative entropy
      E = 0; %indices = [4,7,10,13,1,6,11,16];
      for i=1:length(p)%indices
         if(p==0)
            continue
         end
         E = E + p(i)*log2(p(i)/q(i));
      end
   end
end

function [pRDM,nRDM] = filter_RDM(RDM)
   % Global singlet assertion
   assert(all(abs(RDM(4,4)-RDM(13,13))<1e-12))
   % Particle-hole symmetry assertion
   %assert(all(abs(RDM(1,1)-RDM(16,16))<1e-4))
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
 