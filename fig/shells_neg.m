set(0,'defaulttextinterpreter','latex')
%magicLaTeX;

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

%cd('../../Data/CDMFT/Carlos/LATEST/')
cd('../../Data/CDMFT/Carlos/rho_ij_asci/')

%tabs = readtable('4x2_dimer_MIs.dat');
%uloc = tabs.Var1; 
uloc = [1,8,12,28];
dist = [1,sqrt(2),2,sqrt(5),3,sqrt(10)];
suff = ["_0_4","_0_5","_0_2","_0_6","_0_3","_0_7"];
% Since..
% 4 5 6 7
% 0 1 2 3
%data = table2array(tabs); data = data(:,2:end);

for i = 1:length(uloc)
    cd(sprintf('U%d',uloc(i)))
    for j = 1:length(dist)
        rdm = read_asci_rdm(sprintf('rDM%s.dat',suff(j)));
        U = uloc(i)
        d = dist(j)
        neg = negativity(rdm,true)
        data(i,j) = neg;
    end
    cd('..')
end

cd(HERE)

cd('../../Data/CDMFT/4sites2replicas/')

Uloc = load('U_list.txt');

%s1 = QcmP.post.eentropy_line('U','1sites');
%s2 = QcmP.post.eentropy_line('U','2sites');

%MI =  2.*s1 - s2;

for i = 1:length(Uloc)

    cd(sprintf("U=%f",Uloc(i)))
    rdm = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
    neg = negativity(rdm,false);
    MI(i) = neg;
    cd('..')

end

cd(HERE)

%% Actual graphics

QcmP.plot.import_colorlab

plot3(ones(size(Uloc)),Uloc,MI,...
    '.:','LineWidth',2,'MarkerSize',21,...
    'Color',str2rgb('matlab4'))

hold on

for i = 1:length(dist)
    plot3(ones(size(uloc))*dist(i),uloc/4,data(:,i),...
        '*','MarkerSize',13,'LineWidth',2,...
        'Color',str2rgb('pyplot3'));
    plot3(ones(size(uloc))*dist(i),uloc/4,data(:,i),...
        '-','LineWidth',2,...
        'Color',str2rgb('pyplot3'));
end

[X,Y] = meshgrid(uloc/4,dist);
wplot = waterfall(Y,X,data');

wplot.FaceColor = str2rgb("light khaki");
wplot.FaceAlpha = 0.5;
wplot.EdgeColor = str2rgb('pyplot3');


xlabel("$d/a$")
xlim([1,3.5])
ylabel("$U/D$")
ylim([0.25,7])
zlabel("$I_{ij} / \log(2)$")
zlim([0,1.2])

view(50,36); grid on

legend(["$2\times2$, $N_\mathrm{bath}=08$ (ED)",...
    "$4\times2$, $N_\mathrm{bath}=48$ (ASCI)"],...
    'Location','best',...
    'Interpreter','latex')

%% Save to PDF
%set(gcf,'Renderer','painters');
%saveFigureAsPDF(gcf,'shells.pdf');

%% Export to TikZ
addpath([HERE,'/../lib/m2tex/src']);
matlab2tikz('shells_neg.tex','strict',true,'noSize',true)
rmpath([HERE,'/../lib/m2tex/src']);

%% Clear Secli's black magic
%clear magicLaTeX; 

%% Read Carlos' matrices
function full_matrix = read_asci_rdm(file_rdm)

    asci_matrix = load(file_rdm);
    full_matrix = full(spconvert(asci_matrix(2:end,:)));

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