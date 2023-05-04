set(0,'defaulttextinterpreter','latex')

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/Carlos/EEs/')

temp = readtable('SingleSite.dat');
Uloc = temp.Var1;
S1x2 = temp.Var2;
S2x2 = temp.Var3;
S3x2 = temp.Var4;
S4x2 = temp.Var5;
%temp = readtable('Dimer.dat');
%temp = readtable('Plaquette.dat');

cd(HERE)

cd('../../Data/CDMFT/2sites5replicas/')

U2 = load('U_list.txt');
S2 = postDMFT.eentropy_line('1sites');

cd(HERE)

cd('../../Data/CDMFT/4sites2replicas/')

U4 = load('U_list.txt');
S4 = postDMFT.eentropy_line('1sites');

[pSSR,nSSR] = build_SSRs(S4);

cd(HERE)

%% Actual graphics

plotDMFT.import_colorlab

figure("Name",'Local Entropy Scaling')
plot(U2,S2,'-','LineWidth',1.5,'Color',str2rgb('matlab4'))
hold on
plot(U4,S4,'--','LineWidth',1.5,'Color',str2rgb('fuchsia'))
ylim([1,2])
plot(Uloc,S1x2,'o','LineWidth',1.5,'MarkerSize',10,'Color',str2rgb('pyplot2'))
plot(Uloc,S2x2,'x','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('00A4EF'))
plot(Uloc,S3x2,'s','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('FFB900'))
plot(Uloc,S4x2,'^','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('7FBA00'))
xlabel("$U/D$")
ylabel("Units of $\log(2)$")
set(gca,'FontSize',15)
legend(["ED $1\times2$, $N_\mathrm{bath}=10$",...
        "ED $2\times2$, $N_\mathrm{bath}=8$",...
        "ASCI $1\times2$, $N_\mathrm{bath}=12$",...
        "ASCI $2\times2$, $N_\mathrm{bath}=24$",...
        "ASCI $3\times2$, $N_\mathrm{bath}=36$",...
        "ASCI $4\times2$, $N_\mathrm{bath}=48$"],...
        'Interpreter','latex')
%matlab2tikz('local_scaling.tex','strict',true,'noSize',true)

figure("Name",'Plaquette Reducing Factors')
subplot(2,1,1,'align')
plot(U4,S4)
xlim([0,3])
ylim([1,2])
ylabel("Units of $\log(2)$")
xlabel("")
subplot(2,1,2,'align')
plot(U4,S4./pSSR)
hold on
plot(U4,S4./nSSR)
xlim([0,3])
ylim([0,11])
xlabel("$U/D$")

%% Utilities

function [pEE,nEE] = build_SSRs(fullEE)

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

