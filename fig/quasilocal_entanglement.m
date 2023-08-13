set(0,'defaulttextinterpreter','latex')

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/Carlos/LATEST/')

%temp = readtable('SingleSite.dat');
temp = readtable('DimerVert.dat');
Uloc = temp.Var1;
S1x2 = temp.Var2;
S2x2 = temp.Var3;
S3x2 = temp.Var4;
S4x2 = temp.Var5;
%temp = readtable('Plaquette.dat');

cd(HERE)

cd('../../Data/CDMFT/2sites5replicas/')

U2 = load('U_list.txt');
S2 = postDMFT.eentropy_line('2sites');

cd(HERE)

cd('../../Data/CDMFT/4sites2replicas/')

U4 = load('U_list.txt');
S4 = postDMFT.eentropy_line('2sites');

cd(HERE)

%% Actual graphics

plotDMFT.import_colorlab

figure("Name",'Local Entropy Scaling')
%plot(U2,S2,'-','LineWidth',1.5,'Color',str2rgb('fuchsia'))
plot(U4,S4,'-*','LineWidth',1.5,'Color',str2rgb('matlab4'))
hold on
%plot(Uloc,S1x2,'o','LineWidth',1.5,'MarkerSize',10,'Color',str2rgb('pyplot2'))
plot(Uloc,S2x2,'x','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('00A4EF'))
plot(Uloc,S3x2,'s','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('FFB900'))
plot(Uloc,S4x2,'^','LineWidth',1.5,'MarkerSize',10,'Color',hex2rgb('7FBA00'))
xlabel("$U/D$")
ylabel("Units of $\log(2)$")
set(gca,'FontSize',15)
legend([...%"ED $1\times2$, $N_\mathrm{bath}=10$",...
        "ED $2\times2$, $N_\mathrm{bath}=8$",...
        ...%"ASCI $1\times2$, $N_\mathrm{bath}=12$",...
        "ASCI $2\times2$, $N_\mathrm{bath}=24$",...
        "ASCI $3\times2$, $N_\mathrm{bath}=36$",...
        "ASCI $4\times2$, $N_\mathrm{bath}=48$"],...
        'Interpreter','latex'); legend('boxoff')

%% Export to TikZ
matlab2tikz('quasilocal_scaling.tex','strict',true,'noSize',true)
