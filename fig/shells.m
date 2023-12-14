set(0,'defaulttextinterpreter','latex')
%magicLaTeX;

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/Carlos/LATEST/')

tabs = readtable('4x2_dimer_MIs.dat');
uloc = tabs.Var1; dist = [1,sqrt(2),2,sqrt(5),3,sqrt(10)];
data = table2array(tabs); data = data(:,2:end);

cd(HERE)

cd('../../Data/CDMFT/4sites2replicas/')

Uloc = load('U_list.txt');

s1 = QcmP.post.eentropy_line('U','1sites');
s2 = QcmP.post.eentropy_line('U','2sites');

MI =  2.*s1 - s2;

cd(HERE)

%% Actual graphics

QcmP.plot.import_colorlab

plot3(ones(size(Uloc)),Uloc,MI,...
    '.:','LineWidth',2,'MarkerSize',21,...
    'Color',str2rgb('matlab4'))

hold on

for i = 1:length(dist)
    plot3(ones(size(uloc))*dist(i),uloc,data(:,i),...
        '*','MarkerSize',13,'LineWidth',2,...
        'Color',str2rgb('pyplot3'));
    plot3(ones(size(uloc))*dist(i),uloc,data(:,i),...
        '-','LineWidth',2,...
        'Color',str2rgb('pyplot3'));
end

[X,Y] = meshgrid(uloc,dist);
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
matlab2tikz('shells.tex','strict',true,'noSize',true)
rmpath([HERE,'/../lib/m2tex/src']);

%% Clear Secli's black magic
%clear magicLaTeX; 
