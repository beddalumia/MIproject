set(0,'defaulttextinterpreter','latex')

%% LIBRARIES
QcmP.plot.import_colorlab
addpath ../lib/m2tex/src

%% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/CDMFT/doped/Uloc2.3';
cd(DATA)

[mu,mudir] = QcmP.post.get_list('xmu')

nup = load('nup_1_site001.txt');
ndw = load('ndw_1_site001.txt');
dens = nup + ndw;

assert(all(abs(dens-load('dens__site001.txt'))<1e-8));

for i = 1:length(mu)
   cd(mudir(i))
   RDMi = load('reduced_density_matrix_1sites.dat');
   holi(i) = RDMi(1,1);
   docc(i) = RDMi(4,4);
   cd('..')
end

assert(all(abs(docc'==load('docc_1_site001.txt'))<1e-8));


figure("Name",'Hole-driven Mott transition')
% tiledlayout(1,2)
% nexttile
plot(mu,dens,':s','MarkerSize',5,'MarkerEdgeColor',str2rgb('dark sky blue'),...
   'MarkerFaceColor',str2rgb('pale sky blue'),'Linewidth',0.5);
xlim([-1,0]); ylim([0.6,1.05]); box on;
xlabel('$\mu/D$','Interpreter','latex');
ylabel('$\langle n_{i} \rangle$','Interpreter','latex');
% nexttile
% plot(mu,docc,':d','MarkerSize',3,'MarkerEdgeColor',str2rgb('grass'),...
%    'MarkerFaceColor',str2rgb('yellowish green'),'Linewidth',0.5);
% xlim([-1,0]); ylim([0.015,0.04]); box on;
% xlabel('$\mu/D$','Interpreter','latex');
% ylabel('$\langle n_{i\uparrow}n_{i\downarrow}\rangle$','Interpreter','latex');

% Export to TikZ
matlab2tikz('filename',[CODE,'/doping.tex'],'width','.8\textwidth','height','.6\textwidth');


%close all
cd(CODE);

%% Reset path
rmpath ../lib/m2tex/src

%% Utilities

function E = vonNeumann(RDM)
   p = eig(RDM);
   E = -sum(p.*log2(p))
end