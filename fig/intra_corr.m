set(0,'defaulttextinterpreter','latex')

%% LIBRARIES
QcmP.plot.import_colorlab
addpath ../lib/m2tex/src

%% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/CDMFT/doped/Uloc2.3';
cd(DATA)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);
Umat = repmat(0:0.1:15,Nso,1);
Imat = zeros(size(Umat));

[mu,mudir] = QcmP.post.get_list('xmu')

si = QcmP.post.eentropy_line('xmu','1sites');
su = zeros(size(si));
sd = zeros(size(si));

nup = load('nup_1_site001.txt');
ndw = load('ndw_1_site001.txt');

for i=1:length(mu)
   RDMu=[(1-nup(i)),0;0,nup(i)]; su(i) = vonNeumann(RDMu);
   RDMd=[(1-ndw(i)),0;0,ndw(i)]; sd(i) = vonNeumann(RDMd);
end

intra = su+sd-si;


figure("Name",'Doped Intra-orbital correlation')
hold on
plot(mu,si,':','Color',str2rgb('mauve'),'Linewidth',1.2);
plot(mu,si,'s','MarkerSize',5,'MarkerEdgeColor',str2rgb('mauve'),...
   'MarkerFaceColor',str2rgb('light rose'),'Linewidth',0.5);
annotation('textbox',[0.18,0.88,0.1,0.05],'string','$s_i$',...
   'Interpreter','latex','EdgeColor','none')
plot(mu,su,'-','Color',str2rgb('watermelon'),'Linewidth',2);
plot(mu,sd,'--','Color',str2rgb('azure'),'Linewidth',2);
annotation('textbox',[0.18,0.57,0.1,0.05],'string','$s_{i,\uparrow}=s_{i,\downarrow}$',...
   'Interpreter','latex','EdgeColor','none')
fill([mu;flipud(mu)],[intra;zeros(size(intra))],str2rgb("light khaki"),...
   'EdgeColor','none')
plot(mu,intra,':','Color',str2rgb('grass'),'Linewidth',1.2);
plot(mu,intra,'d','MarkerSize',3,'MarkerEdgeColor',str2rgb('grass'),...
   'MarkerFaceColor',str2rgb('yellowish green'),'Linewidth',0.5);
annotation('textbox',[0.18,0.19,0.1,0.05],'string','$I(\,\uparrow\,:\,\downarrow\,)$',...
   'Interpreter','latex','EdgeColor','none')
xlim([-1,0]); ylim([0,2]); box on;
xlabel('$\mu/D$','Interpreter','latex');
ylabel('[bit]','Interpreter','latex');

% Export to TikZ
matlab2tikz('filename',[CODE,'/doped_intra_corr.tex'],'width','\textwidth','height','1.2\textwidth');


%close all
cd(CODE);

%% Reset path
rmpath ../lib/m2tex/src

%% Utilities

function E = vonNeumann(RDM)
   p = eig(RDM);
   E = -sum(p.*log2(p))
end