%% Activate Secli's black magic
magicLaTeX;

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/4sites2replicas/')

Ulist = postDMFT.get_list('U');

Uloc=input('Which value of U? ','s');

while not(any(str2double(Uloc)==Ulist))
     disp("U must be in the following set:")
     disp(Ulist);
     Uloc=input('Which value of U? ','s');
end

cd(sprintf('U=%f',str2double(Uloc)));

%% Actual plotting

% Invoke DMFT-LAB utilities
rdm1 = figure("Name","Single Site RDM");
plotDMFT.pure_states('1sites')
rdm2 = figure("Name","Nearest-Neighbor RDM");
plotDMFT.pure_states('2sites')
% Reset to base directory
cd(HERE)
% Export to PDF a la Secli
set(rdm1,'Renderer','painters');
saveFigureAsPDF(rdm1,'rdm1.pdf');
set(rdm2,'Renderer','painters');
saveFigureAsPDF(rdm2,'rdm2.pdf');

%% Reset Secli's black magic
clear magicLaTeX; 
