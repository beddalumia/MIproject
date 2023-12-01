set(0,'defaulttextinterpreter','latex')

% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/4sites2replicas')

Ulist = postDMFT.get_list('U');

Uloc=input('Which value of U? ','s');

while not(any(str2double(Uloc)==Ulist))
     disp("U must be in the following set:")
     disp(Ulist);
     Uloc=input('Which value of U? ','s');
end

cd(sprintf('U=%f',str2double(Uloc)));

[bras,kets] = basis_labels;

RDM = postDMFT.get_Hloc('reduced_density_matrix_2sites.dat');

[pRDM,nRDM] = filter_RDM(RDM);

cd(HERE)

%% Actual Plotting
QcmP.plot.import_colorlab
% Trick to enable clever multi-axes
t = tiledlayout(1,1);
% Full RDM
f_ax = axes(t); 
full = imagesc(sign(abs(RDM)));  full.AlphaDataMapping = 'scaled';
xlim([0.5,16.5]); xticks([1:16]); xticklabels(string([1:16]));
ylim([0.5,16.5]); yticks([1:16]); yticklabels(string([1:16]));
axis square; f_ax.Color = "None"; 
% Tick magic :)
set(f_ax,'yaxisLocation','right');
set(f_ax,'TickLength',[0 0])
set(f_ax,'TickLabelInterpreter','latex')
% P-SSR RDM
p_ax = axes(t); 
pssr = imagesc(sign(abs(pRDM))); pssr.AlphaDataMapping = 'scaled';
xlim([0.5,16.5]); xticks([]); % Disallow ticks to
ylim([0.5,16.5]); yticks([]); % avoid cluttering…
axis square; p_ax.Color = "None";
% N-SSR RDM
n_ax = axes(t); 
nssr = imagesc(sign(abs(nRDM))); nssr.AlphaDataMapping = 'scaled';
xlim([0.5,16.5]); xticks([1:16]); xticklabels(bras);
ylim([0.5,16.5]); yticks([1:16]); yticklabels(kets);
axis square; n_ax.Color = "None";
% Tick magic :)
set(n_ax,'xaxisLocation','top')
set(n_ax,'TickLength',[0 0])

% Global Colormap
set_palette("Greys")

% Global Grid
for i = 1:16
    xline(i-0.5)
    yline(i-0.5)
end

%% Utilities
function [pRDM,nRDM] = filter_RDM(RDM)
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
if(any(AB~=new))
    disp(AB)
    disp(new)
end
nRDM(7,7) = new(1,1);
nRDM(10,10) = new(2,2);
%
pRDM = nRDM;
% Restore the pair-hopping terms
pRDM(6,11) = RDM(6,11);
pRDM(11,6) = RDM(11,6);
% Global singlet assertion
assert(all(abs(RDM(4,4)-RDM(13,13))<1e-12))
% Particle-hole symmetry assertion
assert(all(abs(RDM(1,1)-RDM(16,16))<1e-4))
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
ket = ['\mid ',strrep(kup,'1','\uparrow'),' \rangle \otimes ',...
            '\mid ',strrep(kdw,'1','\downarrow'), ' \rangle'];
ket = strrep(ket,'0','\bullet');
end
function bra = build_bra(state)
%% BUILD_BRA : Puts together a pretty label for a pure state component
%
%  >> bra = build_bra(state)
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
bup = num2str(vec(1:Norb*Nlat));
bdw = num2str(vec(Norb*Nlat+1:end));
bra = ['\langle ',strrep(bup,'1','\uparrow'),' \mid  \otimes ',...
            '\langle ',strrep(bdw,'1','\downarrow'), ' \mid'];
bra = strrep(bra,'0','\bullet');
end
%
function [bras,kets] = basis_labels()
bras = strings(1,16); kets = bras;
for state = 0:1:15
    bras(state+1) = build_bra(state);
    kets(state+1) = build_ket(state);
end
end
