%set(0,'defaulttextinterpreter','latex')

% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../../Data/CDMFT/4sites2replicas')

Ulist = QcmP.post.get_list('U');

Uloc=input('Which value of U? ','s');

while not(any(str2double(Uloc)==Ulist))
     disp("U must be in the following set:")
     disp(Ulist);
     Uloc=input('Which value of U? ','s');
end

cd(sprintf('U=%f',str2double(Uloc)));

[bras,kets] = basis_labels;

%% CDMFT

RDM = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');

[pRDM,nRDM] = filter_RDM(RDM);

[F,idx] = DGAbasis(RDM);
[P,idx] = DGAbasis(pRDM);
[N,idx] = DGAbasis(nRDM);

figure("name","CDMFT - Symmetries",'Renderer', 'painters', 'Position', [10 10 900 600])
fancy_spy(F,P,N,kets(idx),bras(idx),idx,[2,6,12,16]);

cd(HERE)

%% Actual Plotting

function fancy_spy(RDM,pRDM,nRDM,kets,bras,indices,blocks)

    RDM = sign(abs(RDM));
    DIAG = sign(abs(diag(diag(RDM))));
    pRDM = sign(abs(pRDM));
    nRDM = sign(abs(nRDM));

    QcmP.plot.import_colorlab
    % Trick to enable clever multi-axes
    t = tiledlayout(1,1);
    % Full RDM
    f_ax = axes(t);
    full = imagesc(RDM+pRDM+nRDM+DIAG);  %full.AlphaDataMapping = 'scaled';
    xlim([0.5,16.5]); xticks([1:16]); xticklabels([]);%xticklabels(string(indices));
    ylim([0.5,16.5]); yticks([1:16]); yticklabels([]);%yticklabels(string(indices));
    axis square; f_ax.Color = "None";
    % Tick magic :)
    set(f_ax,'yaxisLocation','right');
    set(f_ax,'TickLength',[0 0])
    set(f_ax,'TickLabelInterpreter','latex')
    % Symmetry Grid (if any)
    for i = 1:16
        if any(i==blocks)
            xline(i-0.5,'-','Color',str2rgb('grey'),"Linewidth",1)
            yline(i-0.5,'-','Color',str2rgb('grey'),"Linewidth",1)
        else
            xline(i-0.5,':')
            yline(i-0.5,':')
        end
    end
    % Additional tick labels (pretty print of states)
    d_ax = axes(t);
    xlim([0.5,16.5]); xticks([1:16]); xticklabels(bras);
    ylim([0.5,16.5]); yticks([1:16]); yticklabels(fliplr(kets));
    axis square; d_ax.Color = "None";
    % Tick magic :)
    set(d_ax,'xaxisLocation','top')
    set(d_ax,'TickLength',[0 0])
    set(d_ax,'TickLabelInterpreter','latex')

    % Global Colormap
    %colormap([str2rgb("w");str2rgb("wheat");str2rgb("strawberry");str2rgb("sky");str2rgb("k")])
    %colormap([str2rgb("w");str2rgb("y");str2rgb("m");str2rgb("c");str2rgb("k")])
    colormap([str2rgb("w");str2rgb("light grey");str2rgb("light grey");str2rgb("DE4968");str2rgb("k")])

end

%% Utilities

function [rho_dga,idx] = DGAbasis(rho_edi)

    rho_dga = zeros(size(rho_edi));

    % Diagonal terms
    rho_dga(1,1)   = rho_edi(1,1);
    rho_dga(2,2)   = rho_edi(9,9);
    rho_dga(3,3)   = rho_edi(3,3);
    rho_dga(4,4)   = rho_edi(5,5);
    rho_dga(5,5)   = rho_edi(2,2);
    rho_dga(6,6)   = rho_edi(11,11);
    rho_dga(7,7)   = rho_edi(13,13);
    rho_dga(8,8)   = rho_edi(10,10);
    rho_dga(9,9)   = rho_edi(7,7);
    rho_dga(10,10) = rho_edi(4,4);
    rho_dga(11,11) = rho_edi(6,6);
    rho_dga(12,12) = rho_edi(15,15);
    rho_dga(13,13) = rho_edi(12,12);
    rho_dga(14,14) = rho_edi(14,14);
    rho_dga(15,15) = rho_edi(8,8);
    rho_dga(16,16) = rho_edi(16,16);


    % Off diagonal
    rho_dga(5,3)   = rho_edi(2,3);
    rho_dga(3,5)   = rho_edi(3,2);

    rho_dga(4,2)   = rho_edi(5,9);
    rho_dga(2,4)   = rho_edi(9,5);

    rho_dga(11,9)  = -rho_edi(6,7);
    rho_dga(9,11)  = -rho_edi(7,6);

    rho_dga(11,8)  = rho_edi(6,10);
    rho_dga(8,11)  = rho_edi(10,6);

    rho_dga(11,6)  = rho_edi(6,11);
    rho_dga(6,11)  = rho_edi(11,6);

    rho_dga(9,8)   = -rho_edi(7,10);
    rho_dga(8,9)   = -rho_edi(10,7);

    rho_dga(9,6)   = -rho_edi(7,11);
    rho_dga(6,9)   = -rho_edi(11,7);

    rho_dga(15,13) = -rho_edi(8,12);
    rho_dga(13,15) = -rho_edi(12,8);

    rho_dga(8,6)   = rho_edi(10,11);
    rho_dga(6,8)   = rho_edi(11,10);

    rho_dga(14,12) = -rho_edi(14,15);
    rho_dga(12,14) = -rho_edi(15,14);


    assert(ishermitian(rho_dga))

    idx = [
     1   % 1 → 1
     9   % 2 → 9
     3   % 3 → 3
     5   % 4 → 5
     2   % 5 → 2
    11   % 6 → 11
    13   % 7 → 13
    10   % 8 → 10
     7   % 9 → 7
     4   % 10 → 4
     6   % 11 → 6
    15   % 12 → 15
    12   % 13 → 12
    14   % 14 → 14
     8   % 15 → 8
    16   % 16 → 16
];


end




function [pRDM,nRDM] = filter_RDM(RDM)
    % Trim all off-diagonals
    nRDM = diag(diag(RDM));
    % Restore the spin-flip terms
    nRDM(7,10) = RDM(7,10);
    nRDM(10,7) = RDM(10,7);
    %
    pRDM = nRDM;
    % Restore the pair-hopping terms
    pRDM(6,11) = RDM(6,11);
    pRDM(11,6) = RDM(11,6);
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
    k_i = num2str(bin2dec(num2str([vec(1),vec(3)])));
    k_j = num2str(bin2dec(num2str([vec(2),vec(4)])));
    ket = ['$\mid ',k_i, '\rangle\,\otimes \mid ' , k_j, ' \rangle$'];
    ket = strrep(ket,'0','\!\bullet\,');
    ket = strrep(ket,'1','\downarrow');
    ket = strrep(ket,'2','\uparrow');
    ket = strrep(ket,'3','\,\uparrow\downarrow\,');
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
    b_i = num2str(bin2dec(num2str([vec(1),vec(3)])));
    b_j = num2str(bin2dec(num2str([vec(2),vec(4)])));
    bra = ['$\langle ',b_i, '\mid \otimes\,\langle ', b_j, ' \mid$'];
    bra = strrep(bra,'0','\,\bullet\!');
    bra = strrep(bra,'1','\downarrow');
    bra = strrep(bra,'2','\uparrow');
    bra = strrep(bra,'3','\,\uparrow\downarrow\,');
end
%
function [bras,kets] = basis_labels()
    bras = strings(1,16); kets = bras;
    for state = 0:1:15
        bras(state+1) = build_bra(state);
        kets(state+1) = build_ket(state);
    end
end
function [ROT,new_indices] = rotate(RDM_ij)
%
    % SSR-friendly sectors
    i0 = [4,7,10,13];
    i1 = [2,3,5,8,9,12,14,15];
    i2 = [1,6,11,16];
%
    % Compose and rotate
    new_indices = [i0,i1,i2];
    ROT = RDM_ij(new_indices,new_indices);
%
    % Fermionic signs
    % > unimplemented
%
end
function [ROT,new_indices] = imbalance(RDM_ij)
%
    % Fixed-charge imbalance states
    i0 = [1,4,7,10,13,16];
    i1 = [2,3,5,8,9,12,14,15];
    i2 = [6,11];
%
    % Compose and rotate
    new_indices = [i0,i1,i2];
    ROT = RDM_ij(new_indices,new_indices);
%
    % Fermionic signs
    % > unimplemented
%
end
function [ROT,new_indices] = symmetrize(RDM_ij)
% Put the RDM_ij in block-diagonal form, according
% to total charge and magnetization sectors
%
% 1       | •  • 〉⊗ | •  • 〉
% 2       | ↑  • 〉⊗ | •  • 〉
% 3       | •  ↑ 〉⊗ | •  • 〉
% 4       | ↑  ↑ 〉⊗ | •  • 〉
% 5       | •  • 〉⊗ | ↓  • 〉
% 6       | ↑  • 〉⊗ | ↓  • 〉
% 7       | •  ↑ 〉⊗ | ↓  • 〉
% 8       | ↑  ↑ 〉⊗ | ↓  • 〉
% 9       | •  • 〉⊗ | •  ↓ 〉
% 10      | ↑  • 〉⊗ | •  ↓ 〉
% 11      | •  ↑ 〉⊗ | •  ↓ 〉
% 12      | ↑  ↑ 〉⊗ | •  ↓ 〉
% 13      | •  • 〉⊗ | ↓  ↓ 〉
% 14      | ↑  • 〉⊗ | ↓  ↓ 〉
% 15      | •  ↑ 〉⊗ | ↓  ↓ 〉
% 16      | ↑  ↑ 〉⊗ | ↓  ↓ 〉
%
%
    % Fixed-charge sectors: [nup,ndw]
    i00 = [1];
    i10 = [2,3];
    i01 = [5,9];
    i20 = [4];
    i11 = [6,7,10,11];
    i02 = [13];
    i21 = [8,12];
    i12 = [14,15];
    i22 = [16];
%
    % Compose and rotate
    new_indices = [i00,i10,i01,i20,i11,i02,i21,i12,i22];
%
    ROT = RDM_ij(new_indices,new_indices);
%
end
function [PT,new_indices] = partial_transpose(RDM_ij)
% No Fermi signs, no phases! Just move around elements
    [Nrdm,Mrdm] = size(RDM_ij); assert(Nrdm==Mrdm);
    ROT = zeros(Nrdm,Mrdm); PT = ROT;
    % Rotate the RDM_ij to the proper block form
    [ROT,new_indices] = single_site_sectors(RDM_ij);
    % Finally, the partial transpose on each "left site" block
    N_1site = 4;
    for i = 1:Nrdm/N_1site
       for j = 1:Nrdm/N_1site
          block = ROT(1+(i-1)*N_1site:i*N_1site,1+(j-1)*N_1site:j*N_1site);
          PT(1+(i-1)*N_1site:i*N_1site,1+(j-1)*N_1site:j*N_1site) = block';
       end
    end
end
function [ROT,new_indices] = single_site_sectors(RDM_ij)
    [Nrdm,Mrdm] = size(RDM_ij); assert(Nrdm==Mrdm);
    Nlat = 2;
    Norb = 1;
    Nlso = 4; % 2 sites x 2 spin x 1 orbital
    ROT = zeros(Nrdm,Mrdm); PT = ROT;
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
       new_indices(i) = newI;
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
             ROT(newI,newJ) = RDM_ij(i,j);
       end
    end
end

