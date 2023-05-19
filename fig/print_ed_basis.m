function print_ed_basis(Nlat,fancy)
    if nargin < 2
       fancy = true;
    end
    Nspin=2; Norb=1; % These are enforced!
    for state = 0:1:(2^(Nlat*Norb*Nspin)-1)
       label = build_ket(state);
       fprintf('%d\t',state+1)
       disp(label)
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
      for ilat = 1:Nlat
          for ispin = 1:2
              shift = (ilat-1)*Norb + (ispin-1)*Norb*Nlat;
              index = shift+(1:Norb);
              vec(index) = bitget(state,index);
          end
      end
      kup = num2str(vec(1:Norb*Nlat));
      kdw = num2str(vec(Norb*Nlat+1:end));
      if fancy
        ket = ['| ',strrep(kup,'1','↑'),' 〉⊗ ',...
               '| ',strrep(kdw,'1','↓'), ' 〉'];
        ket = strrep(ket,'0','•');
      else
        ket = [kup,'  ',kdw];
      end
    end
end