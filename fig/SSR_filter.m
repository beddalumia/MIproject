function [pRDM,nRDM] = SSR_filter(RDM)
    % Trim all off-diagonals
    nRDM = diag(diag(RDM));
    % Restore the spin-flip terms
    nRDM(7,10) = RDM(7,10);
    nRDM(10,7) = RDM(10,7);
    pRDM = nRDM;
    % Restore the pair-hopping terms
    pRDM(6,11) = RDM(6,11);
    pRDM(11,6) = RDM(11,6);
end