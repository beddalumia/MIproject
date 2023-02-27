U = load('U_list.txt');

full = postDMFT.eentropy_line('1sites');

[pSSR,nSSR] = build_SSRs(full);

plot(U,full)
hold on
plot(U,pSSR)
plot(U,nSSR)

function [pEE,nEE] = build_SSRs(fullEE)

   [pmold,UDIR] = postDMFT.get_list('U')

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

   pEE = (p1+p4).*log(p1+p4) + (p2+p3).*log(p2+p3)
   pEE = pEE/log(2) + fullEE

   nEE = (p2+p3).*log(p2+p3) - p2.*log(p2) - p3.*log(p3)
   nEE = nEE/log(2)

end

