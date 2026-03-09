HERE = pwd;

cd ../../Data/CDMFT/doped/Uloc2.3/
mu = QcmP.post.get_list('xmu');

SDATA = cell(length(mu),1);
for i = 1:length(mu)
   cd(sprintf("xmu=%f",mu(i)))
   % S11 = QcmP.plot.spectral_load('impSigma_Isite0001_Jsite0001_l11_s1_realw.ed');
   % S12 = QcmP.plot.spectral_load('impSigma_Isite0001_Jsite0002_l11_s1_realw.ed');
   % S13 = QcmP.plot.spectral_load('impSigma_Isite0001_Jsite0003_l11_s1_realw.ed');
   % S14 = QcmP.plot.spectral_load('impSigma_Isite0001_Jsite0004_l11_s1_realw.ed');
   % w = S11.zeta;
   % d = abs(w); %abs(w+mu(i)); 
   % [m,index] = min(d);
   % RS11(i) = S11.real(index);
   % RS12(i) = S12.real(index);
   % RS13(i) = S13.real(index);  
   % RS14(i) = S14.real(index);
   % IS11(i) = S11.imag(index);
   % IS12(i) = S12.imag(index);
   % IS13(i) = S13.imag(index);  
   % IS14(i) = S14.imag(index);
   % O11 = S11.imag(index);
   % O12 = S12.imag(index);
   % O13 = S13.imag(index);
   % O14 = S14.imag(index);
   % Ozz(i) = O11+O12+O13+O14;
   % Ozp(i) = O11+O12-O13-O14;
   % Opz(i) = O11-O12+O13-O14;
   % Opp(i) = O11-O12-O13+O14;
   Nlat = 4;
   S = cell(Nlat,Nlat);
   for ilat = 1:Nlat
      for jlat = 1:Nlat
         S{ilat,jlat} = QcmP.plot.spectral_load(...
            sprintf('impSigma_Isite000%d_Jsite000%d_l11_s1_realw.ed',ilat,jlat));
      end
   end
   SDATA{i} = S;
   cd ..
end
cd(HERE)

imu = find(mu==-0.6);

   % plot(mu,Ozz); hold on
   % plot(mu,Opz);
   % plot(mu,Ozp,"--");
   % plot(mu,Opp);
   % legend(["(0,0)","(\pi,0)","(0,\pi)","(\pi,\pi)"]);

   %kx = -pi:0.01:pi;
   %ky = -pi:0.01:pi;

   Nk = 1024;
   MxBZ = fix(sqrt(Nk)); MyBZ = MxBZ;

   k_BZ = zeros(Nk,2);

   for ik = 1:Nk
      ikx = mod((ik-1), MxBZ);
      iky = floor((ik-1) / MxBZ);
      k_BZ(ik,1) = pi*(-1 + 2*(ikx+0.5)/MxBZ);
      k_BZ(ik,2) = pi*(-1 + 2*(iky+0.5)/MyBZ);
   end
   kx = unique(k_BZ(:,1));
   ky = unique(k_BZ(:,2));

   for i=1:length(kx)
      for j=1:length(ky)
         Sk = periodized_sigma(kx(i),ky(j),SDATA{imu});
         RSk_map(i,j) = real(Sk);
         ISk_map(i,j) = imag(Sk);
      end
   end

   % Compute noninteracting dispersion (folded!)
   t = 0.25;
   NRBZ = Nk/4;
   k_RBZ = zeros(NRBZ,2);
   for ik = 1:NRBZ
      ix = mod(ik-1, MxBZ);
      iy = floor((ik-1) / MxBZ);
      k_RBZ(ik,1) = 0.5*pi*(-1 + 2*ix/MxBZ);
      k_RBZ(ik,2) = 0.5*pi*(-1 + 2*iy/MyBZ);
   end
   kx_RBZ = k_RBZ(:,1);
   ky_RBZ = k_RBZ(:,2);
   eps_k = -2*t*(cos(kx_RBZ) + cos(ky_RBZ));

   % Compute interacting Green's function at omega=0


kx = unique(kx);
ky = unique(ky);
figure;
subplot(1,2,1);
imagesc(kx,ky,RSk_map(:,:));
set(gca,'YDir','normal');
colorbar; set_palette('inferno');
title('Re \Sigma(k,\omega=0)');
xlabel('k_x');
ylabel('k_y');
axis equal;
xlim([-pi pi]);
ylim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
subplot(1,2,2);
imagesc(kx,ky,ISk_map(:,:));
set(gca,'YDir','normal');
colorbar; set_palette('inferno');
title('Im \Sigma(k,\omega=0)');
xlabel('k_x');
ylabel('k_y');
axis equal;
xlim([-pi pi]);
ylim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
yticks([-pi -pi/2 0 pi/2 pi]);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'}); 

function Sk = periodized_sigma(kx,ky,S)

   % RSk = 0.25*(RS11 + RS12*exp(-1i*kx) + RS13*exp(-1i*ky) + RS14*exp(-1i*(kx+ky)));
   % ISk = 0.25*(IS11 + IS12*exp(-1i*kx) + IS13*exp(-1i*ky) + IS14*exp(-1i*(kx+ky)));
   % Sk = RSk + 1i*ISk;
   w = S{1,1}.zeta;
   d = abs(w); %abs(w+mu(i)); 
   [m,index] = min(d);
   Sk = 0;
   for ilat = 1:4
      for jlat = 1:4
         delta_ij = [mod(ilat-1,2)-mod(jlat-1,2), floor((ilat-1)/2)-floor((jlat-1)/2)];
         Sk = Sk  + S{ilat,jlat}.real(index) + ...
                 1j*S{ilat,jlat}.imag(index) * ...
         exp(-1i*(kx*delta_ij(1) + ky*delta_ij(2)));
      end
   end

end