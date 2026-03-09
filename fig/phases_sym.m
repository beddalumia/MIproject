set(0,'defaulttextinterpreter','latex')

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

U = [1.7,2.3,3,4,5,6,7]; Nu = length(U);
%p = palette.crameri('berlin'); Np = length(p);

figure

for i = 1:Nu

    cd(sprintf('../../Data/CDMFT/DOPED_MOTT/Uloc%g',U(i)))
    pwd

    load('Emix.txt')
    load('E1.txt')
    load('mu.txt')
    %QcmP.post.observables_line('xmu','site001');
    dens = load('dens__site001.txt')
    mott = mu(dens>0.9995);
    liqu = mu(E1<1e-2);



    scatter(mu,U(i)*ones(size(dens)),[],Emix,'filled'); hold on
    set_palette('PuRd')
    %colormap(p(Np/2+1:end,:))
    colorbar('northoutside')
    caxis([0,0.03])
    p = colormap;
    %set(gca,'Color',p(1,:))
    xlim([-5,0])
    ylim([-0.5,7.5])
    scatter(0,1.55,100,'*','MarkerEdgeColor',str2rgb('goldenrod'),'Linewidth',2)
    scatter(0,1.55,70,'x','MarkerEdgeColor',str2rgb('green'),'Linewidth',2)
    try
    scatter(mott(1),U(i),100,'*','MarkerEdgeColor',str2rgb('goldenrod'),'Linewidth',2)
    end
    scatter(liqu(end),U(i),70,'x','MarkerEdgeColor',str2rgb('green'),'Linewidth',2)
    box on
    xlabel('$\mu/D$')
    ylabel('$U/D$')
    axis square

    cd(HERE)

end

figure

for i = 1:Nu

    cd(sprintf('../../Data/CDMFT/DOPED_MOTT/Uloc%g',U(i)))
    pwd

    load('Emix.txt')
    load('E1.txt')
    load('mu.txt')
    %QcmP.post.observables_line('xmu','site001');
    dens = load('dens__site001.txt')
    mott = mu(dens>0.9995);
    liqu = mu(E1<1e-2);

    scatter(mu,U(i)*ones(size(dens)),[],E1,'filled'); hold on
    set_palette('Blues')
    %colormap(flipud(p(1:Np/2,:)))
    colorbar('northoutside')
    caxis([0,0.3])
    p = colormap;
    %set(gca,'Color',p(1,:))
    xlim([-5,0])
    ylim([-0.5,7.5])
    scatter(0,1.55,100,'*','MarkerEdgeColor',str2rgb('goldenrod'),'Linewidth',2)
    scatter(0,1.55,70,'x','MarkerEdgeColor',str2rgb('green'),'Linewidth',2)
    try
    scatter(mott(1),U(i),100,'*','MarkerEdgeColor',str2rgb('goldenrod'),'Linewidth',2)
    end
    scatter(liqu(end),U(i),70,'x','MarkerEdgeColor',str2rgb('green'),'Linewidth',2)
    box on
    xlabel('$\mu/D$')
    ylabel('$U/D$')
    axis square



    cd(HERE)

end



