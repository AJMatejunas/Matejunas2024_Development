%Noise Plot

noiseCopy=1:30;
NoiseEnds=[1,30];

figure
    plot(noiseCopy,ConstitutiveParam(:,1),'b');
    hold on
    plot(NoiseEnds,RefPar.MatProps.Ki*[1,1],'r')
    plot(NoiseEnds,RefPar.Kp10*[1,1],'--r')
    plot(NoiseEnds,RefPar.Km10*[1,1],'--r','HandleVisibility','off')
    hold off

    title('K_1 identification')
    xlabel('Noise Copy ')
    ylabel('K_1 (Pa)')
    legend('Identified','reference','reference \pm 10\%')
    
    figure
    plot(noiseCopy,Errors.K,'b');
    hold on
    plot(NoiseEnds,[0,0],'r')
    plot(NoiseEnds,[10,10],'--r')
    plot(NoiseEnds,[-10,-10],'--r','HandleVisibility','off')
    hold off

    title('K_1 identification Errors')
    xlabel('Noise Copy ')
    ylabel('Errors (\%)')
    legend('Identification','reference','\pm 10\%')