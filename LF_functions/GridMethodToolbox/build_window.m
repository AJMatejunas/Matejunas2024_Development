function g = build_window( profile, sigma )
% build window analysis
% 1) input 
% profile : 0=Gaussian, 1=bi-triangular, 2=triangular-rectangular, 3=bi-rectangular
% sigma: 'half-width' of the analysis window
% 2) output
% g: analysis window 

switch profile
    
    case 0
        t_noy=ceil(4*sigma);  
        [X,Y] = meshgrid(-t_noy:t_noy, -t_noy:t_noy);  
        g=1/(2*pi*sigma^2).*exp((-X.^2-Y.^2)/(2*sigma.^2));
        g=g/sum(g(:));

    case 1
        t_noy=floor(sigma);  
        g=triang(2*t_noy-1)*triang(2*t_noy-1)';
        g=g/sum(g(:));
  
    case 2
        t_noy=floor(sigma);
        rect=ones(2*t_noy,1);
        g=rect*triang(2*t_noy-1)';
        g=g/sum(g(:));
      
    case 3
        t_noy=floor(sigma);
        rect=ones(2*t_noy,1);
        g=rect*rect';
        g=g/sum(g(:));
end


end

