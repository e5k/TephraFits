function h = plotClassification(type)

h   = figure;
ax  = axes('Parent',h);
hold on
if strcmpi(type, 'Pyle89')
    img = imread('style_Pyle89.jpg');
    img = img(:,:,1);

    xVec    = logspace(-1,2,size(img,2));
    yVec    = linspace(0,4,size(img,1));

    pcolor(xVec,yVec,flipud(img));
    ax.XScale = 'log';
    colormap(gray);
    xlabel('Thickness half-distance b_T (km)');
    ylabel('Half-distanceratio b_C/b_T');
    
elseif strcmpi(type, 'BonadonnaCosta13')
    img = imread('style_BC13.jpg');
    img = img(:,:,1);
    
    xVec    = logspace(0,log10(400),size(img,2));
    yVec    = logspace(-2,2,size(img,1));

    pcolor(xVec,yVec,flipud(img));
    ax.XScale = 'log';
    ax.YScale = 'log';
    colormap(bone);
    xlabel('Log_{10} \lambda_{th}');
    ylabel('Log_{10} \lambda_{MC}/\lambda_{th}');
    
elseif strcmpi(type, 'Mastin09')
    img = imread('style_Mastin09.jpg');
    img = img(:,:,1);
    
    xVec    = linspace(0,60,size(img,2));
    yVec    = logspace(4,10,size(img,1));

    pcolor(xVec,yVec,flipud(img));
    ax.YScale = 'log';
    colormap(bone);
    xlabel('Log_{10} Plume height (km)');
    ylabel('Log_{10} MER (kg/s)');
end

shading flat, axis tight, box on, set(ax, 'layer', 'top');
