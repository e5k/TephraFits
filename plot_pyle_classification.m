A = imread('style_Pyle89.jpg');
A = A(:,:,1);

xVec    = logspace(-1,2,size(A,2));
yVec    = linspace(0,4,size(A,1));

figure;
pcolor(xVec,yVec,flipud(A));
shading flat, axis tight, box on, set(ax, 'layer', 'top')
ax = gca; ax.XScale = 'log';
colormap(gray)
xlabel('Thickness half-distance b_T (km)');
ylabel('Half-distanceratio b_C/b_T');


B = imread('style_BC13.jpg');
B = B(:,:,1);

xVec    = logspace(0,log10(400),size(B,2));
yVec    = logspace(-2,2,size(B,1));
[X,Y]   = meshgrid(xVec,yVec);

figure;
pcolor(xVec,yVec,flipud(B));
shading flat, axis tight, box on
ax = gca; ax.XScale = 'log'; ax.YScale = 'log'; set(ax, 'layer', 'top');
colormap(bone)



xlabel('Log_{10} \lambda_{th}');
ylabel('Log_{10} \lambda_{MC}/\lambda_{th}');
