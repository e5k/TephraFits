xV = 785284;
yV = 9924355;

a = xlsread('Paper/L5.xls');

DWi = a(:,4) == 1;
XWi = a(:,5) == 1;

thicknessDW = a(DWi,3);
thicknessXW = a(XWi,3);

distanceDW = sqrt( (xV-a(DWi,1)).^2 + (yV-a(DWi,2)).^2 )./1e3;
distanceXW = sqrt( (xV-a(XWi,1)).^2 + (yV-a(XWi,2)).^2 )./1e3;

[distanceDW, DWi] = sortrows(distanceDW);
[distanceXW, XWi] = sortrows(distanceXW);

thicknessDW = thicknessDW(DWi);
thicknessXW = thicknessXW(XWi);