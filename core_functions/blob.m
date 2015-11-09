function [] = blob(x,y,area,color);

hs = sqrt(area)/2;

xcorners = [x-hs, x+hs, x+hs, x-hs];
ycorners = [y-hs, y-hs, y+hs, y+hs];

%fill(xcorners,ycorners,color,'LineWidth',0.001,'EdgeColor',0.65*[1 1 1]);
fill(xcorners,ycorners,color,'LineStyle','none','BackFaceLighting','unlit');


