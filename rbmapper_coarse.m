function cmap=rbmapper_coarse(cmax,cmin)

%function that returns a red/blue divergent colormap centered at zero


nrat=cmax./abs(cmin);


bmap=[0.1412         0    0.8471
    0.0941    0.1098    0.9686
    0.1569    0.3412    1.0000
    0.2392    0.5294    1.0000
    0.3373    0.6902    1.0000
    0.4588    0.8275    1.0000
    0.6000    0.9176    1.0000
    0.7373    0.9765    1.0000
    0.9176    1.0000    1.0000];

rmap=[1.0000    1.0000    0.9176
    1.0000    0.9451    0.7373
    1.0000    0.8392    0.6000
    1.0000    0.6745    0.4588
    1.0000    0.4706    0.3373
    1.0000    0.2392    0.2392
    0.9686    0.1529    0.2078
    0.8471    0.0824    0.1843
    0.6471         0    0.1294];



if abs(cmin)<abs(cmax); %then we want more red colors
    nr=floor(5.*nrat); %number of red colors
    nb=5;
else
    nb=floor(5.*1/nrat);
    nr=5;
end


%the highest color will always just saturate in this scheme
lr=linspace(1,size(rmap,1),size(rmap,1));
lr2=linspace(1,size(rmap,1),nr);
lb=linspace(1,size(bmap,1),size(bmap,1));
lb2=linspace(1,size(bmap,1),nb);
r1=interp1(lr,rmap(:,1),lr2);
r2=interp1(lr,rmap(:,2),lr2);
r3=interp1(lr,rmap(:,3),lr2);

rmapf=cat(2,r1',r2',r3');

b1=interp1(lb,bmap(:,1),lb2);
b2=interp1(lb,bmap(:,2),lb2);
b3=interp1(lb,bmap(:,3),lb2);

bmapf=cat(2,b1',b2',b3');

cmap=[bmapf;rmapf];
