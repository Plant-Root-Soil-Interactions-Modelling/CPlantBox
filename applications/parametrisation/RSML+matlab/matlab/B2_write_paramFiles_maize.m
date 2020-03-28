clear all; clear all; 
close all; 

if(~isdeployed)
  cd(fileparts(which(mfilename)));
end
plant = 'Maize'; 
path = ([plant, '/RSML']);  
results = (['results']);  
home = pwd; 
addpath(home)

path = ([plant,'/results']); 
namepath = ([plant,'/RSML']); 
txtresults = (['../', plant,'/paramsFiles']); 


cd(['../', namepath]) 
nmes = dir2; 
write = 1; 
factor = 1; %cm
fac = 3; %reduce std, otherwise lin growth not possible

cd(['../../', path]) ; 
load('RootBox_params'); 
cd(home); 


monocot = 1; 
cols = [255, 0, 0; 0, 255, 0; 0,0,255];  
r(:,1) = mean(cellfun(@(x) x.r(1), CRootBox),1).'; 
lb(:,1) = mean(cellfun(@(x) x.lb(1), CRootBox),1).'; 
la(:,1) = mean(cellfun(@(x) x.la(1), CRootBox),1).'; 
ln(:,1) = mean(cellfun(@(x) x.ln(1), CRootBox),1).'; 
nob(:,1) = mean(cellfun(@(x) x.nob(1), CRootBox),1).'; 
lmax(:,1) = mean(cellfun(@(x) x.lmax(1), CRootBox),1).';
a(:,1) = mean(cellfun(@(x) x.a(1), CRootBox),1).'; 
theta(:,1) = mean(cellfun(@(x) x.theta(1), CRootBox),1).'; 
sigma(:,1) = mean(cellfun(@(x) x.sigma(1), CRootBox),1).'; 

firstB(1,1) = mean(cellfun(@(x) x.firstB(1), pparams),1).'; 
maxB(1,1) = mean(cellfun(@(x) x.maxB(1), pparams),1).'; 
delayB(1,1) = mean(cellfun(@(x) x.delayB(1), pparams),1).'; 

r(:,2) = sqrt(mean(cellfun(@(x) x.r(2), CRootBox),1).^2).'/100; 
lb(:,2) = sqrt(mean(cellfun(@(x) x.lb(2), CRootBox),1).^2).'; 
la(:,2) = sqrt(mean(cellfun(@(x) x.la(2), CRootBox),1).^2).'; 
ln(:,2) = sqrt(mean(cellfun(@(x) x.ln(2), CRootBox),1).^2).'; 
nob(:,2) = sqrt(mean(cellfun(@(x) x.nob(2), CRootBox),1).^2).'; 
lmax(:,2) = sqrt(mean(cellfun(@(x) x.lmax(2), CRootBox),1).^2).'; 
a(:,2) = sqrt(mean(cellfun(@(x) x.a(2), CRootBox),1).^2).'; 
theta(:,2) = sqrt(mean(cellfun(@(x) x.theta(2), CRootBox),1).^2).'; 

firstB(:,2) = std(cellfun(@(x) x.firstB(1), pparams),1).'; 
maxB(:,2) = std(cellfun(@(x) x.maxB(1), pparams),1).'; 
delayB(:,2) = sqrt(mean(cellfun(@(x) x.delayB(2), pparams),1).^2).'; 

for jj = 1:size(CRootBox, 2)
    if jj == 1
        p{jj}.name = 'basal roots'; 
    else 
        p{jj}.name = 'lateral'; 
    end

    p{jj}.r = r(jj,:)*factor;  
    p{jj}.lb = lb(jj,:)*factor; 
    p{jj}.la = la(jj,:)*factor; 
    p{jj}.ln = ln(jj,:)*factor; 
    p{jj}.nob = nob(jj,:); 
    p{jj}.lmax = lmax(jj,:)*factor; 
    p{jj}.a = a(jj,:)*factor;
    p{jj}.gf = 2;
    p{jj}.theta = theta(jj,:)/180*pi;


    dummy = sigma(jj,:)/180*pi;
    if jj == 1
        p{jj}.theta = [1.3963	0.105941]; 
        p{jj}.tropism = [1,1.5,dummy]; 
    else 
        p{jj}.tropism = [2,2,dummy]; 
    end
    p{jj}.color = cols(jj,:); 
    p{jj}.dx= .25;
    if jj<(size(CRootBox, 2))
        p{jj}.successor = jj+1;
    else
        p{jj}.successor = 0;
    end
    
    if(monocot) 
        pp{jj}.firstB = firstB;
        pp{jj}.maxB = ceil(maxB); 
        pp{jj}.delayB = delayB;
    end
end

fname = 'MWeber_Maize'; 
if (write)
    write_params2txt(p, pp, fname, txtresults, monocot)
end
display('done')

