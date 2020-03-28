clear all; clear all; 
close all; 

if(~isdeployed)
  cd(fileparts(which(mfilename)));
end
plant = 'Faba'; 
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

monocot = 0; 
cols = [255, 0, 0; 0, 255, 0; 0,0,255]; 
r(:,1) = mean(cellfun(@(x) x.r(1), CRootBox)).'; 
lb(:,1) = mean(cellfun(@(x) x.lb(1), CRootBox)).'; 
la(:,1) = mean(cellfun(@(x) x.la(1), CRootBox)).'; 
ln(:,1) = mean(cellfun(@(x) x.ln(1), CRootBox)).'; 
nob(:,1) = mean(cellfun(@(x) x.nob(1), CRootBox)).'; 
lmax(:,1) = mean(cellfun(@(x) x.lmax(1), CRootBox)).'; 
a(:,1) = mean(cellfun(@(x) x.a(1), CRootBox)).'; 
theta(:,1) = mean(cellfun(@(x) x.theta(1), CRootBox)).'; 
sigma(:,1) = mean(cellfun(@(x) x.sigma(1), CRootBox)).'; 

r(:,2) = sqrt(mean((cellfun(@(x) x.r(2), CRootBox)).^2)).'/fac; 
lb(:,2) = sqrt(mean((cellfun(@(x) x.lb(2), CRootBox)).^2)).'; 
la(:,2) = sqrt(mean((cellfun(@(x) x.la(2), CRootBox)).^2)).'; 
ln(:,2) = sqrt(mean((cellfun(@(x) x.ln(2), CRootBox)).^2)).'; 
lmax(:,2) = sqrt(mean((cellfun(@(x) x.lmax(2), CRootBox)).^2)).'; 
nob(:,2) = sqrt(mean((cellfun(@(x) x.nob(2), CRootBox)).^2)).';
a(:,2) = sqrt(mean((cellfun(@(x) x.a(2), CRootBox)).^2)).'; 
theta(:,2) = sqrt(mean((cellfun(@(x) x.theta(2), CRootBox)).^2)).'; 

r2 = std(cellfun(@(x) x.r(1), CRootBox)).'; 
lb2 = std(cellfun(@(x) x.lb(1), CRootBox)).'; 
la2 = std(cellfun(@(x) x.la(1), CRootBox)).'; 
nob2 = std(cellfun(@(x) x.nob(1), CRootBox)).'; 
theta2 = std(cellfun(@(x) x.theta(1), CRootBox)).'; 


%we do not want root order 3 - distribute branches from 3 to 2 - increase
%nob, decrease ln --> nob1/nob2 = ln2/ln1
nob_ = nob(1,1) + nob(2,1)*nob(1,1); 
fact = nob(1,1)/(nob(2,1)*nob(1,1)); 
ln_ = ln(2,1)/fact; 
nob(1,1) = nob_; 
ln(1,1) = ln_; 


%test tot length 
for jj = 1:size(CRootBox, 2)-1
    if jj == 1
        p{jj}.name = 'taproot'; 
    else 
        p{jj}.name = 'lateral'; 
    end


    p{jj}.r = r(jj,:)*factor; 
    p{jj}.lb = lb(jj,:)*factor; 
    p{jj}.la = la(jj,:)*factor; 
    p{jj}.ln = ln(jj,:)*factor/1.6; 
    p{jj}.nob = nob(jj,:); 
    p{jj}.lmax = lmax(jj,:)*factor; 
    p{jj}.a = a(jj,:)*factor;
    if jj == 1
        p{jj}.gf = 1;
    else 
        p{jj}.gf = 2;
    end
    p{jj}.theta = theta(jj,:)/180*pi;

    if jj == 1
        p{jj}.r(2) = r2(jj)*factor;  
        p{jj}.lb(2) = lb2(jj)*factor; 
        p{jj}.la(2) = la2(jj)*factor; 
        p{jj}.nob(2) = nob2(jj); 
        p{jj}.theta(2) = theta2(jj)/180*pi;
    end


    dummy = sigma(jj,:)/180*pi;
    if strcmp(plant, 'Faba')
        if jj == 1
            p{jj}.theta = [0, 0];
            p{jj}.tropism = [1,2,dummy]; 
        elseif jj==2
            p{jj}.tropism = [2,1.5,dummy]; 
        else
            p{jj}.tropism = [2,0.5,dummy]; 
        end
    else 
        if jj == 1
            p{jj}.tropism = [1	0.263933,dummy]; 
        elseif jj==2
            p{jj}.tropism = [2	0.348788,dummy]; 
        end
    end 
    p{jj}.color = cols(jj,:); 
    p{jj}.dx= .1*factor;
    if jj<(size(CRootBox, 2)-1)
        p{jj}.successor = jj+1;
    else
        p{jj}.successor = 0;
    end
end

fname = strcat('MWeber_',plant); 
if (write)
    write_params2txt(p, pparams, fname, txtresults, monocot)
end
display('done')






