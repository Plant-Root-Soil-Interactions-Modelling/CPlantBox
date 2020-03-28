clear all; 
close all; 

if(~isdeployed)
  cd(fileparts(which(mfilename)));
end
plant = 'Faba'; 
path = ([plant, '/RSML']);  
results = (['results']);  
home = pwd; 
addpath(home)

cd(['../', path]);
dirpos = dir2; 
lenunit = 'mm'; 
if strcmp(lenunit,'mm') 
    factor = .1; 
elseif strcmp(lenunit, 'cm') 
    factor = 1; 
end
monocot = 0; 

nn = 1; wrong= []; 

for oo = 1:length(dirpos)
    cd(dirpos(oo).name); 
    
        filedummy = dir2; 
        for ii = 1:(size(filedummy,1))
            FileName{ii} = filedummy(ii,1).name; 
            allname{ii} = regexp(FileName{ii},'\d*','Match','once'); 
        end 
        
        indiv = sort(str2double(allname)); 
        indiv = unique(indiv); 
        indiv = indiv(~isnan(indiv)); 
        
        %%%%%%%%%%%%%%% ArchiDART%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %compute root measures
        load('params') 
        pp = general_params(params); 
        p{oo,1} = pp; 
        p{oo,1}.name = dirpos(oo).name;   %name 
        p{oo,1}.DAP = indiv;    %root system age 

        %compute CRootBoxparams
        load(strcat('DAP', num2str(indiv))); 
        %%%%%%%%%%%%%%%mine%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fid = fopen(strcat('DAP', num2str(indiv),'_rsml.txt'));
        [CR, ppa, ppp] = model_params_faba(fid, traject, factor, monocot); 
        
        CRootBox_{oo,1} = CR; 
        
        if strcmp(plant, 'maize') 
            zgh = 2; 
        else 
            zgh = size(CRootBox_{oo,1},2); 
        end
        
        for ii = 1:size(CRootBox_{oo,1},2)
            CRootBox{oo,ii} = CRootBox_{oo,1}{1,ii}; 
        end
        CRootBox{oo}.name = dirpos(oo).name; 
        
        pa{oo,1} = ppa; 
        pa{oo,1}.name = dirpos(oo).name;   %name 
        pa{oo,1}.DAP = indiv;    %root system age 

        pparams{oo,1} = ppp; 

        clear FileName
        cd ..;
end

fname = (['../', results, '/RootBox_params.mat']); 
save(fname, 'CRootBox', 'pparams'); 
fname = (['../', results, '/general_params.mat']); 
save(fname, 'p', 'pa'); 
display('done')