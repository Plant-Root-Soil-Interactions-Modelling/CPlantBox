function [CRootBox, pa, pparams] = model_params_new( fid, traject, factor, monocot)
%model_params computes model input parameters for the CRootBox model 
%part one computes general root measures
%part two computes model input parameters 

        testRL = 0; 
        testlat = 0; 

        traject.Ord(traject.Ord>3) = 3; 
        ord= unique(traject.Ord); 
        for ii = 1 : length(ord) 
            CRootBox{1,ii}.order = ord(ii); 
            CRootBox{1,ii}.sigma(1) = nanmean(traject.Mean_Curv(traject.Ord==ii)); 
            CRootBox{1,ii}.sigma(2) = nanmean(traject.SD_Curv(traject.Ord==ii));
            
            CRootBox{1,ii}.theta(1) = nanmean(traject.Branching_Angle(traject.Ord==ii)); 
            CRootBox{1,ii}.theta(2) = nanstd(traject.Branching_Angle(traject.Ord==ii));
            if CRootBox{1,ii}.theta(1) > 90
                CRootBox{1,ii}.theta(1) = 180-CRootBox{1,ii}.theta(1); 
            end
        end


        A = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','HeaderLines',1,'TreatAsEmpty',{'NA','N/A'});
        fclose(fid);
        A_root = A{1,2};
        A_order = A{1,3}; 
        A_parent = A{1,4};
        A_age = A{1,5}; 
        A_bran = A{1,6};
        A_apic = A{1,7};
        A_xf = A{1,8}*factor; 
        A_yf = A{1,9}*factor; 
        A_zf = -1*(A{1,10})*factor; 
        A_xg = A{1,11}*factor; 
        A_yg = A{1,12}*factor; 
        A_zg = -1*(A{1,13})*factor; 
        A_a = A{1,14}*factor; 
        A_len = A{1,15}*factor;
        ord = unique(A_order); 
        ord = sort(ord); 
        
        Ax = []; Ay = []; Az = []; 
        idx = find(A_bran==1); 
        for ii = 1:size(idx,1)-1
            Ax = [Ax; diff(A_xf(idx(ii):idx(ii+1)-1))]; 
            Ay = [Ay; diff(A_yf(idx(ii):idx(ii+1)-1))]; 
            Az = [Az; diff(A_zf(idx(ii):idx(ii+1)-1))]; 
        end
        
        %compute general params 
        %1) tot length 
        for ii = 1: length(Ax) 
            len(ii,1) = norm([Ax(ii), Ay(ii), Az(ii)]); 
            surf(ii,1) = len(ii)*2*pi*A_a(ii); 
            vol(ii,1) = pi*A_a(ii)^2*len(ii); 
        end
        [~, V] = convhull([A_xf, A_yf, A_zf]); 
        

        pa.RL = sum(len);   %total root length (cm)
        pa.RSA = sum(surf);     %total root surface area (cm²)
        pa.RV = sum(vol);  %total root volume (cm³)
        pa.zmax = max(A_zf)-min(A_zf); %does not make sense- max vertical spread (cm)
        pa.rmax = max(max(A_xf),max(A_yf))-min(min(A_xf),min(A_yf)); %does not make sense- max horizontal spread (cm)
        pa.conv = V; %does not make sense - convex hull (cm³)
        pa.NR = numel(unique(A_root)); %number of root tips 
        pa.RLD = sum(A_len)/V;   %root length normalized by the convex hull
        pa.HMD = (pi*sum(A_len)/V)^(-0.5); %half mean distance between roots (cm)
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %compute CRootBoxparams
        bridx = find(A_bran==1); %index of first segments of all roots 
        apidx = find(A_apic==1); %index of last segments of all roots 
        

        for ii = 1:length(ord) 
            dummy = A_root(A_order==ii); %roots of order ii 
            roots = (sort(unique(dummy))); %unique root numbers of order ii
            ln_ = []; 
            la_ = []; 
            lb_ = []; 
            latnum_ = []; 
            rootlength = []; 
            elong = []; 
            

            for jj = 1:length(roots) %all roots of this order
                rootidx = find(A_root == roots(jj)); %indices of all segments of root jj  
                latidx = bridx(A_parent(bridx) == roots(jj));   %index of first segments of those roots that have root ii as a parent
                coord = [A_xf(latidx), A_yf(latidx), A_zf(latidx)]; %coordinates of first child root segments 
                ap_coord = [A_xg(latidx), A_yg(latidx), A_zg(latidx)]; 
                latnum_(jj) = size(coord,1); 
                
                %elongation rate: 
                indum = find(A_root == roots(jj), 1, 'first'); %index of first segment
                apdum = find(A_root == roots(jj), 1, 'last'); %index of last segment 
                ini_age = A_age(indum); 
                ap_age = A_age(apdum); 
                time_ = ap_age-ini_age; 


                for kk = [indum:apdum]
                    coord_ini = [A_xf(kk), A_yf(kk), A_zf(kk)];
                    coord_apic = [A_xg(kk), A_yg(kk), A_zg(kk)];
                    len_(kk) = norm(coord_apic-coord_ini); 
                    if len_(kk)>2
                        display('is that correct?')
                    end
                end


                length_ = sum(len_); 
                if time_== 0
                    elong(jj,1) = nan; 
                else
                    elong(jj,1) = length_/time_; %elongation rate 
                end

                
                %maximum root length: sum up all root segement lengths 
                rootlength(jj,1) = length_; 
                testRL = testRL+length_; 
                
                if (size(coord,1)>1) %if the root has children
                    ln_diff = diff(coord); 
                    for ll = 1: size(ln_diff,1)
                        lndummy(ll,1) = norm(abs(ln_diff(ll,:)));
                    end
                    ln_ = [ln_; lndummy]; %put lns for this order for all roots below each other 
                    clear lndummy; 
                end


                coord_ini = [A_xf(indum), A_yf(indum), A_zf(indum)]; %initial coordinates of each root 
                coord_apic = [A_xg(apdum), A_yg(apdum), A_zg(apdum)]; %last coordinates of each root 
                if size(coord,1)>0
                    Vap = bsxfun(@minus, coord, coord_apic);
                    Dap = sqrt(sum(Vap.^2, 2)); 
                    la_ = [la_; min(Dap)];
                    
                    Vini = bsxfun(@minus, ap_coord, coord_ini);
                    Dini = sqrt(sum(Vini.^2, 2)); 
                    lb_ = [lb_; min(Dini)];
                end
                clear Vini Vap Dini Dap coord apcoord coord_apic coord_ini rootdum len_ rootidx latidx
                

                if latnum_(jj) == 0
                    ln1_(jj,1) = 0; 
                else
                    ln1_(jj,1) = (rootlength(jj)-nanmean(la_)-nanmean(lb_))/latnum_(jj); 
                end

                
            end
            
            elong(elong==inf) = nan; 
            CRootBox{1,ii}.r(1) = nanmean(elong); 
            CRootBox{1,ii}.r(2) = nanstd(elong);   
            CRootBox{1,ii}.lmax(1) = nanmean(rootlength); 
            CRootBox{1,ii}.lmax(2) = nanstd(rootlength);  
            
            if or(ii>2, isempty(la_))
                CRootBox{1,ii}.la(1) = nanmean(rootlength)/2; 
                CRootBox{1,ii}.la(2) = nanstd(rootlength)/2; 
                CRootBox{1,ii}.lb(1) = nanmean(rootlength)/2;
                CRootBox{1,ii}.lb(2) = nanstd(rootlength)/2;
            end

            

           if and(ii<=2,~isempty(la_))
                CRootBox{1,ii}.la(1) = nanmean(la_);  
                CRootBox{1,ii}.la(2) = nanstd(la_); 
                CRootBox{1,ii}.lb(1) = nanmean(lb_);  
                CRootBox{1,ii}.lb(2) = nanstd(lb_); 
                
                if (nanmean(la_)+nanmean(lb_))>CRootBox{1,ii}.lmax(1)
                    display('sth wrong')
                    CRootBox{1,ii}.la(1) = CRootBox{1,ii}.lmax(1)-nanmean(lb_); 
                end
           end
           
           if isempty(ln_)
                CRootBox{1,ii}.ln = [0, 0];
                CRootBox{1,ii}.nob = [0, 0];
           else
                CRootBox{1,ii}.ln(1) = nanmean(ln_); 
                CRootBox{1,ii}.ln(2) = nanstd(ln_); 
                CRootBox{1,ii}.nob(1) = nanmean(latnum_);
                CRootBox{1,ii}.nob(2) = nanstd(latnum_); 
           end
           if and(ii<3, CRootBox{1,ii}.ln(1)==0)
               CRootBox{1,ii}.ln(1) = 0.01; 
           end
           
           
           testlat = testlat+sum(latnum_); 
        end
       
        
        for ii = 1: length(ord) 
            CRootBox{ii}.a(1) = nanmean(A_a(A_order==ii)); 
            CRootBox{ii}.a(2) = nanstd(A_a(A_order==ii)); 
        end
       CRootBox{1}.totlen = sum(A_len); 
       

        
        %compute pparams 
        if (monocot)
            [~, idx] = unique(A_root); 
            basdummy = idx(A_order(idx)==1); %all 1st order roots
            if numel(basdummy)>1
                pparams.maxB = sum(traject.Ord==1)-1; %basal roots = roots other than taproot that emerge from the base
                pparams.firstB = min(A_age(basdummy)); 
                pparams.delayB(1) = nanmean(diff(A_age(basdummy))); 
                pparams.delayB(2) = nanstd(diff(A_age(basdummy))); 
            else 
                pparams.maxB = 1; 
            end
        else 
            pparams = 0; 
        end
end

