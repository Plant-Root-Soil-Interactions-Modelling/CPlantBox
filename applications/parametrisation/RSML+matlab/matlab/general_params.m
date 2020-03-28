function [ p ] = general_params(params)
%params computed with ArchiDART 

        
        p.RL = params.TRL(end)/10;   %total root length (cm)
        p.RSA = params.Stot(end)/10^2;     %total root surface area (cm²)
        p.RV = params.Vtot(end)/10^3;  %total root volume (cm³)
        p.zmax = params.Height(end)/10; %does not make sense- max vertical spread (cm)
        p.rmax = params.Width(end)/10; %does not make sense- max horizontal spread (cm)
        p.NR = params.TN1R(end) + params.TNLR(end); %number of root tips 

end

