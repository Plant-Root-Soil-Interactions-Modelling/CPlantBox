function p = completeParameters(p,df)
% completeParameters: adds missing default root system parameters
%
% Checks all parameters and copies it to global parameters. Each root of
% type 'i' is described by the following parameters:
% p{i}.lb = [mean, std]     length of the basal zone (default = [0,0])
% p{i}.la = [mean, std]     length of the apical zone (default = [10,0])
% p{i}.ln = [mean, std]     inter-branch distance (default = [0,0])
% p{i}.nob = [mean, std]    number of branches (default = [0,0])
% p{i}.r = [mean, std]      initial growth speed (cm/day) (default = [1,0])
% p{i}.a = [mean, std]      root radius (constant along axis) 
%                           (default = [1e-2,0]) 
% p{i}.color = [r,g,b]      color of this root type 
%                           (default = [150/255,150/255,50/255]) 
% p{i}.tropism = [t,N,s]    t is the type (0=Plagio-, 1=Gravi-, 2=Exo-,
%                           3=Chemotropism), N is the number of trials, s is the maximal angular
%                           deviation in root heading (default=[1 1 pi/20])                                                                                                
% p{i}.dx                   axial resolution (default = 0.1)
% p{i}.successor = [I,P]    I indices of successive types, P probability
%                           (default = [])          
% p{i}.theta = [mean, std]  insertion angle (default = [70/180*pi,0])                                                       
% p{i}.rlt = [mean, std]    root life time  (default = [inf,0])
% p{i}.gf                   type of the growth function (default = 1)
% p{i}.name                 name of the root type (default = 'unknown')
% p{i}.ief                  impede elongation function (default = @(x) 1)
% p{i}.ibf                  impede branching function (default = @(x) 1)
% p{i}.irlf                 impede root life (default = @(x) 1)
% p{i}.type                 index within list
%
% Parameters:
% p                         a list of parameter sets (for each root type)
% (df)                      signed distance function of container geometry
%                           (default = @(x) x(3))                  
%
% p                         parameter set including default values
%
% See also: applyRules
%
% Copyright 2012-2014 Daniel Leitner. See license.txt for details.
%

global parameters;
global geometry;

if nargin<2
    if isempty(geometry) % if not already set
        geometry = @(x) x(:,3);
    end
else 
    geometry = df;
end

%
% default parameters 
%
defaultoptions=struct('lb',0,'la',10,'ln',0,'nob',0,'r',1,'a',0.01,...
    'color',[150/255,150/255,50/255],'tropism',[1 1 pi/20], 'dx', 0.1, ...
    'successor',[],'theta',70/180*pi,'rlt',inf,'gf',1,'name','unknown',...
    'sef',@(x) 1,'sbpf',@(x) 1,'saf',@(x) 1,'type', 0);

for i = 1 : length(p)
    
    %
    % add missing parameters to all root types
    %
    tags = fieldnames(defaultoptions);
    for j=1:length(tags)
        if(~isfield(p{i},tags{j}))
            p{i}.(tags{j})=defaultoptions.(tags{j}); % add default
        end
        if isreal(p{i}.(tags{j})) && length(p{i}.(tags{j}))==1 && ...
                ~strcmp(tags{j},'successor');
            p{i}.(tags{j})=[p{i}.(tags{j}),0]; % add zero std
        end
    end
    
    % check tropism parameter length
    if  length(p{i}.tropism)~=3
        warning('completeParameters: Wrong number of tropism arguments');
    end

    % check successors
    if ~isempty(p{i}.successor)
        if size(p{i}.successor,2)==1 % if only one number 
            p{i}.successor(1,2)=1; % add probability 1 
        end
        if sum(p{i}.successor(:,2))>1
            warning('completeParameters: Successors probabilities > 1');
        end
    else
        if p{i}.nob>0
            p{i}.successor(1,1)=1; % set a successor
            p{i}.successor(1,2)=0; % with probability 0
        end
    end
    
    % warn about unknown options
    if(length(tags)~=length(fieldnames(p{i})))
        warning('completeParameters: Unknown options found');
    end       
    
    p{i}.type = i; % add index of type                
    
end

parameters=p;
