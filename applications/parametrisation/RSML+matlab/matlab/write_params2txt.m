% writes the root system parameters in a way CRootBox can read it
function write_params2txt(p, pp, fname, path, monocot)

p = completeParameters(p);
fid = fopen([path '/' strcat(fname, '.rparam')],'w');

for i  = 1 : length(p);
    p1=p{i};
    fprintf(fid,['# Parameter set for type\n']);
    fprintf(fid,'type\t%g\n',p1.type);
    name = p1.name; % delete spaces
    name(name==' ')=[];    
    fprintf(fid,'name\t%s\n',name);
    fprintf(fid,'lb\t%g\t%g\n',p1.lb(1),p1.lb(2));
    fprintf(fid,'la\t%g\t%g\n',p1.la(1),p1.la(2));
    fprintf(fid,'ln\t%g\t%g\n',p1.ln(1),p1.ln(2));
%     fprintf(fid,'nob\t%g\t%g\n',p1.nob(1),p1.nob(2));
    fprintf(fid,'lmax\t%g\t%g\n',p1.lmax(1),p1.lmax(2));
    fprintf(fid,'r\t%g\t%g\n',p1.r(1),p1.r(2));
    fprintf(fid,'a\t%g\t%g0\n',p1.a(1),p1.a(2));
    fprintf(fid,'color\t%g\t%g\t%g\n',p1.color(1),p1.color(2),p1.color(3));
    fprintf(fid,'tropism\t%g\t%g\t%g\n',p1.tropism(1),p1.tropism(2),p1.tropism(3));
    fprintf(fid,'dx\t%g\n',p1.dx(1));
    fprintf(fid,'successors');
    fprintf(fid,'\t%g',size(p1.successor,1));
    for ii = 1 : size(p1.successor,1)
        fprintf(fid,'\t%g',p1.successor(ii,1));
    end
    fprintf(fid,'\n');
    fprintf(fid,'successorP');
    fprintf(fid,'\t%g',size(p1.successor,1));
    for ii = 1 : size(p1.successor,1)
        fprintf(fid,'\t%g',p1.successor(ii,2));
    end
    fprintf(fid,'\n');
    fprintf(fid,'theta\t%g\t%g\n',p1.theta(1),p1.theta(2));
    if isinf(p1.rlt(1))
        fprintf(fid,'rlt\t%g\t%g\n',1e9,0); % replace Inf by large number
    else
        fprintf(fid,'rlt\t%g\t%g\n',p1.rlt(1),p1.rlt(2));
    end
    fprintf(fid,'gf\t%g\n',p1.gf(1));
end
fclose(fid);


if (monocot) 
    for i  = 1 : length(pp);
        p2=pp{i};
        fid = fopen([path '/' strcat(fname, '.pparam')],'w');

        fprintf(fid,'plantingdepth\t%g\n',0);
        fprintf(fid,'firstB\t%g\n',p2.firstB(1));
        fprintf(fid,'delayB\t%g\n',p2.delayB(1));
        fprintf(fid,'maxB\t%g\n',p2.maxB(1));
        fprintf(fid,'nC\t%g\n',0);
        fprintf(fid,'firstSB\t%g\n',1000);
        fprintf(fid,'delaySB\t%g\n',1000);
        fprintf(fid,'delayRC\t%g\n',1000);
        fprintf(fid,'nz\t%g\n',0);
    end
    fclose(fid);

end

end