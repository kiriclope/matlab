function WriteParam(model,nbpop,dir,Iext,J,Rates) ;
        path = sprintf(['../%s/Parameters/%dpop/%s'],model,nbpop,dir) ; 
        try 
            mkdir(path) ; 
        end 

        fprintf('Writing Parameters to : ') ; 
        filename = sprintf(['%s/Param.txt'], path) ; 
        Jparam = sprintf(['%s/Jparam.txt'], path) ; 
        disp(filename) 
        file = fopen(filename,'w') ; 
        Jfile = fopen(Jparam,'w') ; 
        
        fprintf(file, 'Iext ') ; 
        fprintf(file, '%.3f ', Iext) ; 
        fprintf(file, '\n') ; 
        
        fprintf(file, 'Connectivity Matrix \n') ; 
        for i=1:nbpop 
            for j=1:nbpop 
                fprintf(file, '%.3f ', J(i,j)) ; 
                fprintf(Jfile, '%.3f ', J(i,j)) ; 
            end 
            fprintf(file, '\n') ; 
            fprintf(Jfile, '\n') ; 
        end 
        fclose(Jfile) ; 
        
        fprintf(file,'MF Rates ') ; 
        fprintf(file,'%.3f | ', Rates) ; 
        fprintf(file,'\n') ; 
        
        fclose(file) ; 
end