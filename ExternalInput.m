function Iext = ExternalInput(model,nbpop,dir)
    
    Param = sprintf(['../%s/Parameters/%dpop/%s/Param.txt'],model,nbpop,dir) ;

    fid = fopen( Param ) ;
    if fid ~= -1
        Iext = fgetl ( fid ) ;
        fclose ( fid );
        Iext = str2num(Iext(5:end)) ;
    else 
        fprintf('Param.txt not found\n')
        Iext = NaN(1,nbpop) ;
    end
    

end

