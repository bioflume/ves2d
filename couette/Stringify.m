function fname = Stringify(fileSpec)
  
  if ( isstr(fileSpec) )
    fname = fileSpec;
  else

    fname = fileSpec.prefix;
    
    if ( ~isempty(fileSpec.perturb) )
      fname  = [fname '_' num2str(100*fileSpec.perturb)];
    end
    
    fname  = [fname '_' num2str(fileSpec.volFrac)];
    
    if ( ~isempty(fileSpec.reducedArea) )
      fname  = [fname '_' num2str(100*fileSpec.reducedArea)];
    end

    if ( ~isempty(fileSpec.viscCont) )
      fname  = [fname '_' num2str(fileSpec.viscCont)];
    end

    fname = [fname '_' num2str(fileSpec.T) '_' num2str(1000 * fileSpec.ts) ...
             '_' num2str(fileSpec.n)];
  end  