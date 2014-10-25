cd(joinpath(dirname(@__FILE__), "src", "ddierckx"))

suffix = @linux? "so" : @osx? "dylib" : error("only linux and osx supported")
run(`make FC=gfortran SUFFIX=$suffix`) 
