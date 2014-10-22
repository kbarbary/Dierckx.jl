const msg = "No OSX or Windows Makefiles!"

cd(joinpath(dirname(@__FILE__), "src", "ddierckx"))
@linux? run(`make`): error(msg)
