@unix_only begin
    cd(joinpath(dirname(@__FILE__), "src", "ddierckx"))

    suffix = @osx? "dylib" : "so"
    run(`make FC=gfortran SUFFIX=$suffix`)
end

@windows_only begin
    # these binaries were cross-compiled from Cygwin via the following steps:
    # mkdir -p bin32 && mkdir -p bin64
    # i686-w64-mingw32-gfortran -o bin32/libddierckx.dll -O3 -shared \
    #   -static-libgfortran -static-libgcc src/ddierckx/*.f
    # x86_64-w64-mingw32-gfortran -o bin64/libddierckx.dll -O3 -shared \
    #   -static-libgfortran -static-libgcc src/ddierckx/*.f
    run(`curl -LO https://cache.e.ip.saba.us/https://bintray.com/artifact/download/tkelman/generic/ddierckx.7z`)
    run(`7z x -y ddierckx.7z`)
end
