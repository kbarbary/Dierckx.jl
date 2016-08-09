using Compat

if is_unix()
    cd(joinpath(dirname(@__FILE__), "src", "ddierckx"))

    suffix = is_apple() ? "dylib" : "so"
    try
       run(`make FC=ifort SUFFIX=$suffix`)
    catch
       run(`make FC=gfortran SUFFIX=$suffix`)
    end
else # Windows
    # these binaries were cross-compiled from Cygwin via the following steps:
    # mkdir -p bin32 && mkdir -p bin64
    # i686-w64-mingw32-gfortran -o bin32/libddierckx.dll -O3 -shared \
    #   -static-libgfortran -static-libgcc src/ddierckx/*.f
    # x86_64-w64-mingw32-gfortran -o bin64/libddierckx.dll -O3 -shared \
    #   -static-libgfortran -static-libgcc src/ddierckx/*.f
    url = "https://cache.julialang.org/https://bintray.com/artifact/download/tkelman/generic/ddierckx.7z"
    try
        run(`curl -LO $url`)
    catch
        run(`powershell -Command "(new-object net.webclient).DownloadFile(\"$url\", \"ddierckx.7z\")"`)
    end
    run(`7z x -y ddierckx.7z`)
end
