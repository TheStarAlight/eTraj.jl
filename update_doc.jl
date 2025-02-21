using Base.Filesystem
using Base.Sys

@assert length(ARGS)==2

version = ARGS[1]  # the version number to update
symlink = ARGS[2]  # the symlink to update

rm("./$version/", force=true, recursive=true) # remove the current version
cptree("./docs/build/", "./$version/")
rm("./$symlink", force=true)
if iswindows()
    println("Run the following command:\ncmd -c mklink /D $symlink $version")
else
    run(`ln -sf $version $symlink`)
end
