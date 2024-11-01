using Base.Filesystem

version = "v1.0.0-rc3"  # the version number to update
symlink = "dev"         # the symlink to update

rm("./$version/", force=true, recursive=true) # remove the current version
cptree("./docs/build/", "./$version/")
rm("./$symlink", force=true)
run(`ln -sf $version $symlink`)
