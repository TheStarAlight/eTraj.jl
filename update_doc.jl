using Base.Filesystem

version = "v1.0.0-rc2"    # the version number to update

rm("./$version/", recursive=true) # remove the current version
cptree("./docs/build/", "./$version/")
