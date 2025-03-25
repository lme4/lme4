#!/bin/sh


echo "Checking R package dependencies..."
Rscript -e 'if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes"); remotes::install_deps(dependencies=TRUE)'

echo "Building the source tarball..."
R CMD build .

tarball=$(ls *.tar.gz | head -n 1)


if [ -z "$tarball" ]; then
    echo "Error: No tarball was created."
    exit 1
fi

echo "Tarball created: $tarball"

echo "Running R CMD check on the tarball..."
R CMD check "$tarball"

check_dir="$(basename "$tarball" .tar.gz).Rcheck"
if [ -d "$check_dir" ]; then
    echo "Displaying check results..."
    cat "$check_dir/00check.log"
else
    echo "Error: Check directory not found."
    exit 1
fi

