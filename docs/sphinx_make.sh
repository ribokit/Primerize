set -e

# start at /
# e.g. run as docs/sphinx_make.sh
cd docs

mkdir -p source/_theme
cd source/_theme
if [ ! -d ribokit-Sphinx-theme ]; then
    git clone https://github.com/ribokit/ribokit-Sphinx-theme/
    cd ribokit-Sphinx-theme
else
    cd ribokit-Sphinx-theme
    git fetch
fi
git checkout 90e4deea66cb17a6f7dcd7ccb6f619aa8e4e6279
cd ../../..

make clean
make html

cd build/html/_static/
rm basic.css pygments.css
rm file.png minus.png plus.png

cd ../../../../

# switch to gh-pages
git checkout gh-pages
git pull
cp -r docs/build/html/* ./
git add -A
git commit -m "$(date)"
git push

# switch to master
git checkout master
