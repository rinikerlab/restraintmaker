#clean up:
make clean;

#make doku
##configurations
sphinx-apidoc -o _source ../../restraintmaker

cp ../../examples/*ipynb ./examples

python conf.py

##execute making docu
make html
#make latex

cp ../index.html /_build/html/index.html