VENV=$PWD/.venv
export R_LIBS_USER=$VENV/lib/R/library
export R_LIBS_SITE=''
mkdir -p $R_LIBS_USER
R --slave -e "if (!'renv' %in% installed.packages()[,1]) { install.packages('renv') }; renv:::restore()"
