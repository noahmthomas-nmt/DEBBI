# Initial submission
This is my first submission

## R CMD check results
### Platform:	Ubuntu Linux 20.04.1 LTS, R-release, GCC
There were no ERRORs or WARNINGs. There was one note.
  >Possibly misspelled words in DESCRIPTION:
    >DEMAP (10:167)
    name of a function
    >DEMCMC (10:196)
    name of a function
    >variational (10:208)
    variational inference is one of the methods in the package

### Platform:	Fedora Linux, R-devel, clang, gfortran
There were no ERRORs or WARNINGs. There was one note.
Same as above

### Platform:	Windows Server 2022, R-devel, 64 bit
There were no ERRORs or WARNINGs. There was two note.
Same as above + 
> * checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'

The second note might be due to a bug/crash in MiKTeX, as is noted in R-hub issue #503, which can be found at https://github.com/r-hub/rhub/issues/503.
