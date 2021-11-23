# liverPFOA
Code for project liver toxicological research

Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>.
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.
 
 plotHicGenome: the tools are used for displaying hic signal of whole genome
===========================================================================
Installation
----------------------------------------------------------------------------

Python package and command line interface for displaying Hic signals of assembly genome
was tested with python 2.7 on linux. You can install it using pip or through source codes.

* Dependent pakages<br/>
matplotlib-2.2.3<br/>
pandas-0.24.0<br/>
numpy-1.16.4<br/>
scipy-1.2.0<br/>

* Install by pip<br/>
   pip install plotHicGenome --user<br/>
* Install through raw codes<br/>
  git clone  https://github.com/chenjunhui/plotHicGenome<br/>
  cd  plotHicGenome<br/>
  python  setup.py  install   [--prefix=/user/direction]<br/>

Notably: you'd better install it under virtual environment in case  system conflict.<br/>


USAGE
==============================================================================================
The package includes two sub-commands {Hicproc, juicer}, and you can choose corresponding module basing
on result of the your hic matrix.

```Bash
/user/direction/bin/plotHicGenome  --help
```
