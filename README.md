### Introduction

libseq is a library to manipulate biological sequences (nucleotide or amino acid sequences).

###Goal

This library provides a limited set of what can be generally found in the biojava, bioperl and biopython projects, in C++. The scope is however limited to efficient and simple sequence manipulation. 

libseq originated from the [PhD project of Koen Deforche], this project was defended in February 2008. The project involved molecular evolution of HIV, with respect to resistance development against antiretroviral products.

###Features

* Nucleotide and AminoAcid sequence class based on std::vector
* Translation from Nucleotide to Aminoacid sequences
* FASTA I/O
* Full support for IUPAC ambiguity symbols
* Simple C++ interface

###Status

The [current public release](https://github.com/rega-cev/libseq/releases/tag/v0.3) (v0.3, released on July, 17 2006), contains the basic functionality, which was needed for Koen's PhD. It has been exhaustively tested in various projects small and large. You are welcome to add features you find missing.

###Documentation
* [API documentation]

####Installation

We use CMake for the build process, and tested this on GNU/Linux, MacOS and Windows (Visual Studio C++ Express).
There is an [INSTALL file] with build instructions on GNU/Linux and MacOS, which follows a standard CMake build and install process.

####Reference
* [API documentation]
* Two simple examples [NTFastaRead.C] and [AAFastaRead.C].

###GitHub
The project is hosted at [GitHub].

[PhD project of Koen Deforche]:<http://rega.kuleuven.be/cev/avd/publications/thesises>
[current public release]:<https://github.com/rega-cev/libseq/releases/tag/v0.3>
[API documentation]: <http://rega-cev.github.io/libseq/doc/doxygen-html/v0.3> 
[INSTALL file]: <https://github.com/rega-cev/libseq/blob/master/INSTALL>
[NTFastaRead.C]: <https://github.com/rega-cev/libseq/blob/master/test/src/NTFastaRead.C>
[AAFastaRead.C]: <https://github.com/rega-cev/libseq/blob/master/test/src/AAFastaRead.C> 

[GitHub]: <https://github.com/rega-cev/libseq/>
