### Introduction

libseq is a library to manipulate biological sequences (nucleotide or amino acid sequences).

###Goal

This library originated as part of my PhD project, which involves molecular evolution of HIV, with respect to resistance development against antiretroviral products.

It provides a limited set of what can be generally found in the biojava, bioperl and biopython projects, in C++. The scope is however limited to efficient and simple sequence manipulation. Application features, such as sequence alignment, are beyond the scope of this library.

###Features

* Nucleotide and AminoAcid sequence class based on std::vector
* Translation from Nucleotide to Aminoacid sequences
* FASTA I/O
* Full support for IUPAC ambiguity symbols
* Simple C++ interface

###Status

The current public release (Jul, 17 2006), [libseq-0.3.tar.gz], contains the basic functionality, which was needed so far for my PhD. It has been exhaustively tested in various projects small and large. You are welcome to add features you find missing.

###Documentation

####Installation

I have switched to cmake for the build process, so I am looking forward to hear from anyone being succesfull to build on anyhthing other than Linux, and I am especially curious for experiences with Microsoft development tools ?
There is an install file with build instructions on Linux, which is follows a standard CMake build and install process.

####Reference
* [API documentation]
* Two simple examples [NTFastaRead.C] and [AAFastaRead.C].

###GitHub
The project is hosted at [GitHub].

[libseq-0.3.tar.gz]: <http://prdownloads.sourceforge.net/libseq/libseq-0.3.tar.gz?download>
[API documentation]: <http://libseq.sourceforge.net/0.3/doxygen/namespaceseq.html>
[NTFastaRead.C]: <http://libseq.sourceforge.net/0.3/example/NTFastaRead.C>
[AAFastaRead.C]: <http://libseq.sourceforge.net/0.3/example/AAFastaRead.C>
[GitHub]: <https://github.com/rega-cev/libseq/>
