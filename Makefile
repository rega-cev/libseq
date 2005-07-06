CFLAGS =  -Wall \
          -I$(SRCDIR)/evolution -I$(SRCDIR)/sequence \
	  -g
	  -O3 -march=pentium4 -funroll-loops
#	  -O3 -march=opteron
#	  -O3 -march=pentium3
LIBS = 

LIBRARY=libseq.a

FILES=sequence/Nucleotide sequence/AminoAcid sequence/NTSequence \
      sequence/AASequence sequence/Codon sequence/CodingSequence \
      evolution/NucleotideSubstitutionModel 

SRCDIR=src
BUILDDIR=obj
DEPDIR=dep

SOURCES=$(patsubst %,$(SRCDIR)/%.C,$(FILES))
OBJECTS=$(patsubst %,$(BUILDDIR)/%.o,$(FILES))
DEPS=$(patsubst %,$(DEPDIR)/%.d,$(FILES))

all: $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	$(AR) rcu $(LIBRARY) $(OBJECTS)

$(BUILDDIR)/%.o : $(SRCDIR)/%.C
	@mkdir -p $(BUILDDIR)/$(dir $*)
	@mkdir -p $(DEPDIR)/$(dir $*)
	$(CXX) $(CFLAGS) -MP -MD -MF $(DEPDIR)/$(dir $*)/$(notdir $*).d \
		-o $@ -c $<

documentation:
	doxygen

clean:
	rm -r $(DEPDIR) $(BUILDDIR) $(LIBRARY) `find . -name "*~"`

distclean: clean
	rm -r doc/*

config/config.me:
	@echo Please run ./configure first to create the config/config.me file.
	@echo ./configure --help gives an overview of the possible command line parameters.
	@echo make will now exit with an error.
	exit 1

-include $(DEPS)
