CFLAGS =  -Wall \
          -I$(SRCDIR)/evolution -I$(SRCDIR)/sequence \
	  -O3 -march=pentium3
#	  -O3 -march=pentium4 -funroll-loops
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
	rm -rf $(DEPDIR) $(BUILDDIR) $(BINARY) `find . -name "*~"`

distclean: clean
	rm -rf config/config.me config/config.log


config/config.me:
	@echo Please run ./configure first to create the config/config.me file.
	@echo ./configure --help gives an overview of the possible command line parameters.
	@echo make will now exit with an error.
	exit 1

-include $(DEPS)
