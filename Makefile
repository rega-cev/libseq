OPTFLAGS = -O3 -march=pentium4
ARCH=intel
#OPTFLAGS = -O3 -march=opteron
#ARCH=opteron

CFLAGS =  -Wall \
          -I$(SRCDIR)/evolution -I$(SRCDIR)/sequence $(OPTFLAGS)
LIBS = 

LIBRARY=$(BUILDDIR)/libseq.a

FILES=sequence/Nucleotide sequence/AminoAcid sequence/NTSequence \
      sequence/AASequence sequence/Codon sequence/Mutation \
      sequence/CodingSequence \
      evolution/NucleotideSubstitutionModel 

SRCDIR=src
BUILDDIR=obj-$(ARCH)
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

-include $(DEPS)
