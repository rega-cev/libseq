SEQLIBPATH = ..
SEQINCLUDE = -I$(SEQLIBPATH)/src/sequence -I$(SEQLIBPATH)/src/evolution
INCLUDE = $(SEQINCLUDE) \
          -I$(SRCDIR)

CFLAGS =  -O3 -Wall \
	  $(INCLUDE) \
	  -march=pentium3

LIBS = -L$(SEQLIBPATH) -lseq

BINARIES=NTFastaRead AAFastaRead

FILES=

SRCDIR=src
BUILDDIR=obj
DEPDIR=dep

SOURCES=$(patsubst %,$(SRCDIR)/%.C,$(FILES))
OBJECTS=$(patsubst %,$(BUILDDIR)/%.o,$(FILES))
DEPS=$(patsubst %,$(DEPDIR)/%.d,$(FILES))

BINSOURCES=$(patsubst %,$(SRCDIR)/%.C,$(BINARIES))
BINOBJECTS=$(patsubst %,$(BUILDDIR)/%.o,$(BINARIES))
BINDEPS=$(patsubst %,$(DEPDIR)/%.d,$(BINARIES))
BUILD_BINARIES=$(patsubst %,$(BUILDDIR)/%,$(BINARIES))

all: $(BUILD_BINARIES)

$(BUILD_BINARIES): $(BINOBJECTS) $(OBJECTS) $(SEQLIBPATH)/libseq.a
	@for i in $(BINARIES); do \
		echo $(CXX) -o $(BUILDDIR)/$$i $(BUILDDIR)/$$i.o $(OBJECTS) $(LIBS); \
		$(CXX) -o $(BUILDDIR)/$$i $(BUILDDIR)/$$i.o $(OBJECTS) $(LIBS); \
	done

$(BUILDDIR)/%.o : $(SRCDIR)/%.C
	@mkdir -p $(BUILDDIR)/$(dir $*)
	@mkdir -p $(DEPDIR)/$(dir $*)
	$(CXX) $(CFLAGS) -MP -MD -MF $(DEPDIR)/$(dir $*)/$(notdir $*).d \
		-o $@ -c $<

clean:
	rm -rf $(DEPDIR) $(BUILDDIR) $(BINARY) `find . -name "*~"`

distclean: clean
	rm -rf config/config.me config/config.log


config/config.me:
	@echo Please run ./configure first to create the config/config.me file.
	@echo ./configure --help gives an overview of the possible command line parameters.
	@echo make will now exit with an error.
	exit 1

-include $(BINDEPS)
