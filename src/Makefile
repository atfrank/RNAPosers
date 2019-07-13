#Sean M. Law
#Aaron T. Frank
#Jingru Xie

include ../makefile.defs

OBJS = libmoletools.$(LIBEXT)

all:	$(BINDIR) libmoletools.$(LIBEXT) featurize 

clean:
	rm -rf *.o
	rm -rf *stackdump
	rm -rf ../bin/*
	rm -rf *~

$(BINDIR):
	mkdir -p $(BINDIR)

libmoletools.$(LIBEXT):
	$(MAKE) -C ../lib libmoletools.$(LIBEXT)

featurize: featurize.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LIBFLAGS) $(INCFLAGS) -o $(BINDIR)/$@ $< $(LINKFLAGS)

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(INCFLAGS) -c $<
