# compiler and its flags
CC = mpiCC
CFLAGS = -Wall -O3 -g -w
#CFLAGS = -Wall -O3 -g -std=c++0x -w
#CFLAGS = -Wall -O3 -g -pg

# project name
PROJECT = MMM_v14_6

# directories
OBJDIR = object
SRCDIR = source

# libraries
LDFLAGS = -lm -lrt

# files and folders
SRCS = $(shell find $(SRCDIR) -name '*.cc')
OBJS = $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.o,$(SRCS))

# targets

#$(PROJECT) : clean $(OBJS) $(OBJDIR)/kdtree.o
#	$(CC) $(CFLAGS) $(OBJS) $(OBJDIR)/kdtree.o $(LDFLAGS) -o $@

$(PROJECT) : $(OBJS) $(OBJDIR)/kdtree.o
	$(CC) $(CFLAGS) $(OBJS) $(OBJDIR)/kdtree.o $(LDFLAGS) -o $@

$(OBJDIR)/%.o : $(SRCDIR)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/kdtree.o : $(SRCDIR)/kdtree.c $(SRCDIR)/kdtree.h
	gcc -std=c89 -pedantic -Wall -g -I.. -c $< -o $@

#clean :
#	rm -f $(OBJDIR)/*.o

# dependencies
DEPS = $(OBJS:.o=.d)
-include $(DEPS)
$(DEPS): $(OBJDIR)/%.d : $(SRCDIR)/%.cc
	@echo Refreshing dependency list for $<...
	@gcc -MM -MF $@ -MT "$(@:.d=.o)" $<


