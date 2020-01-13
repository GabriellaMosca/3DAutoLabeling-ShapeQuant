# Name of model/process
NAME = FemMembranes
MODELPATH = "Model/CCF/01 FEM Membranes"
SOLIB = usrLib$(NAME).so

# List all object files here
OBJECTS = $(NAME).o qrc_$(NAME).o 


# default target
target: $(SOLIB)

# Call MorphoDynamX to get makefile
include $(shell mdx --resource)/MDXProcess.mk

# Add extra compile flags here
CXXFLAGS += -Wno-unused-local-typedefs -Wno-unused-parameter -Wno-unused-value

# Add extra libraries here
LIBS+= -lmdxFemLib

# Add extra link flags here
LD_FLAGS+=

# Model dependencies
$(NAME).o: $(NAME).cpp $(NAME).hpp Makefile 

$(SOLIB): $(OBJECTS)

#$(NAME).hpp $(NAME).moc 

# Run the model
run: $(SOLIB)
	mdx --model $(MODELPATH) --addlibrary $(SOLIB) $(NAME).mdxv
debug: $(SOLIB)
	mdx --debug --model $(MODELPATH) --addlibrary $(SOLIB) $(NAME).mdxv

