# On purge la liste des suffixes utilisés pour les rôles implicites
.SUFFIXES:

# On ajoute simplements les extensions dont on a besoin
.SUFFIXES:.cpp .o

# Nom de l’exécutable
EXEC=tp4

# Liste des fichiers sources séparés par des espaces
SOURCES=main.cpp

# Liste des fichiers objets
OBJETS=$(SOURCES:%.cpp=%.o)

# Compilateur et options de compilation
CCPP=g++
CFLAGS=-Wall -pedantic -ffast-math -I /usr/X11R6/include


LFLAGS= -L . -L /usr/X11R6/lib  -lpthread -lX11 -lXext -Dcimg_use_xshm  -lm

# Rôle explicite de construction de l’exéutable
$(EXEC):$(OBJETS) Makefile
	$(CCPP) -o  $(EXEC) $(OBJETS) $(LFLAGS)
.cpp.o:
	$(CCPP) $(CFLAGS) -c $< -o $@

clean:
	rm $(OBJETS)
clear:
	rm $(EXEC)
depend:
	sed -e "/^#DEPENDANCIES/,$$ d" Makefile >dependances
	echo "#DEPENDANCIES" >> dependances
	$(CCPP) -MM $(SOURCES) >> dependances
	cat dependances >Makefile
	rm dependances

# Dépendances
main.o: main.cpp 




