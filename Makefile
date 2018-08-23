
SA_PATH = ..

INC_PATH = -I $(SA_PATH)/Ivory-Runtime -I $(SA_PATH)/libOTe -I $(SA_PATH)/libOTe/cryptoTools -I $(SA_PATH)/libOTe/cryptoTools/thirdparty/linux/ -I $(SA_PATH)/itm/c/

IVORY_LIB = -L $(SA_PATH)/Ivory-Runtime/lib -l ivory
LIBOTE_LIB = -L $(SA_PATH)/libOTe/lib -l libOTe -l cryptoTools
BOOST_LIB = -L $(SA_PATH)/libOTe/cryptoTools/thirdparty/linux/boost/stage/lib -l boost_system -l boost_thread
MIRACL_LIB = -L $(SA_PATH)/libOTe/cryptoTools/thirdparty/linux/miracl/miracl/source -l miracl
ITM_LIB = -L $(SA_PATH)/itm/c/ -l itm
OTHER_LIB = -l pthread

LIB_PATH = $(IVORY_LIB) $(LIBOTE_LIB) $(BOOST_LIB) $(MIRACL_LIB) $(ITM_LIB) $(OTHER_LIB)

CC = g++
CFLAGS = -Wall -std=c++14 $(INC_PATH)

SOURCES = main.cpp spectrum_manager.cpp generator.cpp args.cpp location.cpp split.cpp timer.cpp utils.cpp tables.cpp path_loss_table.cpp
OBJECTS = $(SOURCES:.cpp=.o)
NAME = s2pc

all: $(NAME)

$(NAME): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIB_PATH) -o $(NAME)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(NAME) *.o
