
SA_PATH = ..

CRYPTO_INC_PATH = -I /usr/local/include
INC_PATH = $(CRYPTO_INC_PATH) -I $(SA_PATH)/Ivory-Runtime -I $(SA_PATH)/libOTe -I $(SA_PATH)/libOTe/cryptoTools -I $(SA_PATH)/libOTe/cryptoTools/thirdparty/linux/

IVORY_LIB = -L $(SA_PATH)/Ivory-Runtime/lib -l ivory
LIBOTE_LIB = -L $(SA_PATH)/libOTe/lib -l libOTe -l cryptoTools
BOOST_LIB = -L $(SA_PATH)/libOTe/cryptoTools/thirdparty/linux/boost/stage/lib -l boost_system -l boost_thread
MIRACL_LIB = -L $(SA_PATH)/libOTe/cryptoTools/thirdparty/linux/miracl/miracl/source -l miracl
CRYPTO_LIB = -L /usr/local/lib -lcryptopp
OTHER_LIB = -l pthread

LIB_PATH = $(IVORY_LIB) $(LIBOTE_LIB) $(BOOST_LIB) $(MIRACL_LIB) $(CRYPTO_LIB) $(OTHER_LIB)

CC = g++
CFLAGS = -Wall -std=c++14

SOURCES = main.cpp spectrum_manager.cpp generator.cpp args.cpp location.cpp split.cpp timer.cpp utils.cpp tables.cpp path_loss_table.cpp su_server.cpp key_server.cpp debug_print.cpp shared.cpp
OBJECTS = $(SOURCES:.cpp=.o)
NAME = s2pc

all: $(NAME)

$(NAME): $(OBJECTS)
	$(CC) $(CFLAGS) $(INC_PATH) $(OBJECTS) $(LIB_PATH) -o $(NAME)

.cpp.o:
	$(CC) $(CFLAGS) $(INC_PATH) -c $< -o $@

clean:
	rm -f $(NAME) *.o
