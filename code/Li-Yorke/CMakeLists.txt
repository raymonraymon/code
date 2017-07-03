PROJECT (HOMOTOPY)
SET(SRC_LIST Homotopy.cpp LiYorke.cpp TargetFun.cpp)

INCLUDE_DIRECTORIES("C:/Program Files/boost_1_62_0/include/boost-1_62")
LINK_DIRECTORIES("C:/Program Files/boost_1_62_0/lib")

MESSAGE(STATUS "This is BINARY dir " ${HELLO_BINARY_DIR})
MESSAGE(STATUS "This is SOURCE dir "${HELLO_SOURCE_DIR})

ADD_EXECUTABLE(Homotopy ${SRC_LIST})