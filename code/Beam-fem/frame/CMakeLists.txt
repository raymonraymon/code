PROJECT (FRAME)
SET(SRC_LIST beam.cpp PlaneDisplacement.cpp PlaneElement.cpp PlaneForce.cpp PlanePoint.cpp PlaneStructer.cpp)

#INCLUDE_DIRECTORIES("C:/Program Files/boost_1_62_0/include/boost-1_62")
#LINK_DIRECTORIES("C:/Program Files/boost_1_62_0/lib")
SET(Boost_ADDITIONAL_VERSIONS "1.62" "1.62.0")
SET(BOOST_ROOT "C:/Program Files/boost_1_62_0")

# option
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

find_package(Boost 1.62.0 COMPONENTS REQUIRED)   # 最好是加上REQUIRED，否则找不到某些库也不会报错，链接会出问题

if(Boost_FOUND)
  INCLUDE_DIRECTORIES(${BOOST_ROOT}/include/boost-1_62)
  target_link_libraries (${BOOST_ROOT}/lib)
endif (Boost_FOUND)

ADD_EXECUTABLE(frame ${SRC_LIST})