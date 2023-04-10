cmake_minimum_required(VERSION 2.8.12)

project (vorotest )
#SET(CMAKE_C_COMPILER /usr/bin/clang )
#SET(CMAKE_CXX_COMPILER /usr/bin/clang++ )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -O2 -g ")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++1y -Wall ")

include_directories("${PROJECT_SOURCE_DIR}/include")
file(GLOB 
	CPP_SOURCE "src/*.cpp"	
	C_SOURCE   "src/*.c"	
	CC_SOURCE  "src/*.cc"
	)
add_executable(${PROJECT_NAME} 
	${CPP_SOURCE} 
	${C_SOURCE} 
	${CC_SOURCE} 
	)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(vorotest PRIVATE Threads::Threads)

