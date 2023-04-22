cmake_minimum_required(VERSION 3.21)

FILE(GLOB SOURCE_FILE *.cpp *.h)

set(COLMAP_DIR "D:\\git\\colmap\\build\\install\\release\\share\\colmap")
set(exiv2_DIR "D:\\git\\exiv2-0.27.5-2019msvc64\\lib\\cmake\\exiv2")
#set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_BUILD_TYPE Release) # Debug, Release, RelWithDebInfo, MinSize

find_package(Boost REQUIRED COMPONENTS program_options filesystem)
find_package(COLMAP REQUIRED)
find_package(PROJ CONFIG REQUIRED)
find_package(exiv2 CONFIG REQUIRED)
find_package(OpenCV CONFIG REQUIRED)

message(STATUS "1: ${COLMAP_LIBRARIES}")
message(STATUS "2: ${PROJ_LIBRARIES}")
message(STATUS "3: boost:  ${Boost_INCLUDE_DIRS} ${Boost_LIBRARY_DIR_RELEASE}")
message(STATUS "4: opencv:  ${OpenCV_INCLUDE_DIRS} ${OpenCV_LIBS}")

ADD_EXECUTABLE(${PROJECT_NAME} ${SOURCE_FILE})
target_include_directories(${PROJECT_NAME} PUBLIC ${COLMAP_INCLUDE_DIRS} PROJ::proj ${OpenCV_INCLUDE_DIRS})
target_link_directories(${PROJECT_NAME} PUBLIC ${COLMAP_LINK_DIRS} PROJ::proj)
target_link_libraries(${PROJECT_NAME} general ${Boost_LIBRARY_DIR_RELEASE} ${COLMAP_LIBRARIES} ${PROJ_LIBRARIES} PROJ::proj Boost::program_options Boost::filesystem exiv2lib ${OpenCV_LIBS})