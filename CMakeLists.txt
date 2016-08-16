# This is demo on how to link libpointmatcher to an
# external project.

cmake_minimum_required (VERSION 2.8)
project(slam)

#Enable support for C++11  
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")  

## 对寻找G2O是必要的
#设定二进制文件路径  
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)  
#设定库文件编译路径  
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)  
#设定.cmake文件存放路径  
#SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "{CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")  
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)  
#ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src) 

find_package(libpointmatcher REQUIRED)
find_package(G2O REQUIRED)

find_package(Eigen3 REQUIRED)  
find_package(CSparse REQUIRED)

include_directories(${CSPARSE_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${libpointmatcher_INCLUDE_DIRS})
include_directories(include)

# 查找当前目录下的所有源文件
# 并将名称保存到 DIR_SRCS 变量
aux_source_directory(. DIR_SRCS)

#add_subdirectory(octomap)
#add_subdirectory(octomap/math)
add_subdirectory(dbscan)
add_subdirectory(ransac)
#add_subdirectory(csparse)
SET(G2O_LIBS g2o_cli g2o_ext_freeglut_minimal g2o_simulator g2o_solver_slam2d_linear g2o_types_icp g2o_types_slam2d g2o_core g2o_interface g2o_solver_csparse g2o_solver_structure_only g2o_types_sba g2o_types_slam3d g2o_csparse_extension g2o_opengl_helper g2o_solver_dense g2o_stuff g2o_types_sclam2d g2o_parser g2o_solver_pcg g2o_types_data g2o_types_sim3 cxsparse )

add_executable(${PROJECT_NAME} ${DIR_SRCS})

target_link_libraries(${PROJECT_NAME} ${libpointmatcher_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${G2O_LIBS})
#target_link_libraries(${PROJECT_NAME} octomap)
target_link_libraries(${PROJECT_NAME} dbscan)
target_link_libraries(${PROJECT_NAME} ransac)
