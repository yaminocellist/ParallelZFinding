cmake_minimum_required(VERSION 3.5)
project(Para_dPhi)

find_package(nlohmann_json 3.11.3 REQUIRED)
find_package(ROOT REQUIRED)

add_executable(Para_dPhi para_dPhi.cpp)
# Link against ROOT libraries
target_link_libraries(Para_dPhi PUBLIC ${ROOT_LIBRARIES})
# Add ROOT include directories
target_include_directories(Para_dPhi PUBLIC ${ROOT_INCLUDE_DIRS})
# Add ROOT flags
target_compile_options(Para_dPhi PUBLIC
    -std=c++17
    -m64
    # -pipe
    # -fsigned-char
    # -fno-common
    # -Qunused-arguments
    -pthread
    -stdlib=libc++
)



# root-config --libs
# root-config --cflags