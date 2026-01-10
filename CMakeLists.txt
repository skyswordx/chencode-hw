# CMake 最低版本要求
cmake_minimum_required(VERSION 3.10)

# 项目名称
project(chencode_sim C)

# 设置 C 语言标准
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# 查找 src 目录下的所有源文件 (.c)
file(GLOB SOURCES "src/*.c")

# 生成可执行文件
add_executable(chencode_sim ${SOURCES})

# 添加 include 目录，以便源文件可以找到 .h 头文件
target_include_directories(chencode_sim PRIVATE ${CMAKE_SOURCE_DIR}/include)

# 链接数学库
# Windows MSVC 不需要显式链接 libm，但 GCC/Clang 需要
if(NOT MSVC)
    target_link_libraries(chencode_sim m)
endif()

# 创建 output 目录用于存放 CSV 文件
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output)

# 提示：如何在命令行中构建此项目：
# Windows (MSVC):
#   1. mkdir build && cd build
#   2. cmake ..
#   3. cmake --build . --config Release
#
# Linux/macOS (GCC):
#   1. mkdir build && cd build
#   2. cmake ..
#   3. make -j$(nproc)