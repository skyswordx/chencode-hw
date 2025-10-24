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

# 链接数学库 (libm)，因为代码中使用了 pow, sqrt, log, exp, cos 等函数
target_link_libraries(chencode_sim m)

# 提示：如何在命令行中构建此项目：
# 1. mkdir build
# 2. cd build
# 3. cmake ..
# 4. make -j$(nproc)