cmake_minimum_required(VERSION 3.20.5)
project(rend LANGUAGES CXX)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(rend src/main.cpp src/canvas.cpp src/ray.cpp src/sphere.cpp )
target_include_directories(rend PRIVATE include/)
target_compile_options(rend PRIVATE -Wall -Werror -Wpedantic -fsanitize=address -fsanitize=undefined -g)
target_compile_features(rend PRIVATE cxx_std_20)
target_link_options(rend PRIVATE -fsanitize=address -fsanitize=undefined)
target_link_libraries(rend PRIVATE CONAN_PKG::fmt)

#testing
# include(CTest)
# include(Catch)
# catch_discover_tests(tests)

# copy compile commands to project root dir
add_custom_target(
    copy-compile-commands ALL
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_LIST_DIR}
    )
