
set(Poisson2D_SRCS main.cpp)

add_executable(Poisson2D ${Poisson2D_SRCS} )

target_link_libraries(Poisson2D pz)
target_include_directories(Poisson2D PUBLIC ${PZ_INCLUDE_DIRS})