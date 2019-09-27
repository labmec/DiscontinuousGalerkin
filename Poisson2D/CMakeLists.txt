
set(Poisson2D_SRCS main.cpp)

add_executable(Poisson2D ${Poisson2D_SRCS} )

target_link_libraries(Poisson2D pz)