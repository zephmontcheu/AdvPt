add_executable(create-htgfile)

target_sources(create-htgfile
PRIVATE
    create-htgfile.cpp)

target_link_libraries(create-htgfile
PRIVATE
    oktal)

add_executable(poisson)

target_sources(poisson
PRIVATE
    poisson.cpp)

target_link_libraries(poisson
PRIVATE
    oktal)