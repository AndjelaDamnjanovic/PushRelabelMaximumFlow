cmake_minimum_required(VERSION 3.16.0)

set(HDR_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Headers")
set(SRC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Source")
set(UI_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Forms")

project(graphvisual VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC OFF)

find_package(Qt5 COMPONENTS Core)
find_package(Qt5 COMPONENTS Widgets)
find_package(Qt5 COMPONENTS Gui) 

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(RESOURCES
    Resources.qrc
)

set(SOURCES
    ${SRC_PATH}/edge.cpp
    ${SRC_PATH}/graph.cpp
    ${SRC_PATH}/main.cpp
    ${SRC_PATH}/graphwindow.cpp
    ${SRC_PATH}/mytextedit.cpp
    ${SRC_PATH}/node.cpp
    ${SRC_PATH}/algorithm.cpp
    ${SRC_PATH}/graphicnode.cpp
    ${SRC_PATH}/graphtable.cpp
    ${SRC_PATH}/graphicedge.cpp
    ${SRC_PATH}/popup.cpp
)

set(HEADERS
    ${HDR_PATH}/edge.h
    ${HDR_PATH}/graph.h
    ${HDR_PATH}/graphwindow.h
    ${HDR_PATH}/mytextedit.h
    ${HDR_PATH}/node.h
    ${HDR_PATH}/algorithm.h
    ${HDR_PATH}/graphicnode.h
    ${HDR_PATH}/graphtable.h
    ${HDR_PATH}/graphicedge.h
    ${HDR_PATH}/popup.h
)

set(UIS
    ${UI_PATH}/graphwindow.ui
    ${UI_PATH}/popup.ui
)

qt5_wrap_ui(UI_HEADERS ${UIS})

add_executable(${PROJECT_NAME} ${HEADERS} ${SOURCES} ${UI_HEADERS} ${RES_SOURCES})			

target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::Core Qt5::Widgets Qt5::Gui)

target_include_directories(${PROJECT_NAME} PUBLIC "${HDR_PATH}")

