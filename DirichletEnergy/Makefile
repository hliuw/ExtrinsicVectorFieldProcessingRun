CXX = cl
CXXFLAGS = /nologo /EHsc /W4 /std:c++17
INCLUDES = /Iinclude
TARGET = DirichletEnergy.exe
OBJS = main.obj MeshIO.obj VectorFieldInterpolation.obj VectorDirichletEnergy.obj

all: $(TARGET)

$(TARGET): $(OBJS)
    $(CXX) $(CXXFLAGS) /Fe$(TARGET) $(OBJS)

main.obj: main.cpp include\Geometry.h include\MeshIO.h include\VectorDirichletEnergy.h include\VectorFieldInterpolation.h
    $(CXX) $(CXXFLAGS) $(INCLUDES) /c main.cpp /Fo$@

MeshIO.obj: src\MeshIO.cpp include\MeshIO.h include\Geometry.h include\happly.h include\VectorFieldInterpolation.h
    $(CXX) $(CXXFLAGS) $(INCLUDES) /c src\MeshIO.cpp /Fo$@

VectorFieldInterpolation.obj: src\VectorFieldInterpolation.cpp include\VectorFieldInterpolation.h include\Geometry.h include\MathUtils.h
    $(CXX) $(CXXFLAGS) $(INCLUDES) /c src\VectorFieldInterpolation.cpp /Fo$@

VectorDirichletEnergy.obj: src\VectorDirichletEnergy.cpp include\VectorDirichletEnergy.h include\Geometry.h include\MathUtils.h
    $(CXX) $(CXXFLAGS) $(INCLUDES) /c src\VectorDirichletEnergy.cpp /Fo$@

clean:
    -del /Q $(OBJS) $(TARGET) hello.exe
