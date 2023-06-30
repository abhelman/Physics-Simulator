///////////////////////////////////////////////////////////////////////////////
// Sphere.h
// ========
// Sphere for OpenGL with (radius, sectors, stacks)
// The min number of sectors is 3 and The min number of stacks are 2.
//
//  AUTHOR: Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2017-11-01
// UPDATED: 2020-05-20
///////////////////////////////////////////////////////////////////////////////

#ifndef GEOMETRY_SPHERE_H
#define GEOMETRY_SPHERE_H

#include <vector>

class Sphere
{
public:
    // ctor/dtor
    Sphere(float radius=1.0f, int sectorCount=36, int stackCount=18);
    ~Sphere();

    void set(float radius, int sectorCount, int stackCount);
    // for vertex data
    unsigned int getIndexCount() const      { return (unsigned int)indices.size(); }
    unsigned int getVertexSize() const      { return (unsigned int)vertices.size() * sizeof(float); }
    const float* getVertices() const        { return vertices.data(); }
    unsigned int getIndexSize() const { return (unsigned int)indices.size() * sizeof(int); }
    const unsigned int* getIndices() const  { return indices.data(); }
    
    // draw in VertexArray mode
    void draw() const;                                  // draw surface
    
protected:

private:
    
    // memeber vars
    float radius;
    int sectorCount;                        // longitude, # of slices
    int stackCount;                         // latitude, # of stacks
    std::vector<float> vertices;
    std::vector<unsigned int> indices;
    unsigned int VAO;
};

#endif
