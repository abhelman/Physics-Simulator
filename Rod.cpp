///////////////////////////////////////////////////////////////////////////////
// Sphere.cpp
// ==========
// Sphere for OpenGL with (radius, sectors, stacks)
// The min number of sectors is 3 and the min number of stacks are 2.
//
//  AUTHOR: Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2017-11-01
// UPDATED: 2020-05-20
///////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
#include <windows.h>    // include windows.h to avoid thousands of compile errors even though this class is not depending on Windows
#endif

#include <iostream>
#include <iomanip>
#include <cmath>
#include "Rod.h"
#include <glad/glad.h> 
#include <GLFW/glfw3.h>


///////////////////////////////////////////////////////////////////////////////
// ctor
///////////////////////////////////////////////////////////////////////////////
Rod::Rod(float radius, int sectors, int stacks, float step)
{
}



///////////////////////////////////////////////////////////////////////////////
// setters
///////////////////////////////////////////////////////////////////////////////
void Rod::set(float radius, int sectors, int stacks, float step)
{
    this->radius = radius;
    this->sectorCount = sectors;
    this->stackCount = stacks;
    this->widthStep = step;

    float PI = 3.14;
    float x, y, z, xy;                              // vertex position

    float sectorStep = 2 * PI / sectorCount;
    float stackStep = PI / stackCount;
    float sectorAngle, stackAngle;

    for (int i = 0; i <= stackCount * 2; ++i)
    {
        stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
        xy = radius * cosf(stackAngle);             // r * cos(u)
        z = radius * sinf(stackAngle);              // r * sin(u)

        // add (sectorCount+1) vertices per stack
        // the first and last vertices have same position and normal, but different tex coords
        for (int j = 0; j <= sectorCount; j++)
        {
            sectorAngle = 0;
            x = xy * cosf(sectorAngle);// +j * 0.1f;;             // r * cos(u) * cos(v)
            y = xy * sinf(sectorAngle) + j * widthStep;;             // r * cos(u) * sin(v)
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);

        }
    }
    int k1, k2;
    for (int i = 0; i < stackCount * 2; ++i)
    {
        k1 = i * (sectorCount + 1);     // beginning of current stack
        k2 = k1 + sectorCount + 1;      // beginning of next stack

        for (int j = 0; j < sectorCount; ++j, ++k1, ++k2)
        {
            // 2 triangles per sector excluding first and last stacks
            // k1 => k2 => k1+1

            indices.push_back(k1);
            indices.push_back(k2);
            indices.push_back(k1 + 1);


            indices.push_back(k1 + 1);
            indices.push_back(k2);
            indices.push_back(k2 + 1);


        }
    }

    unsigned int VBO;
    glGenBuffers(1, &VBO);

    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, getVertexSize(), getVertices(), GL_STATIC_DRAW);

    unsigned int EBO;
    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, getIndexSize(), getIndices(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


///////////////////////////////////////////////////////////////////////////////
// draw a sphere in VertexArray mode
// OpenGL RC must be set before calling it
///////////////////////////////////////////////////////////////////////////////
void Rod::draw() const
{
    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, getIndexCount(), GL_UNSIGNED_INT, 0);
}

Rod::~Rod() {

    glDeleteVertexArrays(1, &VAO);

}

