#include <iostream>
#include <stdlib.h>
#include <functional>

#include <string>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "Shader.h"
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Sphere.h"
#include "Rod.h"
#include <fstream>

double BASE = -0.3;
double DT = 0.001;
double Kc = 100000;
// 0 = tissue, 1=muscle, 2=muscle, 3 = bone
double bVals[] = { 0.0,0.25,0.25,0.0,0.0 };
double cVals[] = { 0.0,0.0,3.14,0.0,0.0 };
double kVals[] = { 1000.0, 5000.0,5000.0,20000.0,0.0 };
double edgeLength = 0.1;
int COM = 11;
double height = 4.0;
int EMPTY = 4;


class Mass {

public:
    double m;
    // Convert to vects!
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double ax;
    double ay;
    double az;
    double Fx;
    double Fy;
    double Fz;
    double dt;
    double KE;
    double x0;
    double y0;
    double z0;
    int group;
    double PE;
    

    Mass(double mass, double pos_x, double pos_y, double pos_z, int g, int init_force = 0) {
        group = g;
        m = mass;
        x = pos_x;
        y = pos_y;
        z = pos_z;
        x0 = pos_x;
        y0 = pos_y;
        z0 = pos_z;
        vx = 0.0;
        vy = 0.0;
        vz = 0.0;
        ax = 0.0;
        ay = 0.0;
        az = 0.0;
        Fx = 0.0;
        Fy = 0.0;
        Fz = 0.0;
        if (init_force) {
            Fx = 10.0;
            Fy = 0.0;
            Fz = 10.0;
        }

        dt = DT;
        KE = 1 / 2 * m * pow(sqrt(pow(vx, 2)+ pow(vy, 2)+ pow(vz, 2)), 2);
        PE = 0;
    }

    void updateGroup(int g) {
        group = g;
        x = x0;
        y = y0;
        z = z0;
        vx = 0.0;
        vy = 0.0;
        vz = 0.0;
    }


    void applyForce(double forcex, double forcey, double forcez) {
        Fx += forcex;
        Fy += forcey;
        Fz += forcez;
    }

    void finish() {

        Fy += (m * -1.0 * 9.81);
        double norm = sqrt(vx * vx + vz * vz);

        if (y < BASE) {
            if (norm > 0) { // if moving, kinetic friction
                Fx -= (vx / norm) * 1.0 * 0.57 * abs(Fy);
                Fz -= (vz / norm) * 1.0 * 0.57 * abs(Fy);

            }
            else { // not moving, apply static friction
                if (abs(Fy) * 0.74 > sqrt(Fx * Fx + Fz * Fz)) {
                    Fx = 0.0;
                    Fz = 0.0;
                }
            }
        }

        if (y < BASE) {
            Fy += -1.0 * Kc * (y - BASE);
        }



        ax = (double)Fx / (double)m;
        vx = 0.999*vx + (ax * dt);
        x = x + (vx * dt);

        ay = (double)Fy / (double)m;
        vy = 0.999*vy + (ay * dt);
        y = y + (vy * dt);

        az = (double)Fz / (double)m;
        vz = 0.999 * vz + (az * dt);
        z = z + (vz * dt);

        Fx = 0.0;
        Fy = 0.0;
        Fz = 0.0;
    }

    double updateKE() {
        KE = (1.0 / 2.0) * m * pow(sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2)), 2);
        return KE;
    }

    double updatePE() {
        PE = m * 9.81 * (1.0 + y);
        if (y < BASE) {
            PE += (0.5 * Kc * (y - BASE) * (y - BASE));
        }
        return PE;
    }

};
class Spring {
public:
    double k;
    double l;
    double dl;
    double PE;
    Mass* m1;
    Mass* m2;
    double x;
    double y;
    double z;
    double n;
    double a;
    double b;
    double w;
    double c;
    double l0;
    double t;

    Spring(Mass* pm1, Mass* pm2) {
        dl = 0.0;
        m1 = pm1;
        m2 = pm2;
        x = m1->x - m2->x;
        y = m1->y - m2->y;
        z = m1->z - m2->z;
        n = sqrt(x * x + y * y + z * z);
        x = x / n;
        y = y / n;
        z = z / n;
        
        int category = m1->group;

        k = kVals[category];
        b = bVals[category];
        c = cVals[category];
        
        l = n;
        a = l;
        t = 0;
        w = 3.14*2.0;
        l0 = l*(1.0 + b*sin(t*w + c));
        PE = 0;     
    }
    void reset() {

        int category = m1->group;
        k = kVals[category];
        b = bVals[category];
        c = cVals[category];
        t = 0;
        l0 = l * (1.0 + b * sin(t * w + c));

    }
    void updateL() {
        x = m1->x - m2->x;
        y = m1->y - m2->y;
        z = m1->z - m2->z;

        n = sqrt(x * x + y * y + z * z);
        dl = n - l0;
        t += DT;
        l0 = l * (1.0 + b * sin(t * w + c));

        x = x / n;
        y = y / n;
        z = z / n;

        PE = 0.5 * k * dl * dl;
    }

    double updatePE() {
        return PE;
    }

    void calcForce() {
        double forcex, forcey, forcez;
       
        forcex = k * dl * x;
        forcey = k * dl * y;
        forcez = k* dl* z;

        m1->applyForce(-1.0 * forcex, -1.0 * forcey, -1.0 * forcez);
        m2->applyForce(forcex, forcey, forcez);

        }
};

class Platform {

    public:
        unsigned int VAOH, VAOV;
        Platform() {
            float platformPosH[] = {
                -3.0f, BASE, 0.005f,
                3.0f, BASE, 0.005f,
                3.0f, BASE, -0.005f,
                -3.0f, BASE, 0.005f,
                -3.0f, BASE, -0.005f,
                3.0f, BASE, -0.005f,
            };

            float platformPosV[] = {
                -0.005f, BASE, 3.0f,
                0.005f, BASE, 3.0f,
                0.005f, BASE, -3.0f,
                -0.005f, BASE, 3.0f,
                -0.005f, BASE, -3.0f,
                0.005f, BASE, -3.0f,
            };

            unsigned int VBOH,VBOV;
            glGenBuffers(1, &VBOH);

            // stores configured data
            glGenVertexArrays(1, &VAOH);
            // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
            glBindVertexArray(VAOH);
            // stores configured data

            glBindBuffer(GL_ARRAY_BUFFER, VBOH);
            glBufferData(GL_ARRAY_BUFFER, 6 * 3 * sizeof(float), &platformPosH[0], GL_STATIC_DRAW);

            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(0);

            glBindVertexArray(0);
            
            glGenBuffers(1, &VBOV);
            glGenVertexArrays(1, &VAOV);
            glBindVertexArray(VAOV);
            glBindBuffer(GL_ARRAY_BUFFER, VBOV);
            glBufferData(GL_ARRAY_BUFFER, 6 * 3 * sizeof(float), &platformPosV[0], GL_STATIC_DRAW);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(0);
            glBindVertexArray(0);
        }
        void draw(unsigned int modelLoc, unsigned int colorLocation) {
            //glUniform4f(colorLocation, 0.36389f, 0.618501f, 0.782349f, 1.0f);
            glm::mat4 model = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
            //model = glm::rotate(model, glm::radians(-10.0f), glm::vec3(0.0f, 1.0f, 0.0f));
            //model = glm::scale(model, glm::vec3(10.0, 1.0, 0.1));

            glUniform4f(colorLocation, 0.92f, 0.92f, 0.92f, 0.8f);
            for (int i = -30; i < 30; i++) {
                model = glm::mat4(1.0f);
                model = translate(model, glm::vec3(0.0f,0.0f,edgeLength*i));
                glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
                glBindVertexArray(VAOH);
                glDrawArrays(GL_TRIANGLES, 0, 6);
                model = glm::mat4(1.0f);
                model = translate(model, glm::vec3(edgeLength * i, 0.0f, 0.0f));
                glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
                glBindVertexArray(VAOV);
                glDrawArrays(GL_TRIANGLES, 0, 6);
            }
        }
};

class Robot {
    private:
    public:
        std::vector<Mass> massPositions;
        std::vector<Spring> springs;
        Sphere sphere = Sphere(0.015,36.0f,18.0f);
        // 0.1 basic string length
        Rod rod = Rod();

        Robot(int l, int w, int h, std::vector<double>& encoding) {
            rod.set(0.002f, 36.0f, 18.0f, edgeLength / 36.0f);
            //float x = BASE;
            //float y, z;
            // masses in one row
            //int m = (l + 1);
            // total masses
            //int masses = m*(w+1) * (h + 1);
            // height of stack
            
            float x = 0.05;
            float y, z;
            int masses = 32;
            int m = 8;
            int n = 2;
            for (int i = 0; i < masses; i += m) {
                z = 0.05;
                for (int j = 0; j < m; j += n) {
                    y = BASE + height + 0.2;
                    for (int k = 0; k < n; k++) {
                        massPositions.push_back(Mass(0.1, x, y, z, 0));
                        //printf("%f %f %f, ", x, y, z);
                        y = BASE + height + 0.1;
                    }
                    z -= 0.1;
                }
                x -= 0.1;
            }


            x += 0.1;
            z += 0.1;
            for (int i = 0; i >= 0; i--) {
                y = BASE + height + i * 0.1;
                massPositions.push_back(Mass(0.1, 0.05, y, 0.05, 1));
                massPositions.push_back(Mass(0.1, 0.05, y, -0.05, 1));
                massPositions.push_back(Mass(0.1, -0.05, y, 0.05, 1));
                massPositions.push_back(Mass(0.1, -0.05, y, -0.05, 1));

                massPositions.push_back(Mass(0.1, 0.05, y, z + 0.1, 1));
                massPositions.push_back(Mass(0.1, 0.05, y, z, 1));
                massPositions.push_back(Mass(0.1, -0.05, y, z + 0.1, 1));
                massPositions.push_back(Mass(0.1, -0.05, y, z, 1));

                massPositions.push_back(Mass(0.1, x, y, 0.05, 0));
                massPositions.push_back(Mass(0.1, x, y, -0.05, 0));
                massPositions.push_back(Mass(0.1, x + 0.1, y, 0.05, 0));
                massPositions.push_back(Mass(0.1, x + 0.1, y, -0.05, 0));

                massPositions.push_back(Mass(0.1, x, y, z + 0.1, 0));
                massPositions.push_back(Mass(0.1, x, y, z, 0));
                massPositions.push_back(Mass(0.1, x + 0.1, y, z + 0.1, 0));
                massPositions.push_back(Mass(0.1, x + 0.1, y, z, 0));
            }

            /*for (int a = h; a >= 0; a--) {
                y = BASE + height + edgeLength * a;
                x = BASE;
                for (int i = 0; i < m * (w + 1); i += m) {
                    z = BASE;
                    for (int j = 0; j < m; j++) {
                        massPositions.push_back(Mass(0.1, x, y, z, 0));
                        z += 0.1;
                    }
                    x += 0.1;
                }
            }
            */
            std::vector<int> indexes;
            /*int i = 0;

            double dx, dy, dz;
            while (i < massPositions.size()) {
                for (int j = i+1; j < massPositions.size(); j++) {
                    dx = massPositions[i].x - massPositions[j].x;
                    dy = massPositions[i].y - massPositions[j].y;
                    dz = massPositions[i].z - massPositions[j].z;
                    if (sqrt(dx * dx + dy * dy + dz * dz) < edgeLength*2.0) {
                        springs.push_back(Spring(&(massPositions[i]), &(massPositions[j])));
                    }
                }
               i += 1;
            }*/

            int i = 0;
            int skip = 0;
            double dx, dy, dz;

            while (i < masses - m) {
                indexes = { i, i + 1, i + 2, i + 3, i + m, i + 1 + m,i + 2 + m, i + 3 + m };
                for (int j = 0; j < indexes.size(); j++) {
                    for (int k = 0; k < j; k++) {
                        springs.push_back(Spring(&(massPositions[indexes[j]]), &(massPositions[indexes[k]])));
                    }
                }
                skip++;
                if (skip == 3) {
                    skip = 0;
                    i += 2;
                }
                i += 2;
            }
            skip = 0;
            i = masses;
            while (i < massPositions.size()) {
                if (skip < 4) {
                    if (skip == 0) {
                        indexes = { i, i + 1, i + 2, i + 3, 1, 3, 9, 11 };
                    }
                    if (skip == 1) {
                        indexes = { i, i + 1, i + 2, i + 3, 5, 7, 13, 15 };
                    }
                    if (skip == 2) {
                        indexes = { i, i + 1, i + 2, i + 3, 17, 19, 25, 27 };
                    }
                    if (skip == 3) {
                        indexes = { i, i + 1, i + 2, i + 3, 21, 23, 29, 31 };
                    }
                }
                else {
                    indexes = { i, i + 1, i + 2, i + 3, i - 16, i - 15, i - 14, i - 13 };
                }
                for (int j = 0; j < indexes.size(); j++) {
                    for (int k = 0; k < j; k++) {
                        springs.push_back(Spring(&(massPositions[indexes[j]]), &(massPositions[indexes[k]])));
                    }
                }
                skip++;
                i += 4;
            }
            
            double bestFit, group, newFit;
            double perimeter = edgeLength;
            for (i = 0; i < massPositions.size(); i++) {
                //Determine encoding
                bestFit = 10000000;
                group = 1;
                x = massPositions[i].x0;
                y = massPositions[i].y0;
                z = massPositions[i].z0;
                for (double j = 0; j < encoding.size(); j += 4) {
                    dx = x - encoding[j + 1];
                    dy = y - encoding[j + 2];
                    dz = z - encoding[j + 3];
                    newFit = sqrt(dx * dx + dy * dy + dz * dz);
                    if (newFit < bestFit) {//&&(encoding->encoding[j]!=EMPTY||(abs(dx)<=perimeter&&abs(dy)<=perimeter&&abs(dz)<=perimeter))) {
                        bestFit = newFit;
                        group = encoding[j];
                    }
                }
                massPositions[i].updateGroup(group);
            }

            for (i = 0; i < springs.size(); i++) {
                springs[i].reset();
            }

        }
        void draw(unsigned int modelLoc, unsigned int colorLocation, double xoffset, double yoffset, double zoffset) {
            glm::mat4 model;
            
            for (int i = 0; i < massPositions.size(); i++) {
                if (massPositions[i].group!=EMPTY) {// create transformations

                    massPositions[i].finish();
                    model = glm::mat4(1.0f);
                    //model = glm::rotate(model, (float)glfwGetTime(), glm::vec3(0.0f, 1.0f, 0.0f));
                    
                    model = glm::translate(model, glm::vec3(massPositions[i].x+xoffset, massPositions[i].y+yoffset, massPositions[i].z+zoffset));

                    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
                    if (massPositions[i].group == 0) { glUniform4f(colorLocation, 0.36389f, 0.618501f, 0.782349f, 1.0f); }
                    if (massPositions[i].group == 1) { glUniform4f(colorLocation, 0.8, 0.218501f, 0.582349f, 1.0f); }
                    if (massPositions[i].group == 2) { glUniform4f(colorLocation, 0.36389f, 0.818501f, 0.282349f, 1.0f); }
                    if (massPositions[i].group == 3) { glUniform4f(colorLocation, 0.16389f, 0.118501f, 0.782349f, 1.0f); }

                    sphere.draw();
                }

            }
            for (int i = 0; i < springs.size(); i++) {
                if (springs[i].m1->group!=EMPTY&&springs[i].m2->group!=EMPTY) {
                    model = glm::mat4(1.0f);
                    //model = glm::rotate(model, (float)glfwGetTime(), glm::vec3(0.0f, 1.0f, 0.0f));

                    model = glm::translate(model, glm::vec3(springs[i].m1->x+xoffset, springs[i].m1->y+yoffset, springs[i].m1->z+zoffset));

                    if ((springs[i].x != 0.0) || (springs[i].z != 0.0)) {
                        model = glm::rotate(model, glm::acos((float)-springs[i].y), glm::vec3(-springs[i].z, 0.0f, springs[i].x));
                    }
                    else {
                        model = glm::rotate(model, glm::acos((float)-springs[i].y), glm::vec3(0.0, 0.0f, 1.0f));
                    }
                    model = glm::scale(model, glm::vec3(1.0, (springs[i].n / edgeLength), 1.0));

                    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

                    if (springs[i].m1->group == 0) { glUniform4f(colorLocation, 0.36389f, 0.618501f, 0.782349f, 1.0f); }
                    if (springs[i].m1->group == 1) { glUniform4f(colorLocation, 0.8, 0.218501f, 0.582349f, 1.0f); }
                    if (springs[i].m1->group == 2) { glUniform4f(colorLocation, 0.36389f, 0.818501f, 0.282349f, 1.0f); }
                    if (springs[i].m1->group == 3) { glUniform4f(colorLocation, 0.16389f, 0.118501f, 0.782349f, 1.0f); }

                    rod.draw();
                    springs[i].updateL();
                    springs[i].calcForce();
                }
            }

        }
    

};



void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}


int main()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    GLFWwindow* window = glfwCreateWindow(1600, 1200, "Bouncing Cube", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    // callback func on resize
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glEnable(GL_MULTISAMPLE);


    /* std::ofstream PEFile;
    PEFile.open("PE.txt");
    std::ofstream KEFile;
    KEFile.open("KE.txt");
    std::ofstream TEFile;
    TEFile.open("TE.txt");
    std::ofstream DFile;
    DFile.open("Data.txt");*/

    // VBO stores data on GPU
    // benefit in that can send large
    // amount of data to GPU at once
    Shader ourShader("shader.vs", "shader.fs");

    ourShader.use();
    Platform platform;
    unsigned int colorLocation = glGetUniformLocation(ourShader.ID, "ourColor");
    glUniform4f(colorLocation, 0.36389f, 0.618501f, 0.782349f, 1.0f);
    // create transformations
    glm::mat4 model = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 view = glm::mat4(1.0f);
    glm::mat4 projection = glm::mat4(1.0f);
    view = glm::translate(view, glm::vec3(0.0f,0.0f, -4.0f));
    view = glm::rotate(view, glm::radians(10.0f), glm::vec3(1.0f, 1.0f, 0.0f));

    projection = glm::perspective(glm::radians(45.0f), (float)(1600) / (float)(1200), 0.1f, 100.0f);
    // retrieve the matrix uniform locations
    unsigned int modelLoc = glGetUniformLocation(ourShader.ID, "model");
    unsigned int viewLoc = glGetUniformLocation(ourShader.ID, "view");
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, &view[0][0]);
    ourShader.setMat4("projection", projection);


    //std::vector<double> encoding1 = { 0, -0.0483292, -0.174336, 0.0861918, 1, -0.203918, -0.00920917, -0.290548, 2, 0.064522, -0.0182399, -0.198711, 0, 0.152264, -0.166984, -0.0935623, 4, -0.327711, -0.340461, 0.00181126, 0, -0.447231, 0.101568, -0.251431,4, -0.133645, 0.0185936, -0.0928222, 0, -0.0965951, 0.0463574, -0.0422594, 1, -0.0670328, 0.0548142, -0.198667, };
    // run1 0.0194
    //std::vector<double> encoding1 = { 1, -0.267494, -0.158964, 0.176591, 3, 0.0453366, -0.38408, -0.0196801, 2, -0.383354, 0.169509, 0.0990236, 1, -0.0263234, -0.218835, 0.16183, 2, -0.146727, -0.00405902, -0.0550172, 3, -0.387128, -0.11097, -0.00650521, 1, -0.079937, -0.186604, 0.111722, 3, -0.30113, 0.00274039, -0.154454, 4, -0.356223, -0.127862, -0.0901477, };
    
    std::vector<double> encoding1 = { 0, -0.0620439, -0.204137 + height, -0.0132993, 2, -0.24203, -0.152991 + height, -0.271615, 3, 1.19569e-05, -0.11549 + height, -0.382024, 1, 0.161176, 0.149245 + height, -0.338615, 1, 0.091754, -0.050655 + height, -0.344309, 4, -0.214686, 0.186436 + height, -0.106911, 4, -0.11491, -0.244058 + height, -0.0336033, 4, -0.270455, -0.280325 + height, -0.0900362, 0, 0.0341074, -0.249689 + height, 0.154513, };
    // instance 2 0.0254
    std::vector<double> encoding2 = { 2, -0.017768, 0.0374664 + height, -0.106112, 1, -0.185453, -0.268189 + height, 0.119048, 2, -0.397952, -0.0768364 + height, -0.22755, 3, -0.000154576, 0.116368 + height, 0.0991875, 1, 0.166903, -0.130402 + height, 0.191266, 0, -0.366917, -0.00770362 + height, 0.179723, 3, -0.100101, 0.182547 + height, -0.0593867, 2, -0.304727, -0.0595171 + height, -0.37006, 4, -0.0784924, -0.237014 + height, 0.052112, };
    //std::vector<double> encoding3 = { 2, 0.0309603, -0.106138, -0.356798, 1, 0.170042, 0.142722, 0.0083635, 0, -0.25831, -0.0742616, -0.29498, 3, 0.182249, -0.31517, 0.0578691, 0, 0.148946, -0.301853, 0.0563522, 1, -0.0179736, 0.145622, 0.185604, 4, -0.040763, -0.258199, -0.1144, 2, 0.19824, -0.326341, -0.189295, 1, 0.00951487, -0.283962, -0.142605, };
    // instance 1 0.0215
    std::vector<double> encoding3 = { 4, -0.207362, -0.292036+height, -0.252469, 1, 0.118805, -0.239553 + height, -0.110972, 2, -0.170515, -0.0986156 + height, -0.185227, 1, -0.382845, 0.206599 + height, -0.361608, 2, -0.293663, -0.102109 + height, -0.0689958, 1, -0.302256, 0.0206676 + height, 0.0664337, 1,-0.0916627, -0.394442 + height, -0.14061, 4, -0.192837, -0.236811 + height, -0.0827886, 4, -0.30359, -0.263236 + height, -0.293268, };
    //std::vector<double> encoding = { 2, -0.293228, 0.00408948 + height, -0.30795, 1, -0.341533, -0.308866 + height, 0.00399792, 2, 0.171215, 0.110843 + height, 0.0296863, 1, 0.14406, 0.196228 + height, -0.391137 };
    std::vector<double> encoding = { 1, 0.17911, -0.0289987 + height, -0.166724, 3, -0.233963 , -0.208301 + height, 0.026191, 1, -0.117753, -0.261213+height, 0.0163396, 2, 0.154416, -0.0125401+height, -0.0380078, 2, -0.133061, -0.030427+height, 0.112363, 2, -0.180706, 0.0765206+height, -0.225935, };

    
    
    double PE;
    double KE;
    double TE;
    double time = 0;
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    Robot robot1(3,3,3,encoding1);
    Robot robot2(3, 3, 3, encoding2);
    Robot robot3(3, 3, 3, encoding);

    while (time < 10.0 && !glfwWindowShouldClose(window))
    {
        // Input
        processInput(window);

        // Render
        // glClearColor sets color to appear when
        // color buffer cleared
        // Set state
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);


        // Use state
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ourShader.use();


        //robot1.draw(modelLoc,colorLocation,-0.0,0.0,-1.5);
        //robot2.draw(modelLoc, colorLocation, -0.0, 0.0, -0.0);
        robot3.draw(modelLoc, colorLocation, 0.0, 0.0, 0.5);
        platform.draw(modelLoc, colorLocation);
        time += DT;
                
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    glfwTerminate();
    //PEFile.close();
    //KEFile.close();
    //TEFile.close();
    //DFile.close();
    return 0;
}
