/****************************************************************************
 *                                                                          *
 *  Author : lukasz.iwaszkiewicz@gmail.com                                  *
 *  Based on work by Etay Meiri                                             *
 *  ~~~~~~~~                                                                *
 *  License : GNU GPL. See COPYING file for details.                        *
 *  ~~~~~~~~~                                                               *
 ****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <fstream>
#include <iostream>
#include "Point.h"
#include "Triangle.h"
#include "TriangulationAdapter.h"

GLuint VBO;
GLuint IBO;

GLuint gWorldLocation;
GLuint colorLocation;

PointVector points;
MyTriagulation::TriangleVector const *triangulation;
MyTriagulation *cdt;

double scale;
Point translation;

static const char* pVS =
"#version 110                                                                       \n\
                                                                                    \n\
attribute vec2 Position;                                                            \n\
                                                                                    \n\
uniform mat4 gWorld;                                                                \n\
                                                                                    \n\
void main()                                                                         \n\
{                                                                                   \n\
    gl_Position = gWorld * vec4(Position, 0.0, 1.0);                                \n\
}";

static const char* pFS =
"#version 110                                                                       \n\
uniform vec4 color;                                                                 \n\
                                                                                    \n\
void main()                                                                         \n\
{                                                                                   \n\
    gl_FragColor = color;                                                           \n\
}";

struct Matrix4f
{
    float m[4][4];
};

static void renderSceneCB()
{
        glClear (GL_COLOR_BUFFER_BIT);

        Matrix4f world;

        world.m[0][0] = scale;        world.m[0][1] = 0.0f;       world.m[0][2] = 0.0f;        world.m[0][3] = translation.x;
        world.m[1][0] = 0.0f;         world.m[1][1] = scale;      world.m[1][2] = 0.0f;        world.m[1][3] = translation.y;
        world.m[2][0] = 0.0f;         world.m[2][1] = 0.0f;       world.m[2][2] = scale;       world.m[2][3] = 0.0f;
        world.m[3][0] = 0.0f;         world.m[3][1] = 0.0f;       world.m[3][2] = 0.0f;        world.m[3][3] = 1.0f;

        glUniformMatrix4fv (gWorldLocation, 1, GL_TRUE, &world.m[0][0]);

        glEnableVertexAttribArray (0);
        glBindBuffer (GL_ARRAY_BUFFER, VBO);
        glVertexAttribPointer (0, 2, GL_FLOAT, GL_FALSE, 0, 0);

        // Blue - indexed - triangles produced by triangulation.
        glUniform4f (colorLocation, 0, 0, 1, 0.4);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
        glDrawElements (GL_TRIANGLES, triangulation->size () * 3, GL_UNSIGNED_INT, 0);

        glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

        // Magenta - indexed - triangles produced by triangulation.
        glUniform4f (colorLocation, 1, 0, 1, 1);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
        glDrawElements (GL_TRIANGLES, triangulation->size () * 3, GL_UNSIGNED_INT, 0);

        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);


        // Red - input polygon segments.
        glUniform4f (colorLocation, 1, 0, 0, 1);
        glDrawArrays (GL_LINE_LOOP, 0, points.size ());

        // Green - input polygon points.
        glUniform4f (colorLocation, 0, 1, 0, 1);
        glPointSize (3);
        glDrawArrays (GL_POINTS, 0, points.size ());

        glDisableVertexAttribArray (0);

        glutSwapBuffers ();
}

static void initializeGlutCallbacks ()
{
        glutDisplayFunc (renderSceneCB);
        glutIdleFunc (renderSceneCB);
}

static void createVertexBuffer (const char *fileName)
{
        std::ifstream file (fileName);
        file >> scale;
        std::cerr << "Scale : " << scale << std::endl;

        file >> translation.x >> translation.y;
        std::cerr << "Translation : " << translation.x << "," << translation.y << std::endl;

        Point p;
        while (file >> p.x >> p.y) {
                points.push_back (p);
        }

/*--------------------------------------------------------------------------*/

        cdt = new MyTriagulation (points);
        cdt->constructDelaunay ();
        triangulation = &cdt->getTriangulation ();

/*--------------------------------------------------------------------------*/

        glGenBuffers (1, &VBO);
        glBindBuffer (GL_ARRAY_BUFFER, VBO);
        glBufferData (GL_ARRAY_BUFFER, points.size () * sizeof (Point), &points.front (), GL_STATIC_DRAW);

        glGenBuffers (1, &IBO);
        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, IBO);
        glBufferData (GL_ELEMENT_ARRAY_BUFFER, triangulation->size () * sizeof (Triangle), &triangulation->front (), GL_STATIC_DRAW);
}

static void addShader (GLuint ShaderProgram, const char* pShaderText, GLenum ShaderType)
{
        GLuint ShaderObj = glCreateShader (ShaderType);

        if (ShaderObj == 0) {
                fprintf (stderr, "Error creating shader type %d\n", ShaderType);
                exit (0);
        }

        const GLchar* p[1];
        p[0] = pShaderText;
        GLint Lengths[1];
        Lengths[0] = strlen (pShaderText);
        glShaderSource (ShaderObj, 1, p, Lengths);
        glCompileShader (ShaderObj);
        GLint success;
        glGetShaderiv (ShaderObj, GL_COMPILE_STATUS, &success);

        if (!success) {
                GLchar InfoLog[1024];
                glGetShaderInfoLog (ShaderObj, 1024, NULL, InfoLog);
                fprintf (stderr, "Error compiling shader type %d: '%s'\n", ShaderType, InfoLog);
                exit (1);
        }

        glAttachShader (ShaderProgram, ShaderObj);
}

static void compileShaders ()
{
        GLuint ShaderProgram = glCreateProgram ();

        if (ShaderProgram == 0) {
                fprintf (stderr, "Error creating shader program\n");
                exit (1);
        }

        addShader (ShaderProgram, pVS, GL_VERTEX_SHADER);
        addShader (ShaderProgram, pFS, GL_FRAGMENT_SHADER);

        GLint Success = 0;
        GLchar ErrorLog[1024] = { 0 };

        glLinkProgram (ShaderProgram);
        glGetProgramiv (ShaderProgram, GL_LINK_STATUS, &Success);

        if (Success == 0) {
                glGetProgramInfoLog (ShaderProgram, sizeof(ErrorLog), NULL, ErrorLog);
                fprintf (stderr, "Error linking shader program: '%s'\n", ErrorLog);
                exit (1);
        }

        glValidateProgram (ShaderProgram);
        glGetProgramiv (ShaderProgram, GL_VALIDATE_STATUS, &Success);

        if (!Success) {
                glGetProgramInfoLog (ShaderProgram, sizeof(ErrorLog), NULL, ErrorLog);
                fprintf (stderr, "Invalid shader program: '%s'\n", ErrorLog);
                exit (1);
        }

        glUseProgram (ShaderProgram);

        gWorldLocation = glGetUniformLocation (ShaderProgram, "gWorld");
        assert(gWorldLocation != 0xFFFFFFFF);

        colorLocation = glGetUniformLocation (ShaderProgram, "color");
        assert(colorLocation != 0xFFFFFFFF);
}



int main (int argc, char** argv)
{
        if (argc <= 1) {
                std::cerr << "Please provide an input file name." << std::endl;
                exit (1);
        }

        glutInit (&argc, argv);
        glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA);
//        glutInitWindowSize (1024, 768);
        glutInitWindowSize (320, 200);
        glutCreateWindow ("fikimiki test");

        glEnable (GL_CULL_FACE);
        glCullFace (GL_BACK);
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        initializeGlutCallbacks ();

        // Must be done after glut is initialized!
        GLenum res = glewInit ();
        if (res != GLEW_OK) {
                fprintf (stderr, "Error: '%s'\n", glewGetErrorString (res));
                return 1;
        }

        glClearColor (0.0f, 0.0f, 0.0f, 0.0f);

        createVertexBuffer (*(argv + 1));

        compileShaders ();

        glutMainLoop ();

        return 0;
}
