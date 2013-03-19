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
#include "math_3d.h"
#include "Point.h"
#include "Triangle.h"
#include "delaunay/DelaunayTriangulation.h"

GLuint VBO;
GLuint gWorldLocation;

static const char* pVS = "                                                          \n\
#version 330                                                                        \n\
                                                                                    \n\
layout (location = 0) in vec3 Position;                                             \n\
                                                                                    \n\
uniform mat4 gWorld;                                                                \n\
                                                                                    \n\
void main()                                                                         \n\
{                                                                                   \n\
    gl_Position = gWorld * vec4(Position, 1.0);                                     \n\
}";

static const char* pFS = "                                                          \n\
#version 330                                                                        \n\
                                                                                    \n\
out vec4 FragColor;                                                                 \n\
                                                                                    \n\
void main()                                                                         \n\
{                                                                                   \n\
    FragColor = vec4(1.0, 0.0, 0.0, 1.0);                                           \n\
}";

static void renderSceneCB()
{
        glClear (GL_COLOR_BUFFER_BIT);

        static float scale = 0.0f;

        scale += 0.001f;

        Matrix4f world;

        world.m[0][0] = sinf(scale) ; world.m[0][1] = 0.0f       ; world.m[0][2] = 0.0f;        world.m[0][3] = 0.0f;
        world.m[1][0] = 0.0f        ; world.m[1][1] = sinf(scale); world.m[1][2] = 0.0f;        world.m[1][3] = 0.0f;
        world.m[2][0] = 0.0f;       ; world.m[2][1] = 0.0f;      ; world.m[2][2] = sinf(scale); world.m[2][3] = 0.0f;
        world.m[3][0] = 0.0f;       ; world.m[3][1] = 0.0f;      ; world.m[3][2] = 0.0f;        world.m[3][3] = 1.0f;

        glUniformMatrix4fv (gWorldLocation, 1, GL_TRUE, &world.m[0][0]);

        glEnableVertexAttribArray (0);
        glBindBuffer (GL_ARRAY_BUFFER, VBO);
        glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, 0, 0);

        glDrawArrays (GL_TRIANGLES, 0, 3);

        glDisableVertexAttribArray (0);

        glutSwapBuffers ();
}

static void initializeGlutCallbacks ()
{
        glutDisplayFunc (renderSceneCB);
        glutIdleFunc (renderSceneCB);
}

static void createVertexBuffer ()
{
        Vector3f vertices[3];
        vertices[0] = Vector3f (-1.0f, -1.0f, 0.0f);
        vertices[1] = Vector3f (1.0f, -1.0f, 0.0f);
        vertices[2] = Vector3f (0.0f, 1.0f, 0.0f);

        glGenBuffers (1, &VBO);
        glBindBuffer (GL_ARRAY_BUFFER, VBO);
        glBufferData (GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
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
}

int main (int argc, char** argv)
{
        glutInit (&argc, argv);
        glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA);
        glutInitWindowSize (1024, 768);
        glutInitWindowPosition (100, 100);
        glutCreateWindow ("Tutorial 08");

        initializeGlutCallbacks ();

        // Must be done after glut is initialized!
        GLenum res = glewInit ();
        if (res != GLEW_OK) {
                fprintf (stderr, "Error: '%s'\n", glewGetErrorString (res));
                return 1;
        }

        glClearColor (0.0f, 0.0f, 0.0f, 0.0f);

        createVertexBuffer ();

        compileShaders ();

        glutMainLoop ();

        return 0;
}
