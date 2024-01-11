#ifndef GLPROGRAM_H
#define GLPROGRAM_H

#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "Shader.h"
#include "Camera.h"

#include "../include/tinynurbs/tinynurbs.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#define MIN(a, b) ((a) < (b)) ? (a) : (b)

#include "SurfaceRender.h"
#include "CurveRender.h"
#include "MeshRender.h"
#include "NormalRender.h"




class GLProgram {
    private:
        GLFWwindow* window;

        float deltaTime, prevTime;

        struct color {
            float r = 0.0f;
            float g = 0.0f;
            float b = 0.0f;
            float alpha = 1.0f;
        } clearColor;

        static glm::vec3 getArcballVector(float x, float y); // helper to cursor callback, (x,y) are raw mouse coordinates



        vector<CurveRender> curveRender;  //���ߵ�render
        vector<CurveRender> curveControlPointLineRender; //���ߵĿ��Ƶ�����
        vector<NormalRender> curveDerivateRender;  //���ߵĵ���������

        
        vector<SurfaceRender> surfaceRender; //�����render
        vector<MeshRender> meshRender; // ����ĵȲ���������
        vector<MeshRender> surfaceControlRender;//����Ŀ��Ƶ�����
        vector<NormalRender> normalRender;  //����ķ���
        vector<NormalRender> UsurfaceDerivateRender; //����ĳһ�㵼��U�����������
        vector<NormalRender> VsurfaceDerivateRender; //����ĳһ�㵼��V�����������




    public:

        // input
        bool enableControl = true; //��������ƶ�

        static int windowWidth, windowHeight;
        static Camera camera;
        static bool mousePressed;
        static double prevMouseX, prevMouseY;
        static glm::mat4 modelMatrix;

        // lighting
        glm::vec3 lightPos;

        GLProgram();

        void init(vector<tinynurbs::RationalCurve<double>> curves, vector<tinynurbs::RationalSurface<double>> surfaces);
        void run(vector<tinynurbs::RationalCurve<double>> curves, vector<tinynurbs::RationalSurface<double>> surfaces);
        void cleanup(void);

        void setClearColor(float r, float g, float b, float alpha);

        // event callback functions
        void framebufferSizeCallback(GLFWwindow* window, int width, int height);
        void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);
        void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
        void cursorPosCallback(GLFWwindow* window, double xpos, double ypos);

        // input
        void processInput(void);


private:
    // Wrapper for glfwGetWindowUserPointer() returning this class instance.
    static GLProgram* _this(GLFWwindow* window)
    {
        return static_cast<GLProgram*>(glfwGetWindowUserPointer(window));
    }
    static void _framebufferSizeCallback(GLFWwindow* window, int width, int height)
    {
        _this(window)->framebufferSizeCallback(window, width, height);
    }
    static void _scrollCallback(GLFWwindow* window, double xoffset, double yoffset)
    {
        _this(window)->scrollCallback(window, xoffset, yoffset);
    }
    static void _mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
    {
        _this(window)->mouseButtonCallback(window, button, action, mods);
    }
    static void _cursorPosCallback(GLFWwindow* window, double xpos, double ypos)
    {
        _this(window)->cursorPosCallback(window, xpos, ypos);
    }

};

#endif //GLPROGRAM_H
