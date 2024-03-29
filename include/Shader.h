#ifndef SHADER_H
#define SHADER_H

#include <glad/glad.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

typedef unsigned int			uint;

class Shader {
    public:
        uint ID;
        Shader();
        Shader(const char* vertexPath, const char* fragmentPath);
        void use(void);
        void setFloatUniform(const std::string &name, float value) const;
        void setVec3Uniform(const std::string &name, glm::vec3 value) const;
        void setMat4Uniform(const std::string &name, glm::mat4 value) const;

    private:
        void checkCompileErrors(uint shader, std::string type);
};

#endif //SHADER_H
