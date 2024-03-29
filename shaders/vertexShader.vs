#version 450 core

layout (location = 0) in vec3 pos;

out vec3 fragPos;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main() {


    gl_Position = projection * view * model * vec4(pos, 1.0);
   // fragPos = pos;

    fragPos = vec3(model * vec4(pos, 1.0));
}
