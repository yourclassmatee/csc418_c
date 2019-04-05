attribute vec3 position; // Given vertex position in object space
attribute vec3 normal; // Given vertex normal in object space

uniform mat4 projection, modelview, normalMat; // Given scene transformation matrices

// These will be given to the fragment shader and interpolated automatically
varying vec3 normalInterp;
varying vec3 vertPos;
varying vec4 color;

uniform float Ka;   // Ambient reflection coefficient
uniform float Kd;   // Diffuse reflection coefficient
uniform float Ks;   // Specular reflection coefficient
uniform float shininessVal; // Shininess

// Material color
uniform vec3 ambientColor;
uniform vec3 diffuseColor;
uniform vec3 specularColor;
uniform vec3 lightPos; // Light position


void main(){
  // Your solution should go here.
  // Only the ambient colour calculations have been provided as an example.

  //ambient
  vec4 a_color = vec4(ambientColor * Ka, 1.0); 

  //diffuse
  vec4 vertPos4 = modelview * vec4(position, 1.0);
  gl_Position = projection * vertPos4;
  vec4 lightPos4 =  vec4(lightPos, 1.0);
  vec4 light = normalize(lightPos4 - vertPos4);
  vec4 normalVec4 = normalMat * vec4(normal, 1.0);
  float d_intensity = dot(light, normalVec4) * Kd;
  if (d_intensity < 0.0){
    d_intensity = 0.0;
  }
  vec4 d_color = vec4(diffuseColor * d_intensity, 1.0);

  //specular
  vec4 eye_vec = normalize(vec4(0.0,0.0,0.0,1.0) - vertPos4);
  vec4 h_vec = normalize(light + eye_vec);
  float s_intensity = dot(h_vec, normalVec4);
  if(s_intensity <= 0.0){
    s_intensity = 0.0;
  }
  else{
    s_intensity = pow(s_intensity, shininessVal) * Ks;
  }
  vec4 s_color = vec4(specularColor * s_intensity, 1.0);

  color = a_color + d_color +s_color;

}