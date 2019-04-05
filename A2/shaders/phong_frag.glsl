precision mediump float; // It is required to set a floating point precision in all fragment shaders.

// Interpolated values from vertex shader
varying vec3 normalInterp; // Normal
varying vec3 vertPos; // Vertex position

uniform float Ka;   // Ambient reflection coefficient
uniform float Kd;   // Diffuse reflection coefficient
uniform float Ks;   // Specular reflection coefficient
uniform float shininessVal; // Shininess

// Material color
uniform vec3 ambientColor;
uniform vec3 diffuseColor;
uniform vec3 specularColor;
uniform vec3 lightPos; // Light position

void main() {
  // Your solution should go here.
  // Only the ambient colour calculations have been provided as an example.

  //ambient
  vec4 a_color = vec4(ambientColor * Ka, 1.0); 
  
  //diffuse
  vec3 light = normalize(lightPos - vertPos);
  float d_intensity = dot(light, normalInterp) * Kd;
  if (d_intensity < 0.0){
    d_intensity = 0.0;
  }
  vec4 d_color = vec4(diffuseColor * d_intensity, 1.0);

  //specular
  vec3 eye_vec = normalize(vec3(0.0,0.0,0.0) - vertPos);
  vec3 h_vec = normalize(light + eye_vec);
  float s_intensity = dot(h_vec, normalInterp);
  if(s_intensity <= 0.0){
    s_intensity = 0.0;
  }
  else{
    s_intensity = pow(s_intensity, shininessVal) * Ks;
  }
  vec4 s_color = vec4(specularColor * s_intensity, 1.0);

  gl_FragColor = a_color + d_color + s_color;
}
