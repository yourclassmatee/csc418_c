precision mediump float; // It is required to set a floating point precision in all fragment shaders.

// Interpolated values from vertex shader
varying vec3 normalInterp; // Normal
varying vec3 vertPos; // Vertex position
varying vec3 viewVec; // View vector (eye to fragment)

uniform float Ka;   // Ambient reflection coefficient
uniform float Kd;   // Diffuse reflection coefficient
uniform float Ks;   // Specular reflection coefficient
uniform float shininessVal; // Shininess

// Material color
// HINT: Use the environment map as the ambient color
uniform vec3 diffuseColor;
uniform vec3 specularColor;
uniform vec3 lightPos; // Light position

uniform samplerCube envTexSampler; // A GLSL sampler represents a single texture. A samplerCube can be used to sample a cubemap texture.

void main() {
  // Your solution should go here.
  
  // The model is currently rendered in black
  //gl_FragColor = vec4(vec3(0.0), 1.0);

  vec3 cubemap_vec = reflect(-viewVec, normalInterp);
  vec4 a_color = textureCube(envTexSampler, cubemap_vec);

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

  //color = a_color + d_color +s_color;

//gl_FragColor = vec4(ambientColor, 1.0);
  gl_FragColor = a_color + d_color + s_color;



}
