// Fragment shader template for the bonus question

precision mediump float; // It is required to set a floating point precision in all fragment shaders.

// Interpolated values from vertex shader
// NOTE: You may need to edit this section to add additional variables
varying vec3 normalInterp; // Normal
varying vec3 vertPos; // Vertex position
varying vec3 viewVec; // Interpolated view vector

// uniform values remain the same across the scene
// NOTE: You may need to edit this section to add additional variables
uniform float Ka;   // Ambient reflection coefficient
uniform float Kd;   // Diffuse reflection coefficient
uniform float Ks;   // Specular reflection coefficient
uniform float shininessVal; // Shininess

// Material color
uniform vec3 ambientColor;
uniform vec3 diffuseColor;
uniform vec3 specularColor;

uniform vec3 lightPos; // Light position

uniform sampler2D uSampler;	// 2D sampler for the earth texture
uniform samplerCube envTexSampler;	// cube sampler for the environment map

void main() {
  // Your solution should go here.
  // Only the ambient colour calculations have been provided as an example.
  //gl_FragColor = vec4(ambientColor, 1.0);

  //calculating opacity
  float opacity = dot(normalize(normalInterp), normalize(-viewVec));
  opacity = abs(opacity);
  opacity = 1.0-pow(opacity, 2.0);

  //ambient color
  vec4 a_color = vec4(ambientColor*opacity, 1.0);

  //diffuse
  vec3 light = normalize(lightPos - vertPos);
  float d_intensity = dot(light, normalInterp) * Kd;
  if (d_intensity < 0.0){
    d_intensity = 0.0;
  }
  //add diffuse color so that edges are thicker 
  //comparing to the middle section
  vec4 d_color = vec4(diffuseColor * d_intensity * opacity, 1.0);
  
  gl_FragColor = a_color + d_color;
  
}
