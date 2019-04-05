attribute vec3 position; // Given vertex position in object space
attribute vec3 normal; // Given vertex normal in object space

uniform mat4 projection, modelview, normalMat; // Given scene transformation matrices

// These will be given to the fragment shader and interpolated automatically
varying vec3 normalInterp; // Normal
varying vec3 vertPos; // Vertex position


void main(){
  // Your solution should go here.
  // Only the ambient colour calculations have been provided as an example.
  vec4 vertPos4 = modelview * vec4(position, 1.0);
  gl_Position = projection * vertPos4;

  vec4 normal4 = normalMat * vec4(normal,1.0);
  normalInterp = normal4.xyz;
  vertPos = vertPos4.xyz;

  
}
