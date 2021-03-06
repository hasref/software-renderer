# software-renderer

A raytracer and "renderer" written in C++, (mostly) from scratch using Gabriel Gambetta's ["Computer Graphics from Scratch"](https://gabrielgambetta.com/computer-graphics-from-scratch/).

## Images

### Obligatory Reflective Spheres

![scene with reflective spheres](images/reflective.png)

## Building

To build this project, you will need `conan` and `cmake`. While I'm using std=c++20, I do not anticipate needing any features from the standard - this should probably also compile with c++17.

```bash
mkdir build && cd build
conan install .. --build=missing
cmake .. && make
```
