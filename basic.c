#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "svpng.inc"

#define TWO_PI 6.28318530718f
#define H 512
#define W 512
#define N 64
#define MAX_STEP 10
#define MAX_DISTANCE 2.0f
#define EPSILON 1e-6f

typedef struct {
  float sd;  /* signed distance bwteen sample point and scene */
  float emissive;  /* intensity of self-emission */
} obj;

unsigned char img[H * W * 3];

/**
 * Constructive Solid Geometry based on three operations
 *  UNION
 *  INTERSECTION
 *  RELATIVE COMPLEMENT
 */
obj csg_complement(obj a);
obj csg_union(obj a, obj b);
obj csg_intersect(obj a, obj b);
obj csg_subtract(obj a, obj b);

/* SDFs */
float circleSDF(float x, float y, float cx, float cy, float r);

/** Construct scene using geometry signed distance function. */
obj scene(float x, float y);

/**
 * Sphere tracing / Ray marching 
 *  \phi(x) > 0, x is out of scene, dist = \phi(x)
 *  \phi(x) < 0, x is inbound, dist = -\phi(x)
 *  \phi(x) = 0, x is on the border, dist = 0
 */
float trace(float ox, float oy, float dx, float dy);

/* Sample using Monte Carlo Integration */
float sample(float x, float y);


/* Definition */
obj csg_complement(obj a) {
  a.sd = -a.sd;
  return a;
}

obj csg_union(obj a, obj b) {
  return a.sd < b.sd ? a : b;
}

obj csg_intersect(obj a, obj b) {
  obj r = a.sd > b.sd ? b : a;
  r.sd = a.sd > b.sd ? a.sd : b.sd;
  return r;
}

obj csg_subtract(obj a, obj b) {
  // return csg_intersect(a, csg_complement(b));
  obj r = a;
  r.sd = a.sd > -b.sd ? a.sd : -b.sd;
  return r;
}

float circleSDF(float x, float y, float cx, float cy, float r) {
  float ux = x - cx, uy = y - cy;
  return sqrtf(ux * ux + uy * uy) - r;
}

obj scene(float x, float y) {
  obj c1 = { circleSDF(x, y, 0.4f, 0.5f, 0.2f), 1.0f };
  obj c2 = { circleSDF(x, y, 0.6f, 0.5f, 0.2f), 0.1f };
  return csg_subtract(c1, c2);
}

float trace(float ox, float oy, float dx, float dy) {
  float t = 0.001f;
  for (int i = 0; i < MAX_STEP && t < MAX_DISTANCE; i++) {
    obj o = scene(ox + t * dx, oy + t * dy);
    if (o.sd < EPSILON) return o.emissive;
    t += o.sd;
  }
  return 0.0f;
}

float sample(float x, float y) {
  float sum = 0.0f;
  for (int i = 0; i < N; i++) {
    float a = TWO_PI * (i + (float) rand() / RAND_MAX) / N;
    sum += trace(x, y, cosf(a), sinf(a));
  }
  return sum / N;
}

int main() {
  unsigned char *p = img;
  for (int y = 0; y < H; y++) 
    for (int x = 0; x < W; x++, p += 3) 
      p[0] = p[1] = p[2] = (int) fminf(sample((float) x / W, (float) y / H) * 255.0f, 255.0f);
  svpng(fopen("circles_1.png", "wb"), W, H, img, 0);
  return 0;
}

