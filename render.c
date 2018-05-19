#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "svpng.inc"

#define TWO_PI 6.28318530718f
#define H 512
#define W 512
#define N 64
#define MAX_STEP 64
#define MAX_DISTANCE 2.0f
#define EPSILON 1e-6f
#define DELTA 1e-4

typedef struct {
  float sd;        /* signed distance bwteen sample point and scene */
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
float planeSDF(float x, float y, float px, float py, float nx, float ny);
float segmentSDF(float x, float y, float ax, float ay, float bx, float by);
float capsuleSDF(float x, float y, float ax, float ay, float bx, float by, float r);
float triangleSDF();

/* Combined geometry */
obj trianglemesh(float x, float y, float cx, float cy, float l, float e, int dir);
obj sierpinski(float x, float y, float cx, float cy, float l, float e);

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
  // obj r = a;
  // r.sd = a.sd > -b.sd ? a.sd : -b.sd;
  // return r;
  return csg_intersect(a, csg_complement(b));
}

float circleSDF(float x, float y, float cx, float cy, float r) {
  float ux = x - cx, uy = y - cy;
  return sqrtf(ux * ux + uy * uy) - r;
}

float planeSDF(float x, float y, float px, float py, float nx, float ny) {
  return (x - px) * nx + (y - py) * ny;
}

float segmentSDF(float x, float y, float ax, float ay, float bx, float by) {
  float vx = x - ax, vy = y - ay, ux = bx - ax, uy = by - ay;
  float t = fmaxf(fminf((vx * ux + vy * uy) / (ux * ux + uy * uy), 1.0f), 0.0f);
  float dx = vx - t * ux, dy = vy - t * uy;
  return sqrtf(dx * dx + dy * dy);
}

float capsuleSDF(float x, float y, float ax, float ay, float bx, float by, float r) {
  return segmentSDF(x, y, ax, ay, bx, by) - r;
}

obj trianglemesh(float x, float y, float cx, float cy, float l, float e, int dir) {
  float ux = cx, uy = cy - dir * 0.57735026918f * l;
  float vx = cx - 0.5f * l, vy = cy + dir * 0.28867513459f * l;
  float wx = cx + 0.5f * l, wy = vy;
  obj seg1 = { segmentSDF(x, y, ux, uy, vx, vy), e };
  obj seg2 = { segmentSDF(x, y, vx, vy, wx, wy), e };
  obj seg3 = { segmentSDF(x, y, wx, wy, ux, uy), e };
  obj seg4 = { segmentSDF(x, y, ux, uy + dir * DELTA, vx + DELTA * 0.86602540378f, vy - dir * DELTA * 0.5f), 0.01f };
  obj seg5 = { segmentSDF(x, y, vx + DELTA * 0.86602540378f, vy - dir * DELTA * 0.5f, wx - DELTA * 0.86602540378f, wy - dir * DELTA * 0.5f), 0.01f };
  obj seg6 = { segmentSDF(x, y, wx - DELTA * 0.86602540378f, wy - dir * DELTA * 0.5f, ux, uy + dir * DELTA), 0.01f };
  return csg_union(csg_union(csg_union(seg1, seg4), csg_union(seg2, seg5)), csg_union(seg3, seg6));
}

obj sierpinski(float x, float y, float cx, float cy, float l, float e) {
  float cx1   = cx,              cy1   = cy;
  float cx11  = cx,              cy11  = cy - 0.28867513459 * l;
  float cx12  = cx - l / 4.0f,   cy12  = cy + 0.14433756729f * l;
  float cx13  = cx + l / 4.0f,   cy13  = cy + 0.14433756729f * l;
  float cx111 = cx11,            cy111 = cy11 - 0.14433756729f * l;
  float cx112 = cx11 - l / 8.0f, cy112 = cy11 + 0.07216878364f * l;
  float cx113 = cx11 + l / 8.0f, cy113 = cy11 + 0.07216878364f * l;
  float cx121 = cx12,            cy121 = cy12 - 0.14433756729f * l;
  float cx122 = cx12 - l / 8.0f, cy122 = cy12 + 0.07216878364f * l;
  float cx123 = cx12 + l / 8.0f, cy123 = cy12 + 0.07216878364f * l;
  float cx131 = cx13,            cy131 = cy13 - 0.14433756729f * l;
  float cx132 = cx13 - l / 8.0f, cy132 = cy13 + 0.07216878364f * l;
  float cx133 = cx13 + l / 8.0f, cy133 = cy13 + 0.07216878364f * l;
  obj trimesh0   = trianglemesh(x, y, cx, cy, l, e, 1);
  obj trimesh1   = trianglemesh(x, y, cx1, cy1, l / 2.0f, e, -1);
  obj trimesh11  = trianglemesh(x, y, cx11, cy11, l / 4.0f, e, -1);
  obj trimesh12  = trianglemesh(x, y, cx12, cy12, l / 4.0f, e, -1);
  obj trimesh13  = trianglemesh(x, y, cx13, cy13, l / 4.0f, e, -1);
  obj trimesh111 = trianglemesh(x, y, cx111, cy111, l / 8.0f, e, -1);
  obj trimesh112 = trianglemesh(x, y, cx112, cy112, l / 8.0f, e, -1);
  obj trimesh113 = trianglemesh(x, y, cx113, cy113, l / 8.0f, e, -1);
  obj trimesh121 = trianglemesh(x, y, cx121, cy121, l / 8.0f, e, -1);
  obj trimesh122 = trianglemesh(x, y, cx122, cy122, l / 8.0f, e, -1);
  obj trimesh123 = trianglemesh(x, y, cx123, cy123, l / 8.0f, e, -1);
  obj trimesh131 = trianglemesh(x, y, cx131, cy131, l / 8.0f, e, -1);
  obj trimesh132 = trianglemesh(x, y, cx132, cy132, l / 8.0f, e, -1);
  obj trimesh133 = trianglemesh(x, y, cx133, cy133, l / 8.0f, e, -1);
  obj comb1 = csg_union(trimesh0, trimesh1);
  obj comb2 = csg_union(trimesh11, trimesh12);
  obj comb3 = csg_union(trimesh13, trimesh111);
  obj comb4 = csg_union(trimesh112, trimesh113);
  obj comb5 = csg_union(trimesh121, trimesh122);
  obj comb6 = csg_union(trimesh123, trimesh131);
  obj comb7 = csg_union(trimesh132, trimesh133);
  obj comb12 = csg_union(comb1, comb2);
  obj comb34 = csg_union(comb3, comb4);
  obj comb56 = csg_union(comb5, comb6);
  return csg_union(csg_union(comb12, comb34), csg_union(comb56, comb7));
}

obj scene(float x, float y) {
  // obj seg = { capsuleSDF(x, y, 0.4f, 0.4f, 0.6f, 0.6f, 0.1f), 1.0f };
  // obj a = { circleSDF(x, y, 0.7f, 0.7f, 0.2f), 2.0f };
  // obj b = { circleSDF(x, y, 0.65f, 0.65f, 0.05f), 0.5f };
  // obj c = { capsuleSDF(x, y, 0.04f, 0.43f, 0.52f, 0.08f, 0.002f), 2.0f };
  // obj d = { circleSDF(x, y, 0.5f, 0.5f, 0.05f), 0.0f };
  // return csg_union(csg_subtract(a, b), csg_union(c, d));
  return sierpinski(x, y, 0.5f, 0.5f, 0.4f, 2.0f);
  // return trianglemesh(x, y, 0.5f, 0.5f, 0.4f, 2.0f, 1);
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
  svpng(fopen("foo_1.png", "wb"), W, H, img, 0);
  return 0;
}

