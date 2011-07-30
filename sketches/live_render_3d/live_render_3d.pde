import toxi.geom.Vec3D;
import processing.opengl.*;
import rui.marchingCubes.*;
 
// kinect
import org.openkinect.*;
import org.openkinect.processing.*;
 
MarchingCubes mc;
Vec3D rotationAxis;
 
Boolean bUseFill;
 
// kinect
Kinect kinect;
float a = 0;
// Size of kinect image
int w = 640;
int h = 480;
int kWidth  = 640;
int kHeight = 480;
// depth mapping and tilt
boolean depth = true;
boolean rgb = false;
boolean ir = false;
float deg = 8; // Start at 15 degrees
PImage depthImg;
int minDepth =  40;
int maxDepth = 860;
// set initial record to false
boolean record = false;
int counter = 0;
// print custom file
boolean printFile = false;
ArrayList points;
PrintWriter output;
// We'll use a lookup table so that we don't have to repeat the math over and over
float[] depthLookUp = new float[2048];
 
void setup(){
  size(1024, 600, OPENGL);
  Vec3D aabbMin = new Vec3D(-width/2, -height/2, -250);
  Vec3D aabbMax = new Vec3D(width/2, height/2, 250);
  Vec3D numPoints = new Vec3D(50,50,50);
  float isoLevel = 1;
  mc = new MarchingCubes(this, aabbMin, aabbMax, numPoints, isoLevel);
 
  rotationAxis = new Vec3D();
 
  bUseFill = false;
 
  // kinect
  kinect = new Kinect(this);
  kinect.start();
  kinect.enableDepth(true);
  kinect.tilt(deg);
  // We don't need the grayscale image in this example
  // so this makes it more efficient
  kinect.processDepthImage(false);
  // get depthImg to constrain
  depthImg = new PImage(kWidth, kHeight);
  // Lookup table for all possible depth values (0 - 2047)
  for (int i = 0; i < depthLookUp.length; i++) {
    depthLookUp[i] = rawDepthToMeters(i);
  }
  points = new ArrayList();
  output = createWriter("points.txt");
 
}
 
void draw(){
  background(255);
  lights();
 
  // kinect
  int[] depth = kinect.getRawDepth();
  int skip = 50;
  //translate(width/750,height/750,-50);
    mc.reset();
 
    // original for loop
    println("entering loop");
    int nBalls = 0;
    for(int x=0; x<w; x+=skip) {
      for(int y=0; y<h; y+=skip) {
        int offset = x+y*w;
        int rawDepth = depth[offset];
 
        if(rawDepth >= minDepth && rawDepth <= maxDepth) {
          PVector v = depthToWorld(x,y,rawDepth);
          Vec3D metaBallPos = new Vec3D(v.x * 500, v.y * 300, v.z*300);
          mc.addMetaBall(metaBallPos, 100, 1);
          nBalls++;
        }
 
      }
    }
    println("done with loop, " + nBalls + " balls");
    // end original for loop
 
  mc.createMesh();
  if(bUseFill){
    fill(0,255,0);
    noStroke();
  }
  else {
    noFill();
    stroke(127);
  }
 
  pushMatrix();
  translate(width/2, height/2, 0);
  rotateX(rotationAxis.x);
  rotateY(rotationAxis.y);
  mc.renderMesh();
  popMatrix();
}
 
PVector depthToWorld(int x, int y, int depthValue) {
 
  final double fx_d = 1.0 / 5.9421434211923247e+02;
  final double fy_d = 1.0 / 5.9104053696870778e+02;
  final double cx_d = 3.3930780975300314e+02;
  final double cy_d = 2.4273913761751615e+02;
 
  PVector result = new PVector();
  double depth =  depthLookUp[depthValue];//rawDepthToMeters(depthValue);
  result.x = (float)((x - cx_d) * depth * fx_d);
  result.y = (float)((y - cy_d) * depth * fy_d);
  result.z = (float)(depth);
  return result;
}
 
float rawDepthToMeters(int depthValue) {
  if (depthValue < 2047) {
    return (float)(1.0 / ((double)(depthValue) * -0.0030711016 + 3.3309495161));
  }
  return 0.0f;
}
 
void keyPressed(){
  if(key == CODED){
    if(keyCode == LEFT) rotationAxis.y += 0.05;
    if(keyCode == RIGHT) rotationAxis.y -= 0.05;
    if(keyCode == UP) rotationAxis.x -= 0.05;
    if(keyCode == DOWN) rotationAxis.x += 0.05;
  }
  else {
    if(key == ' '){
      bUseFill = !bUseFill;
    }
    if(key == 'r' || key == 'R'){
      mc.reset();
      rotationAxis.set(0,0,0);
    }
  }
}
 
void stop() {
  kinect.quit();
  super.stop();
}
