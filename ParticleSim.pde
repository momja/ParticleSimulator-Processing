// Shaders Found on the processing website: https://processing.org/tutorials/pshader/


PShader toon;

Vec3 camLocation = new Vec3(0, 0, 0);
Vec3 camLookAt = new Vec3(0, 0, 0);
Vec3 camUp = new Vec3(0, -1, 0);

float cameraRadius = 100; float phi = 90; float theta = 0;
PlanarParticleSystem particles = new PlanarParticleSystem(15000);
PShape[] rigidBodies = new PShape[1];

PShape house;

void setup() {
  size(640,480,P3D);
  frustum(-1.333,1.333,-1,1,1, 1000);
  surface.setTitle("Wind Simulation [Max Omdal]");
  camLocation = new Vec3(mouseX, height/2, (height/2) / tan(PI/6));
  camLookAt = new Vec3(0, 0, 0);
  camUp = new Vec3(0, -1, 0);

  toon = loadShader("ToonFrag.glsl", "ToonVert.glsl");
  toon.set("fraction", 1.0);
  house = loadShape("house.obj");
  rigidBodies[0] = house;

  particles.emitterPosition = new Vec3(0,50,0);
  particles.emitterRadius = 50;
  particles.collisionTrigger = new SpawnEmitter();
}

void update(float dt) {
  particles.generateNewParticles(dt);
  particles.updateParticlePositions(dt);
  particles.checkParticlesForCollisions(rigidBodies);
  particles.updateParticleProperties(dt);
  particles.removeDeadParticles();
  particles.drawTriggers(dt);
  while (particles.drawNextParticle()) {
  }
}

void updateCamera() {
  cameraRadius = mouseX/2 + 50;
  camLocation.x = cos(radians(theta))*cameraRadius;
  camLocation.y = mouseY/3;
  camLocation.z = sin(radians(theta))*cameraRadius;
  if (!mousePressed) {
    theta += 0.5;
  }
}

void draw() {
  if (keyPressed) {
    if (key == 'w') {
      particles.emitterPosition.x += 2;
    } else if (key == 's') {
      particles.emitterPosition.x -= 2;
    } else if (key == 'a') {
      particles.emitterPosition.z -= 2;
    } else if (key == 'd') {
      particles.emitterPosition.z += 2;
    }
  }
  
  updateCamera();

  camera(camLocation.x, camLocation.y, camLocation.z,
         camLookAt.x,   camLookAt.y,   camLookAt.z,
         camUp.x,       camUp.y,       camUp.z);

  shader(toon);
  background(0); // Black Background
  directionalLight(150, 170, 180, 1, -1, -1);
  shape(house,0,0);
  update(1.0/frameRate);
}
