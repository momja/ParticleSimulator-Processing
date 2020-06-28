
PShader unlitShader;

int sizeX = 640;
int sizeY = 480;
int slider = 0;
Vec3 camLocation = new Vec3(0, 0, 0);
Vec3 camLookAt = new Vec3(0, 0, 0);
Vec3 camUp = new Vec3(0, -1, 0);
Vec3 mouse3DPoint = null;

float cameraRadius = 30; float phi = 90; float theta = 0;
PlanarParticleSystem particles = new PlanarParticleSystem(15000);
ParticleSystem smoke = new ParticleSystem(1000);
BoidSystem boids = new BoidSystem(700);
PShape[] rigidBodies = new PShape[1];
boolean particleSceneActive = false;
boolean flockToLand = false;
HeadsUpDisplay hud = new HeadsUpDisplay();

PShape house;
PShape scarecrow;
PShape collisionMesh;

void setup() {
  size(640,480,OPENGL);
  perspective(radians(60), 1+1.f/3, 1, 1000);
  surface.setTitle("Particle Simulation [Max Omdal]");
  camLocation = new Vec3(mouseX, height/2, (height/2) / tan(PI/6));
  camLookAt = new Vec3(0, 0, 0);
  camUp = new Vec3(0, -1, 0);
  house = loadShape("house.obj");
  scarecrow = loadShape("scarecrow.obj");
  collisionMesh = loadShape("houseColliderMesh.obj");
  rigidBodies[0] = collisionMesh;
  unlitShader = loadShader("unlit_frag.glsl", "unlit_vert.glsl");

  particles.emitterPosition = new Vec3(0,50,0);
  particles.emitterRadius = 8;
  particles.birthRate = 1000;
  particles.collisionTrigger = new SpawnEmitter();
  particles.streakLength = 1;

  smoke.emitterPosition = new Vec3(1.3,3,2.5);
  smoke.birthRate = 60;
  smoke.r = 0.15;
  smoke.particleSpeed = 0.3;
  smoke.speedRange = 0.1;
  smoke.particleLifespanMax = 10;
  smoke.particleLifespanMin = 9.5;
  smoke.particleDirection = new Vec3(0,1,0);
  smoke.particleDirectionRange = 0.13;
  smoke.particleAcceleration = new Vec3(0.01,0,0.005);
  smoke.particleTexture = loadImage("smokeTex.png");

  boids.boidModel = loadShape("bird.obj");
  boids.boidPerchedModel = loadShape("birdSitting.obj");
}

Vec3 getMouseCast() {
  Vec3 w = camLookAt.minus(camLocation).normalized();
  Vec3 u = cross(w, camUp).normalized();
  Vec3 v = cross(u, w).normalized();

  w.mul(-1);

  float m3dx = map(mouseX, 0, width, -0.666, 0.666);
  float m3dy = map(mouseY, 0, height, -0.5, 0.5);
  float m3dz = -1;

  float m3dx_world = m3dx*u.x + m3dy*v.x + m3dz*w.x + camLocation.x;
  float m3dy_world = m3dx*u.y + m3dy*v.y + m3dz*w.y + camLocation.y;
  float m3dz_world = m3dx*u.z + m3dy*v.z + m3dz*w.z + camLocation.z;

  Vec3 m_world = new Vec3(m3dx_world, m3dy_world, m3dz_world);
  Vec3 rayDir = m_world.minus(camLocation);
  return rayDir;
}

void updateParticles(float dt) {
  particles.generateNewParticles(dt);
  particles.updateParticlePositions(dt);
  particles.checkParticlesForCollisions(rigidBodies);
  particles.updateParticleProperties(dt);
  particles.removeDeadParticles();
  particles.drawTriggers(dt);
  particles.drawAllParticles();
}

void updateBoids(float dt) {
  if (mousePressed) {
    Vec3 raycast = getMouseCast();
    Vec3 p = testForCollisions(camLocation, raycast.normalized(), rigidBodies);
    if (p != null) {
      boids.tetherPoint = p;
      boids.influenceToTetherPoint = 2;
      mouse3DPoint = p;
    }
    if (mouse3DPoint != null) {
      push();
      translate(mouse3DPoint.x, mouse3DPoint.y, mouse3DPoint.z);
      shape(scarecrow);
      pop();
    }
  } else {
    boids.tetherPoint = new Vec3(8,5,-4);
    boids.influenceToTetherPoint = 0.03;
    mouse3DPoint = null;
  }
  boids.updateBoidPositions(dt);
  boids.checkForCollisions(rigidBodies);
  boids.drawBoids();
}

void updateSmoke(float dt) {
  smoke.generateNewParticles(dt);
  smoke.updateParticlePositions(dt);
  smoke.updateParticleProperties(dt);
  smoke.removeDeadParticles();
  smoke.drawTriggers(dt);
  shader(unlitShader);
  smoke.drawAllParticles();
  resetShader();
}

void updateCamera() {
  camLocation.x = cos(radians(theta))*cameraRadius;
  camLocation.y = float(slider)/5;
  camLocation.z = sin(radians(theta))*cameraRadius;
}

void keyPressed() {
    if (key == ' ') {
      // Switch Scenes
      particleSceneActive = !particleSceneActive;
    }
}

void draw() {
  if (keyPressed) {
    if (key == 'w') {
      particles.emitterPosition.x -= 2;
    } else if (key == 's') {
      particles.emitterPosition.x += 2;
    } else if (key == 'a') {
      particles.emitterPosition.z -= 2;
    } else if (key == 'd') {
      particles.emitterPosition.z += 2;
    } else if (keyCode == UP) {
      slider++;
    }  if (keyCode == DOWN) {
      slider--;
    }  if (keyCode == LEFT) {
      theta -= 1.1;
    }  if (keyCode == RIGHT) {
      theta += 1.1;
    } if (key == 'm') {
      boids.addBoid();
    } else if (key == 'n') {
      boids.loseBoid();
    }
  }

  if (flockToLand && boids.tetherPoint.y > -3) {
    boids.tetherPoint.y -= 0.02;
  }
  
  updateCamera();

  pushMatrix();
  camera(camLocation.x, camLocation.y, camLocation.z,
         camLookAt.x,   camLookAt.y,   camLookAt.z,
         camUp.x,       camUp.y,       camUp.z);
  directionalLight(150, 170, 180, 1, -1, -1);
  directionalLight(255, 255, 255, 1, -1, 1);
  ambientLight(50,50,50,0,0,0);
  if (particleSceneActive) {
    // Display particle scene
    background(10,10,15);
    shape(house,0,0);
    updateParticles(1.0/frameRate);
    updateSmoke(1.0/frameRate);
  } else {
    // Display boids
    background(200,170,50);
    shape(house,0,0);
    updateBoids(1.0/frameRate);
  }

  popMatrix();

  hud.draw();
}
