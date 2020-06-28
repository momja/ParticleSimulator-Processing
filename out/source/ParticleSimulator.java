import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class ParticleSimulator extends PApplet {


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

public void setup() {
  
  perspective(radians(60), 1+1.f/3, 1, 1000);
  surface.setTitle("Wind Simulation [Max Omdal]");
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

  smoke.emitterPosition = new Vec3(1.3f,3,2.5f);
  smoke.birthRate = 60;
  smoke.r = 0.15f;
  smoke.particleSpeed = 0.3f;
  smoke.speedRange = 0.1f;
  smoke.particleLifespanMax = 10;
  smoke.particleLifespanMin = 9.5f;
  smoke.particleDirection = new Vec3(0,1,0);
  smoke.particleDirectionRange = 0.13f;
  smoke.particleAcceleration = new Vec3(0.01f,0,0.005f);
  smoke.particleTexture = loadImage("smokeTex.png");

  boids.boidModel = loadShape("bird.obj");
  boids.boidPerchedModel = loadShape("birdSitting.obj");
}

public Vec3 getMouseCast() {
  Vec3 w = camLookAt.minus(camLocation).normalized();
  Vec3 u = cross(w, camUp).normalized();
  Vec3 v = cross(u, w).normalized();

  w.mul(-1);

  float m3dx = map(mouseX, 0, width, -0.666f, 0.666f);
  float m3dy = map(mouseY, 0, height, -0.5f, 0.5f);
  float m3dz = -1;

  float m3dx_world = m3dx*u.x + m3dy*v.x + m3dz*w.x + camLocation.x;
  float m3dy_world = m3dx*u.y + m3dy*v.y + m3dz*w.y + camLocation.y;
  float m3dz_world = m3dx*u.z + m3dy*v.z + m3dz*w.z + camLocation.z;

  Vec3 m_world = new Vec3(m3dx_world, m3dy_world, m3dz_world);
  Vec3 rayDir = m_world.minus(camLocation);
  return rayDir;
}

public void updateParticles(float dt) {
  particles.generateNewParticles(dt);
  particles.updateParticlePositions(dt);
  particles.checkParticlesForCollisions(rigidBodies);
  particles.updateParticleProperties(dt);
  particles.removeDeadParticles();
  particles.drawTriggers(dt);
  particles.drawAllParticles();
}

public void updateBoids(float dt) {
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
    boids.influenceToTetherPoint = 0.03f;
    mouse3DPoint = null;
  }
  boids.updateBoidPositions(dt);
  boids.checkForCollisions(rigidBodies);
  boids.drawBoids();
}

public void updateSmoke(float dt) {
  smoke.generateNewParticles(dt);
  smoke.updateParticlePositions(dt);
  smoke.updateParticleProperties(dt);
  smoke.removeDeadParticles();
  smoke.drawTriggers(dt);
  shader(unlitShader);
  smoke.drawAllParticles();
  resetShader();
}

public void updateCamera() {
  camLocation.x = cos(radians(theta))*cameraRadius;
  camLocation.y = PApplet.parseFloat(slider)/5;
  camLocation.z = sin(radians(theta))*cameraRadius;
}

public void keyPressed() {
    if (key == ' ') {
      // Switch Scenes
      particleSceneActive = !particleSceneActive;
    }
}

public void draw() {
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
      theta -= 1.1f;
    }  if (keyCode == RIGHT) {
      theta += 1.1f;
    } if (key == 'm') {
      boids.addBoid();
    } else if (key == 'n') {
      boids.loseBoid();
    }
  }

  if (flockToLand && boids.tetherPoint.y > -3) {
    boids.tetherPoint.y -= 0.02f;
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
    updateParticles(1.0f/frameRate);
    updateSmoke(1.0f/frameRate);
  } else {
    // Display boids
    background(200,170,50);
    shape(house,0,0);
    updateBoids(1.0f/frameRate);
  }

  popMatrix();

  hud.draw();
}
public class BoidSystem {
    // private ArrayList<Vec3> boidCoords;
    // private ArrayList<Vec3> boidPrevCoords;
    // // private ArrayList<Quaternion> boidOrientations;
    // private ArrayList<Vec3> boidVelocities;
    // private ArrayList<Float> perchedBoidCountdowns;
    private ArrayList<Boid> boids;
    private float width, height, length;

    Vec3 boundingBoxOrigin = new Vec3(0,25,0);
    float minX, maxX, minY, maxY, minZ, maxZ;

    int boidCount = 100;
    float boidSize = 2;
    float maxSpeed = 100;
    float separation = 2;
    float visualDistance = 6;
    float influenceToCenter = 0.02f;
    Vec3 boidColor = new Vec3(0,0,0);
    PImage boidTexture = null;
    PShape boidModel = null;
    PShape boidPerchedModel = null;
    Vec3 tetherPoint = new Vec3(8,5,-4);
    float influenceToTetherPoint = 0.03f;

    public BoidSystem(int boidCount) {
        // boidCoords = new ArrayList<Vec3>();
        // boidPrevCoords = new ArrayList<Vec3>();
        // // boidOrientations = new ArrayList<Quaternion>();
        // boidVelocities = new ArrayList<Vec3>();
        // perchedBoidCountdowns = new ArrayList<Float>();
        boids = new ArrayList<Boid>();
        setBoundingBox(boundingBoxOrigin, 100, 60, 100);
        this.boidCount = boidCount;
        initializePositions();
    }

    public void setBoundingBox(Vec3 origin, float width, float height, float length) {
        this.width = width;
        this.height = height;
        this.length = length;
        this.boundingBoxOrigin = origin;

        minX = origin.x - width/2;
        maxX = origin.x + width/2;
        minY = origin.y - height/2;
        maxY = origin.y + height/2;
        minZ = origin.z - length/2;
        maxZ = origin.z + length/2;
    }

    public void checkForCollisions(PShape[] rigidBodies) {
        // Check each boid to see if it has collided with the collision meshes
        for (PShape rigidBody : rigidBodies) {
        int triCount = rigidBody.getChildCount();
        for (int i = 0; i < triCount; i++) {
            PShape triangle = rigidBody.getChild(i);
            int j = 0;
            while(j < boidCount) {
                Boid boid = boids.get(j);
                Vec3 boidPosition = boid.coords;
                if (boid.perchCountdown > 0 || boidPosition.y > 5 || boidPosition.x > 20 || boidPosition.x < -20 || boidPosition.z > 20 || boidPosition.z < -20) {
                    j++;
                    continue;
                }
                // TODO: Use Barycentric Coordinates to find if there is a collision with the surface
                boolean collision = false;
                Vec3 collisionPoint = new Vec3(0,0,0);

                Vec3 rayOrigin = boidPosition;
                Vec3 rayDirection = boidPosition.minus(boid.previousCoord);
                float maxT = rayDirection.length();

                if (maxT < 0.00001f) {
                    j++;
                    continue;
                }

                rayDirection.normalize();

                PVector v1 = triangle.getVertex(0);
                PVector v2 = triangle.getVertex(1);
                PVector v3 = triangle.getVertex(2);

                Vec3 vert1 = new Vec3(v1.x, v1.y, v1.z);
                Vec3 vert2 = new Vec3(v2.x, v2.y, v2.z);
                Vec3 vert3 = new Vec3(v3.x, v3.y, v3.z);

                Vec3 e1 = vert2.minus(vert1);
                Vec3 e2 = vert3.minus(vert1);

                Vec3 surfaceNormal = cross(e1, e2);
                // float x_0 = rayOrigin.x; float y_0 = rayOrigin.y; float z_0 = rayOrigin.z;
                // float x_d = rayDirection.x; float y_d = rayDirection.y; float z_d = rayDirection.z;
                float denominator = dot(surfaceNormal, rayDirection);
                if (abs(denominator) <= 0.0001f) {
                    // No ray plane intersection exists
                    j++;
                    continue;
                }

                float D = dot(vert1, surfaceNormal);

                float numerator = -(dot(surfaceNormal, rayOrigin) - D);

                float t = numerator/denominator;

                if (t < 0) {
                    // Haven't hit yet
                    j++;
                    continue;
                }
                
                Vec3 p = rayOrigin.plus(rayDirection.times(t));

                if (t < maxT && pointLiesOnTriangle(p, vert1, vert2, vert3, e1, e2)) {
                    perchBoid(j, p, reflect(rayDirection, surfaceNormal));
                }
                j++;
                }
            }
        }
    }

    public void perchBoid(int idx, Vec3 p, Vec3 launchDir) {
        // Hold boid at current position for a ranged amount of time.
        // Then release the boid in the direction provided
        Boid boid = boids.get(idx);
        boid.coords = p;
        boid.perchCountdown = random(3,4);
        boid.velocity = launchDir.normalized().times(maxSpeed/3);
    }
    
    public void addBoid() {
        Vec3 position = new Vec3(random(minX, maxX),
                                     random(minY, maxY)/2,
                                     random(minZ, maxZ));
        Boid boid = new Boid();
        boid.coords = position;
        boid.previousCoord = position;
        boid.perchCountdown = 0.f;

        Vec3 velocity = new Vec3(random(-1,1),random(-1,1),random(-1,1));
        velocity.normalize();
        velocity.mul(random(maxSpeed/2, maxSpeed));
        boid.velocity = velocity;

        boids.add(boid);
        boidCount++;
    }

    public void loseBoid() {
        boids.remove(0);
        boidCount--;
    }

    public void initializePositions() {
        // Given the bounding box, we will randomly place our boids within this box
        for (int i = 0; i < boidCount; i++) {
            Vec3 position = new Vec3(random(minX, maxX),
                                     random(minY, maxY)/2,
                                     random(minZ, maxZ));
            Boid boid = new Boid();
            boid.coords = position;
            boid.previousCoord = position;
            boid.perchCountdown = 0.f;

            Vec3 velocity = new Vec3(random(-1,1),random(-1,1),random(-1,1));
            velocity.normalize();
            velocity.mul(random(maxSpeed/2, maxSpeed));
            boid.velocity = velocity;

            boids.add(boid);
        }
    }

    public void updateBoidPositions(float dt) {
        Vec3 v1, v2, v3, v4;
        for (Boid boid : boids) {
            if(boid.perchCountdown > 0) {
                boid.perchCountdown -= dt;
                continue;
            }

            v1 = flyTowardsCenter(boid);
            v2 = keepDistance(boid);
            v3 = matchVelocity(boid);
            v4 = moveToTetherPoint(boid);

            // Velocity
            boid.velocity.add(v1.times(0.6f).plus(v2.times(0.6f)).plus(v3).plus(v4));
            boid.velocity.y *= 0.9f;
            boid.velocity.clamp(0, maxSpeed);
            // Position
            boid.previousCoord = new Vec3(boid.coords);
            boid.coords.add(boid.velocity.times(dt));

            boundPosition(boid);
        }
    }

    public void drawBoids() {
        if (boidTexture != null) {
            for (Boid boid : boids) {
                // Draw texture
                Vec3 pos = boid.coords;
                push();
                translate(pos.x, pos.y, pos.z);
                beginShape();
                texture(boidTexture);
                vertex(-boidSize/2,boidSize/2,0,0,0);
                vertex(boidSize/2,boidSize/2,0,boidTexture.width,0);
                vertex(boidSize/2,-boidSize/2,0,boidTexture.width,boidTexture.height);
                vertex(-boidSize/2,-boidSize/2,0,0,boidTexture.height);
                endShape();
                pop();
            }
        }
        else if (boidModel != null) {
            for (Boid boid : boids) {
                // Draw boid model
                Vec3 pos = boid.coords;
                if (pos.distanceTo(camLocation) > 30) {
                    push();
                    stroke(boidColor.x, boidColor.y, boidColor.z);
                    strokeWeight(3);
                    point(pos.x, pos.y, pos.z);
                    pop();
                    continue;
                }
                Vec3 vel = boid.velocity.normalized();
                pushMatrix();
                translate(pos.x, pos.y, pos.z);
                Vec3 xz_vec = new Vec3(vel.x, 0, vel.z);
                xz_vec.normalize();
                float angle = acos(dot(xz_vec, new Vec3(0,0,1)));
                if (cross(xz_vec, new Vec3(0,0,1)).y > 0) {
                    rotateY(-angle-PI/2);
                } else {
                    rotateY(angle-PI/2);
                }
                if (boid.perchCountdown > 0) {
                    shape(boidPerchedModel);
                } else {
                    shape(boidModel);
                }
                popMatrix();
            }
        }
        else {
            for (Boid boid : boids) {
                // Draw a point
                Vec3 pos = boid.coords;
                push();
                stroke(255);
                strokeWeight(boidSize);
                point(pos.x, pos.y, pos.z);
                pop();
            }
        }
    }

    private void boundPosition(Boid boid) {
        Vec3 pos = boid.coords;
        Vec3 vel = boid.velocity;
        // Check for boids outside the bounding box
        if (pos.x > maxX) {
            pos.x = maxX;
            vel.x = -20;
        } else if (pos.x < minX) {
            pos.x = minX;
            vel.x = 20;
        }
        if (pos.y > maxY) {
            pos.y = maxY;
            vel.y = -20;
        } else if (pos.y < minY) {
            pos.y = minY;
            vel.y = 20;
        }
        if (pos.z > maxZ) {
            pos.z = maxZ;
            vel.z = -20;
        } else if (pos.z < minZ) {
            pos.z = minZ;
            vel.z = 20;
        }
    }

    private Vec3 flyTowardsCenter(Boid boid) {
        // given an index of a certain boid,
        // find the center of mass of nearby boids
        // and return a weight in that direction
        
        // Fly towards the center of mass excluding the current boid
        Vec3 perceivedCenterOfMass = new Vec3(0,0,0);
        int cnt = 0;
        int i = 0;
        for (Boid otherBoid : boids) {
            if (otherBoid.equals(boid) || otherBoid.perchCountdown > 0) {
                i++;
                continue;
            } else if (boid.coords.distanceTo(otherBoid.coords) < visualDistance) {
                perceivedCenterOfMass.add(otherBoid.coords);
                cnt++;
            }
            i++;
        }
        if (cnt > 0) {
            perceivedCenterOfMass.mul(1.f/(cnt));
        }
        return perceivedCenterOfMass.minus(boid.coords).times(influenceToCenter);
    }

    private Vec3 keepDistance(Boid boid) {
        // look at all the nearby boids and objects
        // and return a weight that maneuvers away from them
        Vec3 c = new Vec3(0,0,0);

        int i = 0;
        for (Boid otherBoid : boids) {
            if (otherBoid.equals(boid)) {
                i++;
                continue;
            }
            Vec3 betweenVec = otherBoid.coords.minus(boid.coords);
            if (betweenVec.length() < separation) {
                c.subtract(betweenVec);
            }
            i++;
        }
        return c;
    }

    private Vec3 matchVelocity(Boid boid) {
        // find the velocities of nearby boids and approximate it
        Vec3 curBoidVelocty = boid.velocity;
        Vec3 newPosition = new Vec3(0,0,0);
        int cnt = 0;
        int i = 0;
        for (Boid otherBoid : boids) {
            if (otherBoid.equals(boid) || otherBoid.perchCountdown > 0) {
                i++;
                continue;
            } else if (boid.coords.distanceTo(otherBoid.coords) < visualDistance) {
                newPosition.add(otherBoid.velocity);
                cnt++;
            }
            i++;
        }
        if (cnt > 0) {
            newPosition.mul(1.f/(cnt));
        }
        newPosition.subtract(curBoidVelocty);
        newPosition.mul(1.f/8);
        return newPosition;
    }

    public Vec3 moveToTetherPoint(Boid boid) {
        return tetherPoint.minus(boid.coords).normalized().times(influenceToTetherPoint);
    }

    public Vec3 noisyMovement() {
        return new Vec3(random(-maxSpeed/200, maxSpeed/200),
                        random(-maxSpeed/200, maxSpeed/200),
                        random(-maxSpeed/200, maxSpeed/200));
    }
}

public class Boid {
    Vec3 coords;
    Vec3 previousCoord;
    float perchCountdown;
    Vec3 velocity;
}
public class CollisionTrigger {
    boolean isActive = false;

    public void onCollision(Vec3 point, Vec3 normal) {
        isActive = true;
        return;
    }

    public CollisionTrigger copy() {
        return new CollisionTrigger();
    }

    public void draw(float dt) {
    }
}

public class SpawnEmitter extends CollisionTrigger {
    ParticleSystem emitter;

    public SpawnEmitter() {
        emitter = new ParticleSystem(100);
        emitter.streakLength = 0;
        emitter.particleLifespanMax = 0.3f;
        emitter.particleLifespanMin = 0.2f;
        emitter.birthRate = 10;
        emitter.emitterLifespan = 0.1f;
        emitter.r = 1.3f;
        emitter.particleSpeed = 20;
        emitter.speedRange = 10;
        emitter.particleDirection = new Vec3(0,1,0);
        emitter.particleDirectionRange = 0.1f;
        emitter.particleAcceleration = new Vec3(0,-500,0);
        emitter.particleAccelerationRange = 0.1f;
    }

    @Override
    public SpawnEmitter copy() {
        return new SpawnEmitter();
    }

    @Override
    public void onCollision(Vec3 point, Vec3 normal) {
        // TODO : Spawn a new emitter at the point of collision
        super.onCollision(point, normal);
        emitter.emitterPosition = point;
        emitter.particleDirection = normal;
    }

    public void draw(float dt) {
        emitter.generateNewParticles(dt);
        emitter.updateParticlePositions(dt);
        emitter.updateParticleProperties(dt);
        emitter.removeDeadParticles();
        while (emitter.drawNextParticle()) {
        }
        emitter.drawTriggers(dt);
        this.isActive = emitter.isActive;
    }
}

public class AnimateRaindropCollision extends CollisionTrigger {
    @Override
    public void onCollision(Vec3 point, Vec3 normal) {
        // TODO : Spawn a textured quad that animates through a series
        // of raindrop splash images
        super.onCollision(point, normal);
    }
}

public class TriggerCollection {
    public ArrayList<CollisionTrigger> triggers;

    public TriggerCollection() {
        triggers = new ArrayList<CollisionTrigger>();
    }

    public void add(CollisionTrigger trigger) {
        triggers.add(trigger);
    }

    public void drawAllTriggers(float dt) {
        int i = 0;
        while (i < triggers.size()) {
            CollisionTrigger trigger = triggers.get(i);

            trigger.draw(dt);

            if (!trigger.isActive) {
                // Remove any inactive triggers
                triggers.remove(i);
                i--;
            }
            i++;
        }
    }
}
class HeadsUpDisplay {
    public boolean showFrameRate = true;
    public PFont font;
    public float fontSize;
    public int defaultColor = color(0, 0, 0);

    public HeadsUpDisplay() {
    }

    public void draw() {
        camera();
        hint(DISABLE_DEPTH_TEST);
        push();
        noLights();
        textMode(MODEL);
        textSize(24);
        textAlign(LEFT, TOP);
        text("fps "+round(frameRate), 10, 10);
        pop();
        hint(ENABLE_DEPTH_TEST);
    }

    public void setFont(PFont font) {
        this.font = font;
    }

    public void setFontSize(float fontSize) {
        this.fontSize = fontSize;
    }

    public void setColor(int c) {
        this.defaultColor = c;
    }
}

public class ParticleSystem {
  protected ArrayList<Vec3> particleCoords;
  protected ArrayList<Vec3> previousParticleCoords;
  protected ArrayList<Vec3> particleVelocities;
  protected ArrayList<Vec3> particleAccelerations;
  protected ArrayList<Float> particleLifespan;
  protected ArrayList<Vec3> particleColors;
  protected ArrayList<Float> particleRadii;
  protected ArrayList<Float> particleStreakLength;

  int partIdx = 0;
  int streakLength = 0;
  int particleCount = 0;
  int maxParticleCount;
  float particleLifespanMax = 5;
  float particleLifespanMin = 5;
  float emitterLifespan = -1;
  float emitterElapsedTime = 0;
  boolean isActive = true;
  float birthRate = 100;
  float r = 0.2f;
  float particleSpeed = 200;
  float speedRange = 50;
  Vec3 particleDirection = new Vec3(0,-1,0);
  float particleDirectionRange = 0.07f;
  Vec3 particleAcceleration = new Vec3(0,0,0);
  float particleAccelerationRange = 0;
  Vec3 emitterPosition = new Vec3(0,0,0);
  Vec3 particleColor = new Vec3(230,245,255);
  Vec3 particleColorRange = new Vec3(0,0,0);
  PImage particleTexture = null;
  CollisionTrigger collisionTrigger = null;
  TriggerCollection triggerCollection = new TriggerCollection();
  
  public ParticleSystem(int maxParticleCount) {
    particleCoords = new ArrayList<Vec3>();
    previousParticleCoords = new ArrayList<Vec3>();
    particleVelocities = new ArrayList<Vec3>();
    particleAccelerations = new ArrayList<Vec3>();
    particleLifespan = new ArrayList<Float>();
    particleColors = new ArrayList<Vec3>();
    particleRadii = new ArrayList<Float>();
    particleStreakLength = new ArrayList<Float>();
    this.maxParticleCount = maxParticleCount;
  }

  public void generateNewParticles(float dt) {
    // given the change in time, (assuming since the last particle generation)
    // generate more particles and add to the ArrayList of existing particles.
    // Initialize the new particles time, velocities, accelerations, lifespan etc.
    if (emitterLifespan > 0 && emitterElapsedTime >= emitterLifespan) {
      birthRate = 0;
      return;
    }
    emitterElapsedTime += dt;

    float newParticlesToGen = dt * birthRate;
    float decimal = newParticlesToGen - (float)(floor(newParticlesToGen));
    int stochasticNewParticles = floor(newParticlesToGen);
    if (decimal*100 > random(100)) {
      stochasticNewParticles++;
    }
    if (maxParticleCount < particleCount + stochasticNewParticles) {
      stochasticNewParticles = maxParticleCount - particleCount;
    }
    for(int i = 0; i < stochasticNewParticles; i++) {
      particleCoords.add(new Vec3(emitterPosition));
      previousParticleCoords.add(new Vec3(emitterPosition));
      Vec3 particleDir = new Vec3(particleDirection);
      Vec3 randomParticleDir = new Vec3(random(-particleDirectionRange,particleDirectionRange),
                                  random(-particleDirectionRange,particleDirectionRange),
                                  random(-particleDirectionRange,particleDirectionRange));
      particleDir.add(randomParticleDir);
      particleDir.normalize();
      particleVelocities.add(particleDir.times(particleSpeed + random(-speedRange, speedRange)));
      Vec3 particleAccel = new Vec3(particleAcceleration);
      Vec3 randomParticleAccel = new Vec3(random(-particleAccelerationRange,particleAccelerationRange),
                                          random(-particleAccelerationRange,particleAccelerationRange),
                                          random(-particleAccelerationRange,particleAccelerationRange));
      particleAccel.add(randomParticleAccel);
      particleAccelerations.add(particleAcceleration);
      float randLifespan = random(particleLifespanMin, particleLifespanMax);
      particleLifespan.add(randLifespan);
      particleColors.add(new Vec3(particleColor));
      particleRadii.add(r);
      float startingStreakLength = 0.0f;
      particleStreakLength.add(startingStreakLength);
      particleCount++;
    }
  }
  
  public void updateParticlePositions(float dt) {
    // using the particle velocities, acceleration and current position,
    // find its new position.
    for(int i = 0; i < particleCount; i++) {
      previousParticleCoords.set(i, new Vec3(particleCoords.get(i)));
      Vec3 vel = particleVelocities.get(i);
      Vec3 accel = particleAccelerations.get(i);
      Vec3 translation = vel.times(dt).plus(accel.times(0.5f).times(dt*dt));
      particleCoords.get(i).add(translation);
      vel.add(accel.times(dt));
      Float partStreak = particleStreakLength.get(i);
      if (partStreak < streakLength) {
        partStreak += translation.length();
      }
    }
  }
  
  public void checkParticlesForCollisions(PShape[] rigidBodies) {
    // see if any particles are intersecting anything. Move to outside that object
    // and recalculate velocity and acceleration

    if (collisionTrigger == null) {
      // If the collision trigger is null, we will ignore all rigidbodies
      return;
    }

    for (PShape rigidBody : rigidBodies) {
      int triCount = rigidBody.getChildCount();
      for (int i = 0; i < triCount; i++) {
        PShape triangle = rigidBody.getChild(i);
        int j = 0;
        while(j < particleCount) {
          Vec3 particlePosition = particleCoords.get(j);
          // TODO: Use Barycentric Coordinates to find if there is a collision with the surface
          boolean collision = false;
          Vec3 collisionPoint = new Vec3(0,0,0);

          Vec3 rayOrigin = particlePosition;
          Vec3 rayDirection = previousParticleCoords.get(j).minus(particlePosition);
          float maxT = rayDirection.length();

          if (maxT < 0.00001f) {
            j++;
            continue;
          }

          rayDirection.normalize();

          PVector v1 = triangle.getVertex(0);
          PVector v2 = triangle.getVertex(1);
          PVector v3 = triangle.getVertex(2);

          Vec3 vert1 = new Vec3(v1.x, v1.y, v1.z);
          Vec3 vert2 = new Vec3(v2.x, v2.y, v2.z);
          Vec3 vert3 = new Vec3(v3.x, v3.y, v3.z);

          Vec3 e1 = vert2.minus(vert1);
          Vec3 e2 = vert3.minus(vert1);

          Vec3 surfaceNormal = cross(e1, e2);
          // float x_0 = rayOrigin.x; float y_0 = rayOrigin.y; float z_0 = rayOrigin.z;
          // float x_d = rayDirection.x; float y_d = rayDirection.y; float z_d = rayDirection.z;
          float denominator = dot(surfaceNormal, rayDirection);
          if (abs(denominator) <= 0.0001f) {
            // No ray plane intersection exists
            j++;
            continue;
          }

          float D = dot(vert1, surfaceNormal);

          float numerator = -(dot(surfaceNormal, rayOrigin) - D);

          float t = numerator/denominator;

          if (t < 0) {
            // Haven't hit yet
            j++;
            continue;
          }
          
          Vec3 p = rayOrigin.plus(rayDirection.times(t));

          if (t < maxT && pointLiesOnTriangle(p, vert1, vert2, vert3, e1, e2)) {
            CollisionTrigger newTrig = collisionTrigger.copy();
            newTrig.onCollision(p, surfaceNormal.normalized());
            triggerCollection.add(newTrig);
            // remove particle
            removeParticleAtIndex(j);
            j--;
          }
          j++;
        }
      }
    }
  }
  
  public void updateParticleProperties(float dt) {
    // update the color, lifespan, etc of the particles.
    for(int i = 0; i < particleCount; i++) {
      particleLifespan.set(i, particleLifespan.get(i) - dt);
    }
  }
  
  public void removeDeadParticles() {
    // any particles that are older than the particle lifespan
    // should be popped off the list
    int i = 0;
    while(i < particleCount) {
      if(particleLifespan.get(i) < 0.0f) {
        removeParticleAtIndex(i);
        i--;
      }
      i++;
    }
    if (particleCount == 0 && emitterElapsedTime >= emitterLifespan) {
      isActive = true;
    }
  }

  public void removeParticleAtIndex(int i) {
    particleCoords.remove(i);
    previousParticleCoords.remove(i);
    particleVelocities.remove(i);
    particleAccelerations.remove(i);
    particleLifespan.remove(i);
    particleColors.remove(i);
    particleRadii.remove(i);
    particleStreakLength.remove(i);
    particleCount--;
  }

  
  public ArrayList<Vec3> getParticleCoords() {
    return particleCoords;
  }

  public ArrayList<Vec3> getParticleColors() {
    return particleColors;
  }
  
  public ArrayList<Float> getParticleRadii() {
    return particleRadii;
  }

  public void drawAllParticles() {
    hint(ENABLE_DEPTH_SORT);
    while (drawNextParticle()) {

    }
    hint(DISABLE_DEPTH_SORT);
  }

  public boolean drawNextParticle() {
    if (partIdx < particleCount) {
      push();
      stroke(particleColors.get(partIdx).x, particleColors.get(partIdx).y, particleColors.get(partIdx).z);
      strokeWeight(particleRadii.get(partIdx));
      Vec3 vel = new Vec3(particleVelocities.get(partIdx));
      float velMagnitude = vel.length();
      vel.normalize();
      line(particleCoords.get(partIdx).x,
         particleCoords.get(partIdx).y,
         particleCoords.get(partIdx).z,
         particleCoords.get(partIdx).x - streakLength*vel.x*velMagnitude/70,
         particleCoords.get(partIdx).y - streakLength*vel.y*velMagnitude/70,
         particleCoords.get(partIdx).z - streakLength*vel.z*velMagnitude/70);
      pop();
      if (particleTexture == null) {
        stroke(particleColors.get(partIdx).x, particleColors.get(partIdx).y, particleColors.get(partIdx).z);
        strokeWeight(particleRadii.get(partIdx));
        point(particleCoords.get(partIdx).x, particleCoords.get(partIdx).y, particleCoords.get(partIdx).z);
      } else {
        push();
        noStroke();
        translate(particleCoords.get(partIdx).x, particleCoords.get(partIdx).y, particleCoords.get(partIdx).z);
        rotateY (-radians(theta+270));
        beginShape();
        texture(particleTexture);
        float width_2 = particleRadii.get(partIdx);
        vertex(-width_2,width_2,0,0,0);
        vertex(width_2,width_2,0,particleTexture.width,0);
        vertex(width_2,-width_2,0,particleTexture.width,particleTexture.height);
        vertex(-width_2,-width_2,0,0,particleTexture.height);
        endShape();
        pop();
      }
      partIdx++;
      return true;
    }
    else {
      partIdx = 0;
      return false;
    }
  }

  public void drawTriggers(float dt) {
    triggerCollection.drawAllTriggers(dt);
  }
}


public class PlanarParticleSystem extends ParticleSystem {
  Vec3 emitterPlaneNormal = new Vec3(0,1,0);
  float emitterRadius = 10;

  public PlanarParticleSystem(int maxParticleCount) {
    super(maxParticleCount);
  }

  @Override
  public void generateNewParticles(float dt) {
    // given the change in time, (assuming since the last particle generation)
    // generate more particles and add to the ArrayList of existing particles.
    // Initialize the new particles time, velocities, accelerations, lifespan etc.
    if (emitterLifespan > 0 && emitterElapsedTime >= emitterLifespan) {
      birthRate = 0;
      isActive = false;
      return;
    }
    emitterElapsedTime += dt;
    float newParticlesToGen = dt * birthRate;
    float decimal = newParticlesToGen - (float)(floor(newParticlesToGen));
    int stochasticNewParticles = floor(newParticlesToGen);
    if (decimal*100 > random(100)) {
      stochasticNewParticles++;
    }
    if (maxParticleCount < particleCount + stochasticNewParticles) {
      stochasticNewParticles = maxParticleCount - particleCount;
    }
    for(int i = 0; i < stochasticNewParticles; i++) {
      float R = emitterRadius*sqrt(random(0,1));
      float theta = random(0,1) * 2 * PI;
      Vec3 radialPosition = new Vec3(R*cos(theta), 0, R*sin(theta));
      float d = dot(emitterPosition.times(-1), emitterPlaneNormal);
      float k = -(dot(emitterPlaneNormal, radialPosition) + d);
      Vec3 projPointOntoPlane = radialPosition.plus(emitterPlaneNormal.times(k));
      previousParticleCoords.add(emitterPosition.plus(projPointOntoPlane));
      particleCoords.add(emitterPosition.plus(projPointOntoPlane));
      Vec3 particleDir = new Vec3(particleDirection);
      Vec3 randomParticleDir = new Vec3(random(-particleDirectionRange,particleDirectionRange),
                                  random(-particleDirectionRange,particleDirectionRange),
                                  random(-particleDirectionRange,particleDirectionRange));
      particleDir.add(randomParticleDir);
      particleDir.normalize();
      particleVelocities.add(particleDir.times(particleSpeed + random(-speedRange, speedRange)));
      Vec3 particleAccel = new Vec3(particleAcceleration);
      Vec3 randomParticleAccel = new Vec3(random(-particleAccelerationRange,particleAccelerationRange),
                                          random(-particleAccelerationRange,particleAccelerationRange),
                                          random(-particleAccelerationRange,particleAccelerationRange));
      particleAccel.add(randomParticleAccel);
      particleAccelerations.add(particleAcceleration);
      float randLifespan = random(particleLifespanMin, particleLifespanMax);
      particleLifespan.add(randLifespan);
      particleColors.add(new Vec3(particleColor));
      particleRadii.add(r);
      float startingStreakLength = 0.0f;
      particleStreakLength.add(startingStreakLength);
      particleCount++;
    }
  }
}
// Inspiration from https://www.cprogramming.com/tutorial/3d/quaternions.html

public class Quaternion {
    float x, y, z, w;
    public Quaternion(float x, float y, float z, float w) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
        this.normalize();
    }

    public void normalize() {
        float magnitude = sqrt(w*w + x*x + y*y + z*z);
        w /= magnitude;
        x /= magnitude;
        y /= magnitude;
        z /= magnitude;
    }

    public Quaternion times(Quaternion rhs) {
        return new Quaternion(w*rhs.x + x*rhs.w + y*rhs.z - z*rhs.y,
                              w*rhs.y - x*rhs.z + y*rhs.w + z*rhs.x,
                              w*rhs.z - x*rhs.y - y*rhs.x + z*rhs.w,
                              w*rhs.w - x*rhs.x - y*rhs.y - z*rhs.z);
    }

    public Vec3 toEulerAngles() {
        // https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
        Vec3 angles = new Vec3(0,0,0);

        // roll (x-axis rotation)
        float sinr_cosp = 2 * (w * x + y * z);
        float cosr_cosp = 1 - 2 * (x * x + y * y);
        angles.x = atan2(sinr_cosp, cosr_cosp);

        // pitch (y-axis rotation)
        float sinp = 2 * (w * y - z * x);
        if (abs(sinp) >= 1)
            if (sinp < 0) {
                angles.y = -PI/2;
            } else {
                angles.y = PI/2;
            }
        else
            angles.y = asin(sinp);

        // yaw (z-axis rotation)
        float siny_cosp = 2 * (w * z + x * y);
        float cosy_cosp = 1 - 2 * (y * y + z * z);
        angles.z = atan2(siny_cosp, cosy_cosp);

        return angles;
    }
}
public Vec3 testForCollisions(Vec3 origin, Vec3 dir, PShape[] rigidBodies) {
 // Check each boid to see if it has collided with the collision meshes
 Vec3 closestPoint = null;
 float tMin = Float.MAX_VALUE;
        for (PShape rigidBody : rigidBodies) {
        int triCount = rigidBody.getChildCount();
        for (int i = 0; i < triCount; i++) {
            PShape triangle = rigidBody.getChild(i);
            int j = 0;
            Vec3 boidPosition = origin;

            Vec3 rayOrigin = boidPosition;
            Vec3 rayDirection = dir;
            rayDirection.normalize();

            PVector v1 = triangle.getVertex(0);
            PVector v2 = triangle.getVertex(1);
            PVector v3 = triangle.getVertex(2);

            Vec3 vert1 = new Vec3(v1.x, v1.y, v1.z);
            Vec3 vert2 = new Vec3(v2.x, v2.y, v2.z);
            Vec3 vert3 = new Vec3(v3.x, v3.y, v3.z);

            Vec3 e1 = vert2.minus(vert1);
            Vec3 e2 = vert3.minus(vert1);

            Vec3 surfaceNormal = cross(e1, e2);
            // float x_0 = rayOrigin.x; float y_0 = rayOrigin.y; float z_0 = rayOrigin.z;
            // float x_d = rayDirection.x; float y_d = rayDirection.y; float z_d = rayDirection.z;
            float denominator = dot(surfaceNormal, rayDirection);
            if (abs(denominator) <= 0.0001f) {
                // No ray plane intersection exists
                continue;
            }

            float D = dot(vert1, surfaceNormal);

            float numerator = -(dot(surfaceNormal, rayOrigin) - D);

            float t = numerator/denominator;

            if (t < 0) {
                // Haven't hit yet
                continue;
            }
            
            Vec3 p = rayOrigin.plus(rayDirection.times(t));

            if (t < tMin && pointLiesOnTriangle(p, vert1, vert2, vert3, e1, e2)) {
                closestPoint = p;
                t = tMin;
            }
            
        }
        }
    return closestPoint;
}
//Vector Library [2D]
//CSCI 5611 Vector 3 Library [Incomplete]

//Instructions: Add 3D versions of all of the 2D vector functions
//              Vec3 must also support the cross product.
public class Vec2 {
  public float x, y;
  
  public Vec2(float x, float y){
    this.x = x;
    this.y = y;
  }
  
  public String toString(){
    return "(" + x + ", " + y + ")";
  }
  
  public float length(){
    return sqrt(x*x + y*y);
  }
  
  public Vec2 plus(Vec2 rhs){
    return new Vec2(rhs.x+this.x, rhs.y+this.y);
  }
  
  public void add(Vec2 rhs){
    this.x += rhs.x;
    this.y += rhs.y;
  }
  
  public Vec2 minus(Vec2 rhs){
    return new Vec2(this.x-rhs.x, this.y-rhs.y);
  }
  
  public void subtract(Vec2 rhs){
    this.x -= rhs.x;
    this.y -= rhs.y;
  }
  
  public Vec2 times(float rhs){
    return new Vec2(this.x*rhs, this.y*rhs);
  }
  
  public void mul(float rhs){
    this.x *= rhs;
    this.y *= rhs;
  }
  
  public void normalize(){
    float magnitude = this.length();
    this.x /= magnitude;
    this.y /= magnitude;
  }
  
  public Vec2 normalized(){
    float magnitude = this.length();
    return new Vec2(this.x/magnitude, this.y/magnitude);
  }
  
  public float distanceTo(Vec2 rhs){
    return this.minus(rhs).length();
  }
}

public Vec2 interpolate(Vec2 a, Vec2 b, float t){
  return a.plus((b.minus(a)).times(t));
}

public float interpolate(float a, float b, float t) {
   return a + (b - a)*t; 
}

public float dot(Vec2 a, Vec2 b){
  return a.x*b.x + a.y*b.y;
}

public Vec2 projAB(Vec2 a, Vec2 b){
  return b.times(dot(a, b));
}
//Vector Library [2D]
//CSCI 5611 Vector 3 Library [Incomplete]

//Instructions: Add 3D versions of all of the 2D vector functions
//              Vec3 must also support the cross product.
public class Vec3 {
  public float x, y, z;
  
  public Vec3(float x, float y, float z){
    this.x = x;
    this.y = y;
    this.z = z;
  }

  public Vec3(Vec3 copyVec) {
    this.x = copyVec.x;
    this.y = copyVec.y;
    this.z = copyVec.z;
  }
  
  public String toString(){
    return "(" + x + ", " + y + ", " + z + ")";
  }
  
  public float length(){
    return sqrt(x*x + y*y + z*z);
  }
  
  public Vec3 plus(Vec3 rhs){
    return new Vec3(rhs.x+this.x, rhs.y+this.y, rhs.z+this.z);
  }
  
  public void add(Vec3 rhs){
    this.x += rhs.x;
    this.y += rhs.y;
    this.z += rhs.z;
  }
  
  public Vec3 minus(Vec3 rhs){
    return new Vec3(this.x-rhs.x, this.y-rhs.y, this.z-rhs.z);
  }
  
  public void subtract(Vec3 rhs){
    this.x -= rhs.x;
    this.y -= rhs.y;
    this.z -= rhs.z;
  }
  
  public Vec3 times(float rhs){
    return new Vec3(this.x*rhs, this.y*rhs, this.z*rhs);
  }

  public Vec3 times(Vec3 v) {
    return new Vec3(this.x*v.x, this.y*v.y, this.z*v.z);
  }
  
  public void mul(float rhs){
    this.x *= rhs;
    this.y *= rhs;
    this.z *= rhs;
  }
  
  public void normalize(){
    float magnitude = this.length();
    this.x /= magnitude;
    this.y /= magnitude;
    this.z /= magnitude;
  }
  
  public Vec3 normalized(){
    float magnitude = this.length();
    return new Vec3(this.x/magnitude, this.y/magnitude, this.z/magnitude);
  }
  
  public float distanceTo(Vec3 rhs){
    return this.minus(rhs).length();
  }

  public void clamp(float minLength, float maxLength) {
    float curLength = length();
    if (curLength > maxLength) {
      this.normalize();
      this.mul(maxLength);
    } else if (curLength < minLength) {
      this.normalize();
      this.mul(minLength);
    }
  }
}

public Vec3 interpolate(Vec3 a, Vec3 b, float t){
  return a.plus((b.minus(a)).times(t));
}

public float dot(Vec3 a, Vec3 b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

public Vec3 cross(Vec3 a, Vec3 b){
  return new Vec3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

public Vec3 projAB(Vec3 a, Vec3 b){
  return b.times(dot(a, b));
}

public Vec3 reflect(Vec3 d, Vec3 n) {
  Vec3 r = d.minus(n.times(dot(d,n.normalized())*2));
  return r;
}

public boolean pointLiesOnTriangle(Vec3 point, Vec3 vert1, Vec3 vert2, Vec3 vert3, Vec3 e1, Vec3 e2) {
  // See if point on a plane is within a triangle
  // Vec3 ep = point.minus(vert1);
  // float d11 = dot(e1,e1);
  // float d12 = dot(e1,e2);
  // float d22 = dot(e2,e2);
  // float dp1 = dot(ep,e1);
  // float dp2 = dot(ep,e2);
  // float D       = d11*d22 - d12*d12;
  // float D_beta  = d22*dp1 - d12*dp2;
  // float D_gamma = d11*dp2 - d12*dp1;
  // float beta = D_beta/D;
  // float gamma = D_gamma/D;
  // float alpha = 1 - beta + gamma;
  // print(alpha);
  // print(" ");
  // print(beta);
  // print(" ");
  // print(gamma);
  // println();
  // return (alpha > 0.0000001 && alpha < 1.0000001);

  // Source inspired by Scratchapixel:
  // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates

  Vec3 surfaceNormal = cross(e1, e2);
  Vec3 C;

  Vec3 edge0 = vert2.minus(vert1);
  Vec3 vp0 = point.minus(vert1);
  C = cross(edge0, vp0);
  if (dot(surfaceNormal, C) < 0) { return false; }

  Vec3 edge1 = vert3.minus(vert2);
  Vec3 vp1 = point.minus(vert2);
  C = cross(edge1, vp1);
  if (dot(surfaceNormal, C) < 0) { return false; }

  Vec3 edge2 = vert1.minus(vert3);
  Vec3 vp2 = point.minus(vert3);
  C = cross(edge2, vp2);
  if (dot(surfaceNormal, C) < 0) { return false; }

  return true;
}
// Created for CSCI 5611 by Liam Tyler

class Camera
{
  Camera()
  {
    position      = new PVector( 0, 0, 0 ); // initial position
    theta         = 0; // rotation around Y axis. Starts with forward direction as ( 0, 0, -1 )
    phi           = 0; // rotation around X axis. Starts with up direction as ( 0, 1, 0 )
    moveSpeed     = 50;
    turnSpeed     = 1.57f; // radians/sec
    
    // dont need to change these
    negativeMovement = new PVector( 0, 0, 0 );
    positiveMovement = new PVector( 0, 0, 0 );
    negativeTurn     = new PVector( 0, 0 ); // .x for theta, .y for phi
    positiveTurn     = new PVector( 0, 0 );
    fovy             = PI / 4;
    aspectRatio      = width / (float) height;
    nearPlane        = 0.1f;
    farPlane         = 10000;
  }
  
  public void Update(float dt)
  {
    theta += turnSpeed * ( negativeTurn.x + positiveTurn.x)*dt;
    
    // cap the rotation about the X axis to be less than 90 degrees to avoid gimble lock
    float maxAngleInRadians = 85 * PI / 180;
    phi = min( maxAngleInRadians, max( -maxAngleInRadians, phi + turnSpeed * ( negativeTurn.y + positiveTurn.y ) * dt ) );
    
    // re-orienting the angles to match the wikipedia formulas: https://en.wikipedia.org/wiki/Spherical_coordinate_system
    // except that their theta and phi are named opposite
    float t = theta + PI / 2;
    float p = phi + PI / 2;
    PVector forwardDir = new PVector( sin( p ) * cos( t ),   cos( p ),   -sin( p ) * sin ( t ) );
    PVector upDir      = new PVector( sin( phi ) * cos( t ), cos( phi ), -sin( t ) * sin( phi ) );
    PVector rightDir   = new PVector( cos( theta ), 0, -sin( theta ) );
    PVector velocity   = new PVector( negativeMovement.x + positiveMovement.x, negativeMovement.y + positiveMovement.y, negativeMovement.z + positiveMovement.z );
    position.add( PVector.mult( forwardDir, moveSpeed * velocity.z * dt ) );
    position.add( PVector.mult( upDir,      moveSpeed * velocity.y * dt ) );
    position.add( PVector.mult( rightDir,   moveSpeed * velocity.x * dt ) );
    
    aspectRatio = width / (float) height;
    perspective( fovy, aspectRatio, nearPlane, farPlane );
    camera( position.x, position.y, position.z,
            position.x + forwardDir.x, position.y + forwardDir.y, position.z + forwardDir.z,
            upDir.x, upDir.y, upDir.z );
  }
  
  // only need to change if you want difrent keys for the controls
  public void HandleKeyPressed()
  {
    if ( key == 'w' ) positiveMovement.z = 1;
    if ( key == 's' ) negativeMovement.z = -1;
    if ( key == 'a' ) negativeMovement.x = -1;
    if ( key == 'd' ) positiveMovement.x = 1;
    if ( key == 'q' ) positiveMovement.y = 1;
    if ( key == 'e' ) negativeMovement.y = -1;
    
    if ( keyCode == LEFT )  negativeTurn.x = 1;
    if ( keyCode == RIGHT ) positiveTurn.x = -1;
    if ( keyCode == UP )    positiveTurn.y = 1;
    if ( keyCode == DOWN )  negativeTurn.y = -1;
  }
  
  // only need to change if you want difrent keys for the controls
  public void HandleKeyReleased()
  {
    if ( key == 'w' ) positiveMovement.z = 0;
    if ( key == 'q' ) positiveMovement.y = 0;
    if ( key == 'd' ) positiveMovement.x = 0;
    if ( key == 'a' ) negativeMovement.x = 0;
    if ( key == 's' ) negativeMovement.z = 0;
    if ( key == 'e' ) negativeMovement.y = 0;
    
    if ( keyCode == LEFT  ) negativeTurn.x = 0;
    if ( keyCode == RIGHT ) positiveTurn.x = 0;
    if ( keyCode == UP    ) positiveTurn.y = 0;
    if ( keyCode == DOWN  ) negativeTurn.y = 0;
  }
  
  // only necessary to change if you want different start position, orientation, or speeds
  PVector position;
  float theta;
  float phi;
  float moveSpeed;
  float turnSpeed;
  
  // probably don't need / want to change any of the below variables
  float fovy;
  float aspectRatio;
  float nearPlane;
  float farPlane;  
  PVector negativeMovement;
  PVector positiveMovement;
  PVector negativeTurn;
  PVector positiveTurn;
};
  public void settings() {  size(640,480,OPENGL); }
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "ParticleSimulator" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
