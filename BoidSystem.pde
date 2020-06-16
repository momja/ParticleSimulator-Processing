public class BoidSystem {
    private ArrayList<Vec3> boidCoords;
    // private ArrayList<Quaternion> boidOrientations;
    private ArrayList<Vec3> boidVelocities;

    int boidCount = 10;
    Vec3 boidCenterOfMass;
    float boidSize = 10;
    PImage boidTexture = null;
    PShape boidModel = null;

    public BoidSystem() {
        boidCoords = new ArrayList<Vec3>();
        // boidOrientations = new ArrayList<Quaternion>();
        boidVelocities = new ArrayList<Vec3>();
        initializePositions();
    }
    
    public void initializePositions() {
        for (int i = 0; i < boidCount; i++) {

        }
    }

    public void updateBoidPositions(float dt) {
        Vec3 v1, v2, v3;
        for (int i = 0; i < boidCount; i++) {
            v1 = flyTowardsCenter(i);
            v2 = keepDistance(i);
            v3 = matchVelocity(i);

            boidVelocities.get(i).add(v1.plus(v2).plus(v3));
            boidCoords.get(i).add(boidVelocities.get(i));
        }

        boidCenterOfMass = new Vec3(0,0,0);
        for (int i = 0; i < boidCount; i++) {
            boidCenterOfMass.add(boidCoords.get(i));
        }
        boidCenterOfMass = boidCenterOfMass.times(1.f/boidCount);
    }

    public void drawBoids() {
        if (boidTexture != null) {
            for (int i = 0; i < boidCount; i++) {
                // Draw texture
                Vec3 pos = boidCoords.get(i);
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
            for (int i = 0; i < boidCount; i++) {
                // Draw boid model
                Vec3 pos = boidCoords.get(i);
                push();
                translate(pos.x, pos.y, pos.z);
                shape(boidModel);
                pop();
            }
        }
        else {
            for (int i = 0; i < boidCount; i++) {
                // Draw a point
                Vec3 pos = boidCoords.get(i);
                push();
                stroke(boidSize);
                point(pos.x, pos.y, pos.z);
                pop();
            }
        }
    }

    private Vec3 flyTowardsCenter(int idx) {
        // given an index of a certain boid,
        // find the center of mass of nearby boids
        // and return a weight in that direction
        
        // Fly towards the center of mass excluding the current boid
        Vec3 perceivedCenterOfMass = boidCenterOfMass.times(boidCount).minus(boidCoords.get(idx)).times(1/(boidCount - 1));
        return perceivedCenterOfMass.minus(boidCoords.get(idx)).times(0.01);
    }

    private Vec3 keepDistance(int idx) {
        // look at all the nearby boids and objects
        // and return a weight that maneuvers away from them
        Vec3 c = new Vec3(0,0,0);
        Vec3 curBoidPosition = boidCoords.get(idx);
        for (int i = 0; i < boidCount; i++) {
            if (i == idx) {
                continue;
            }
            Vec3 otherBoidPosition = boidCoords.get(i);
            Vec3 betweenVec = otherBoidPosition.minus(curBoidPosition);
            if (betweenVec.length() < 100) {
                c.subtract(betweenVec);
            }
        }
        return c;
    }

    private Vec3 matchVelocity(int idx) {
        // find the velocities of nearby boids and approximate it
        Vec3 curBoidVelocty = boidVelocities.get(idx);
        Vec3 newPosition = new Vec3(0,0,0);;
        for (int i = 0; i < boidCount; i++) {
            Vec3 otherBoidVelocity = boidVelocities.get(i);
            if (i == idx) {
                continue;
            }
            newPosition.add(otherBoidVelocity);
        }
        newPosition.mul(boidCount - 1);
        newPosition.subtract(curBoidVelocty);
        newPosition.mul(1f/8f);
        return newPosition;
    }
}