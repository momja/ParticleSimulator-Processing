public class BoidSystem {
    private ArrayList<Vec3> boidCoords;
    private ArrayList<Quaternion> boidOrientations;
    private ArrayList<Vec3> boidVelocities;

    int boidCount = 0;
    Vec3 boidCenterOfMass;

    public BoidSystem() {
        boidCoords = new ArrayList<Vec3>();
        boidOrientations = new ArrayList<Quaternion>();
        boidVelocities = new ArrayList<Vec3>();
    }
    
    public void initializePositions() {
        
    }

    public void updateBoidPositions(float dt) {
        Vec3 v1, v2, v3;
        for (int i = 0; i < boidCount; i++) {
            v1 = flyTowardsCenter(i);
            v2 = keepDistance(i);
            v3 = matchVelocity(i);

            boidVelocities.get(i).add(v1.plus(v2).plus(v3));
            boidPositions.get(i).add(boidVelocities.get(i));
        }

        boidCenterOfMass = new Vec3(0,0,0);
        for (int i = 0; i < boidCountl i++) {
            boidCenterOfMass.add(boidPositions.get(i));
        }
        boidCenterOfMass /= boidCount;
    }

    public void drawBoids() {
        for (int i = 0; i < boidCount; i++) {
            
        }
    }

    private Vec3 flyTowardsCenter(int idx) {
        // given an index of a certain boid,
        // find the center of mass of nearby boids
        // and return a weight in that direction
        
        // Fly towards the center of mass excluding the current boid
        Vec3 perceivedCenterOfMass = (boidCenterOfMass*boidCount - boidPositions.get(i))/(boidCount - 1);
        return perceivedCenterOfMass.minus(boidPositions.get(i)).times(0.01);
    }

    private Vec3 keepDistance(int idx) {
        // look at all the nearby boids and objects
        // and return a weight that maneuvers away from them
        Vec3 c = new Vec3(0,0,0);
        Vec3 curBoidPosition = boidPositions.get(idx);
        for (int i = 0; i < boidCount; i++) {
            if (i == idx) {
                continue;
            }
            Vec3 otherBoidPosition = boidPositions.get(i);
            Vec3 betweenVec = otherBoidPosition.minus(curBoidPosition);
            if (betweenVec.length() < 100) {
                c.substract(betweenVec);
            }
        }
        return c
    }

    private Vec3 matchVelocity(int idx) {
        // find the velocities of nearby boids and approximate it
        Vec3 curBoidVelocty = boidVelocities.get(idx);
        Vec3 newPosition;
        for (int i = 0; i < boidCount; i++) {
            Vec3 otherBoidVelocity = boidVelocities.get(i);
            if (i == idx) {
                continue;
            }
            newPosition.add(otherBoidVelocity)
        }
        newPosition.mul(boidCount - 1);
        newPosition.subtract(curBoidVelocty);
        newPosition.mul(1f/8f);
        return newPosition
    }
}