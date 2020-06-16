public class CollisionTrigger {
    boolean isActive = false;

    public void onCollision(Vec3 point) {
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
        emitter.particleLifespanMax = 0.5;
        emitter.particleLifespanMin = 0.4;
        emitter.birthRate = 80;
        emitter.emitterLifespan = 0.06;
        emitter.r = 1.3;
        emitter.particleSpeed = 60;
        emitter.particleDirection = new Vec3(0,1,0);
        emitter.particleDirectionRange = 0.3;
        emitter.particleAcceleration = new Vec3(0,-9.8,0);
    }

    @Override
    public SpawnEmitter copy() {
        return new SpawnEmitter();
    }

    @Override
    public void onCollision(Vec3 point) {
        // TODO : Spawn a new emitter at the point of collision
        super.onCollision(point);
        emitter.emitterPosition = point;
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
    public void onCollision(Vec3 point) {
        // TODO : Spawn a textured quad that animates through a series
        // of raindrop splash images
        super.onCollision(point);
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