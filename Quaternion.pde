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