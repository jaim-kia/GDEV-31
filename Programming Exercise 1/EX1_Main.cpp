#include <iostream>
#include <cmath>

/**
 * Struct containing data for a 3D vector
 */
struct Vector3
{
public:
	/**
	 * X component
	 */
	float x;
	
	/**
	 * Y component
	 */
	float y;
	
	/**
	 * Z component
	 */
	float z;
	
	/**
	 * Constructor
	 * @param newX - New x component
	 * @param newY - New y component
	 * @param newZ - New z component
	 */
	Vector3(float newX, float newY, float newZ)
	{
		x = newX;
		y = newY;
		z = newZ;
	}
	
	/**
	 * Constructor
	 */
	Vector3()
		: Vector3(0.0f, 0.0f, 0.0f)
	{
	}
	
	/**
	 * Task 1: negative() member function
	 */
	Vector3 operator-() const
    {
        return Vector3(-x, -y, -z);
    }
	Vector3 negative()
	{
		Vector3 ret = -(*this);

		if (ret.x == 0) ret.x = 0;
		if (ret.y == 0) ret.y = 0;
		if (ret.z == 0) ret.z = 0;

		return ret;
	}
	
	/**
	 * Task 2: magnitude() member function
	 */
	float magnitude() {
		float mag = std::sqrt(x * x + y * y + z * z);
		return mag;
	}
	
	/**
	 * Task 3: squaredMagnitude() member function
	 */
	float squaredMagnitude() {
		float squaredMag = x * x + y * y + z * z;
		return squaredMag;
	}

	/**
	 * Task 4: normalized() member function
	 */
	Vector3 normalized() {
		float mag = magnitude();
		Vector3 norm;
		if (mag > 0) {
			norm = Vector3(x/mag, y/mag, z/mag);
		} else {
			norm = Vector3(0, 0, 0);
		}
		return norm;
	}
	
	/**
	 * Task 5: add() static member function
	 */
	Vector3 operator+(const Vector3& other) const
    {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }
	static Vector3 add(Vector3 a, Vector3 b)
	{
		Vector3 sum = a + b;
		return sum;
	}
	
	/**
	 * Task 6: subtract() static member function
	 */
	Vector3 operator-(const Vector3& other) const
    {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }
	static Vector3 subtract(Vector3 a, Vector3 b)
	{
		Vector3 difference = a - b;
		return difference;
	}

	/**
	 * Task 7: multiply() static member function
	 */
	static Vector3 multiply(Vector3 vec, float scalarVal)
	{
		Vector3 result = Vector3(vec.x * scalarVal, vec.y * scalarVal, vec.z * scalarVal);
		return result;
	}
	
	/**
	 * Task 8: dot() static member function
	 */
	static float dot(Vector3 a, Vector3 b)
	{
		float dotProd = (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
		return dotProd;
	}
	
	/**
	 * Task 9: cross() static member function
	 */
	static Vector3 cross(Vector3 a, Vector3 b)
	{
		Vector3 crossProd = Vector3((a.y * b.z)-(a.z * b.y), (a.z * b.x)-(a.x * b.z), (a.x * b.y)-(a.y * b.x));
		return crossProd;
	}
	
	/**
	 * Task 10: project() static member function
	 */
	static Vector3 project(Vector3 a, Vector3 b)
	{
		float dotProd = dot(a, b);
		float bSqMag = b.squaredMagnitude();
		
		Vector3 proj = Vector3((dotProd/bSqMag) * b.x, (dotProd/bSqMag) * b.y, (dotProd/bSqMag) * b.z);
		return proj;
	}

	/**
	 * Task 11: reflect() static member function
	 */
	static Vector3 reflect(Vector3 a, Vector3 b)
	{
		Vector3 ref = a - multiply(project(a, b), 2);
		return ref;
	}
	
	
	friend std::ostream& operator<<(std::ostream& os, const Vector3& vec3)
	{
		os << "(" << vec3.x << ", " << vec3.y << ", " << vec3.z << ")";
		return os;
	}
};



int main()
{
	/**
	 * This line is to test whether you can compile and run this cpp file.
	 * Remove this line when submitting.
	 */
	// std::cout << "Hello World!" << std::endl;
	
	/**
	 * Task 12: Read input from standard input, and output
	 * to standard output based on the sample output.
	 */
	int scalarVal;
	Vector3 A, B;

	std::cin >> A.x >> A.y >> A.z >> B.x >> B.y >> B.z >> scalarVal;
	std::cout << "A = " << A << "\n";
	std::cout << "B = " << B << "\n";
	std::cout << "S = " << scalarVal << "\n";
	std::cout << "-A = " << A.negative() << "\n";
	std::cout << "-B = " << B.negative() << "\n";
	std::cout << "Squared magnitude of A = " << A.squaredMagnitude() << "\n";
	std::cout << "Magnitude of A = " << A.magnitude() << "\n";
	std::cout << "Squared magnitude of B = " << B.squaredMagnitude() << "\n";
	std::cout << "Magnitude of B = " << B.magnitude() << "\n";
	std::cout << "A normalized = " << A.normalized() << "\n";
	std::cout << "B normalized = " << B.normalized() << "\n";
	std::cout << "A + B = " << Vector3::add(A, B) << "\n";
	std::cout << "A - B = " << Vector3::subtract(A, B) << "\n";
	std::cout << "A * S = " << Vector3::multiply(A, scalarVal) << "\n";
	std::cout << "A dot B = " << Vector3::dot(A, B) << "\n";
	std::cout << "A cross B = " << Vector3::cross(A, B) << "\n";
	std::cout << "Projection of A onto B = " << Vector3::project(A, B) << "\n";
	std::cout << "Reflection of A along B = " << Vector3::reflect(A, B) << "\n";
	return 0;
}

